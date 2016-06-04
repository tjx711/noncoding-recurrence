library(MASS)
library(boot)
library(bootstrap)
library(DAAG)
library(caret)

setwd("/home/tanjinxu/Project/noncoding")
#setwd("/Users/tanjinxu/Documents/project/noncoding")
source("noncoding_loaddata.R")

## Load the data
load(file="train_test.rda")
train_data <- datlist$train
test_data <- datlist$test

## Sampling subset from the whole data
train_subset <- train_data[sample(nrow(train_data),0.3*nrow(train_data)),1:ncol(train_data)-1]
test_subset <- test_data[sample(nrow(test_data),0.3*nrow(test_data)),1:ncol(train_data)-1]
train_data <- train_subset
test_data <- test_subset

## Here we substract the counts number by 1 to satisfy
## the Negative-Binomial Distribution assumption, which
## we verfied by the counts distribution plots
train_data$counts <- train_data$counts - 1
test_data$counts <- test_data$counts - 1
train_x <- train_data[,1:ncol(train_data)-1]
train_y <- train_data$counts
test_x <- test_data[,1:ncol(test_data)-1]
test_y <- test_data$counts

##
## TO-DO: try what we get if the data is normalized/standarized
##
pre_normalize <- FALSE
if (pre_normalize){
  train_x_normed <- data.frame(sapply(train_x,function(x) (x-mean(x))/sd(x)))
  test_x_normed <- data.frame(sapply(test_x,function(x) (x-mean(x))/sd(x)))
  train_data <- cbind(train_x_normed, train_y)
  test_data <- cbind(test_x_normed, test_y)
  colnames(train_data)[ncol(train_data)] <- "counts"
  colnames(test_data)[ncol(test_data)] <- "counts"
}
  
## Apply Poisson-Regression with 10-fold Cross-validation
kfold <- 10
niters <- 1

res.err.score.df <- data.frame()
res.fit.imp <- data.frame(rep(0,29))

for(i in 1:niters){
  ## Apply Poisson Regression model
  #glm.fit.res <- glm(counts ~ ., data = train_data, family=poisson(link=log))
  
  ## Apply Negative Binomial Regression model
  glm.fit.res <- glm.nb(counts ~ ., data = train_data, link=log)
  
  ## Check if the model coverged or not
  cat("Regression model converged?", glm.fit.res$converged)
  cat("\n")

  ## Extract feature importances
  fit.sum <- summary(glm.fit.res)
  fit.cofs <- fit.sum$coefficients
  fit.zval <- fit.cofs[2:nrow(fit.cofs),4]
  
  ## Call varImp from "caret" package to calculate the importance score
  ## so that we can compare this score with the scores from the skilearn
  ## regression models
  fit.imp <- varImp(glm.fit.res, scale=FALSE)
  fit.imp$Overall <- (fit.imp$Overall - min(fit.imp$Overall))/
                      (max(fit.imp$Overall)-min(fit.imp$Overall))
  #fit.imp.adj <- fit.imp[order(-fit.imp$Overall),,drop=FALSE]

  ## Not sure why here the result returned by residulas is not equal to the real residual
  ## calculated by (y - y_hat)
  #error.train <- mean(residuals(glm.fit.res)^2)
  
  ## Calculate the training error and the cross-validation error
  #error.train <- mean((train_data$counts - glm.fit.res$fitted.values)^2)
  #error.valid <- cv.glm(train_data,glm.fit.res,K=kfold)$delta
  
  ## Check the goodness_of_fit of the learned model
  ## Calculate the spearman coefficient and the pearson coefficient
  #cor.sp <- cor(train_data$counts, glm.fit.res$fitted.values,method="spearman",use="complete")
  #cor.pearson <- cor(train_data$counts, glm.fit.res$fitted.values,method="pearson",use="complete")
  #pseudo.r2 <- 100*(glm.fit.res$null.deviance - glm.fit.res$deviance)/glm.fit.res$null.deviance
  #pchi <- pchisq(glm.fit.res$null.deviance - glm.fit.res$deviance, 
  #               glm.fit.res$df.null - glm.fit.res$df.residual, lower.tail=FALSE)
  
  ## Evaluation metrics used for all regression models in this project
  ## Reference: http://scikit-learn.org/stable/modules/model_evaluation.html
  ## 1> Mean Absolute Error (MUAE)
  ## 2> Mean Squared Error (MSE)
  ## 3> Median Absolute Error (MEAE)
  ## 4> Explained Variance Score (EAS)
  ## 5> R-squared Score (R2)
  error.train.muae <- mean(abs(train_data$counts - glm.fit.res$fitted.values))
  error.train.mse <- mean((train_data$counts - glm.fit.res$fitted.values)^2)
  error.train.meae <- median(abs(train_data$counts - glm.fit.res$fitted.values))
  error.train.eas <- 1 - var(train_data$counts - glm.fit.res$fitted.values) / var(train_data$counts)
  error.train.r2 <- 1 - sum((train_data$counts - glm.fit.res$fitted.values)^2) / 
                        sum((train_data$counts - mean(train_data$counts)^2))
  
  ## Print the train/validation metric values
  #cat(sprintf("iter:%d, train error:%f, spearman:%f, pearson:%f, validation error:%f, pseudo_r2:%f, pchi:%f\n", 
  #            i, error.train, cor.sp, cor.pearson, error.valid[1], pseudo.r2, pchi))
  cat(sprintf("iter:%d, training: mean AE:%f, mean SE:%f, median AE:%f, exp var:%f, r2:%f\n",
              i , error.train.muae, error.train.mse, error.train.meae, error.train.eas, error.train.r2))
  
  ## Now fitting the testing data to the learned model
  test.x <- test_data[,1:ncol(test_data)-1]
  test.y <- test_data[,ncol(test_data)]
  pred.counts <- predict.glm(glm.fit.res,newdata=test.x,type="response",se.fit=FALSE)
  #error.test <- mean((pred.counts - test.y)^2)
  #cor.sp.test <- cor(test.y, pred.counts, method="spearman", use="complete")
  #cor.pearson.test <- cor(test.y, pred.counts, method="pearson", use="complete")
  
  error.test.muae <- mean(abs(test.y - pred.counts))
  error.test.mse <- mean((test.y - pred.counts)^2)
  error.test.meae <- median(abs(test.y - pred.counts))
  error.test.eas <- 1 - var(test.y - pred.counts) / var(test.y)
  error.test.r2 <- 1 - sum((test.y - pred.counts)^2) / sum((test.y - mean(test.y))^2)
  
  ## Print the testing metric values
  #cat(sprintf("iter:%d, test error:%f, spearman:%f, pearson:%f\n",
  #            i, error.test, cor.sp.test, cor.pearson.test))
  cat(sprintf("iter:%d, testing: mean AE:%f, mean SE:%f, median AE:%f, exp var:%f, r2:%f\n",
              i, error.test.muae, error.test.mse, error.test.meae, error.test.eas, error.test.r2))

  res <- c(error.train.muae, error.train.mse, error.train.meae, error.train.eas, error.train.r2,
          error.test.muae, error.test.mse, error.test.meae, error.test.eas, error.test.r2)
  
  res.err.score.df <- rbind(res.err.score.df,res)
  res.fit.imp <- cbind(res.fit.imp, fit.imp)
}

## Write the result to file so we can plot with output from skilearn regression models
write.table(res.err.score.df,"nc.poisson.res.err.score.tsv",sep="\t",row.names=TRUE)
write.table(rowMeans(res.fit.imp), "nc.poisson.res.imp.tsv",sep="\t")
cat("average erros:\n",colMeans(res.err.score.df))

