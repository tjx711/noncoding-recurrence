library(MASS)
library(boot)
library(bootstrap)
library(DAAG)
library(caret)

setwd("/home/tanjinxu/Project/noncoding")
#setwd("/Users/tanjinxu/Documents/project/noncoding")
source("source/noncoding_loaddata.R")

## Load the data
## the last column data is only used for sampling
## not used for training or testing
load(file="data/train_test.sampled.mixed.331.0.3.rda")
train_data <- datlist.sample$train
test_data <- datlist.sample$test
train_data <- train_data[, 1:(ncol(train_data)-1)]
test_data <- test_data[, 1:(ncol(test_data)-1)]

## Here we substract the rawcounts number by 1 to satisfy
## the Negative-Binomial Distribution assumption, which
## we verfied by the rawcounts distribution plots
train_data$rawcounts <- train_data$rawcounts - 1
test_data$rawcounts <- test_data$rawcounts - 1
train_x <- train_data[,1:(ncol(train_data)-2)]
train_y <- train_data$rawcounts
test_x <- test_data[,1:(ncol(test_data)-2)]
test_y <- test_data$rawcounts

train_data <- cbind(train_x, train_y)
test_data <- cbind(test_x, test_y)
colnames(train_data)[ncol(train_data)] <- "rawcounts"
colnames(test_data)[ncol(test_data)] <- "rawcounts"

##
## TO-DO: try what we get if the data is normalized/standarized
##
pre_normalize <- FALSE
if (pre_normalize){
  train_x_normed <- data.frame(sapply(train_x,function(x) (x-mean(x))/sd(x)))
  test_x_normed <- data.frame(sapply(test_x,function(x) (x-mean(x))/sd(x)))
  train_data <- cbind(train_x_normed, train_y)
  test_data <- cbind(test_x_normed, test_y)
  colnames(train_data)[ncol(train_data)] <- "rawcounts"
  colnames(test_data)[ncol(test_data)] <- "rawcounts"
}
  
## Apply Poisson-Regression with 10-fold Cross-validation
kfold <- 3
niters <- 1

res.err.score.df <- data.frame()
res.fit.imp <- data.frame(rep(0,29))

for(i in 1:niters){
  ## Apply Negative Binomial Regression model
  glm.fit.res <- glm.nb(rawcounts ~ ., data = train_data, link=log)
  
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
  
  ## Evaluation metrics used for all regression models in this project
  ## Reference: http://scikit-learn.org/stable/modules/model_evaluation.html
  ## 1> Mean Absolute Error (MUAE)
  ## 2> Mean Squared Error (MSE)
  ## 3> Median Absolute Error (MEAE)
  ## 4> Explained Variance Score (EAS)
  ## 5> R-squared Score (R2)
  error.train.muae <- mean(abs(train_data$rawcounts - glm.fit.res$fitted.values))
  error.train.mse <- mean((train_data$rawcounts - glm.fit.res$fitted.values)^2)
  error.train.meae <- median(abs(train_data$rawcounts - glm.fit.res$fitted.values))
  error.train.eas <- 1 - var(train_data$rawcounts - glm.fit.res$fitted.values) / var(train_data$rawcounts)
  error.train.r2 <- 1 - sum((train_data$rawcounts - glm.fit.res$fitted.values)^2) / 
                        sum((train_data$rawcounts - mean(train_data$rawcounts)^2))
 
  cat(sprintf("iter:%d, training: mean AE:%f, mean SE:%f, median AE:%f, exp var:%f, r2:%f\n",
              i , error.train.muae, error.train.mse, error.train.meae, error.train.eas, error.train.r2))
  
  ## Now fitting the testing data to the learned model
  test.x <- test_data[, 1:ncol(test_data)-1]
  test.y <- test_data[, ncol(test_data)]
  pred.rawcounts <- predict.glm(glm.fit.res, newdata=test.x, type="response", se.fit=FALSE)
  
  error.test.muae <- mean(abs(test.y - pred.rawcounts))
  error.test.mse <- mean((test.y - pred.rawcounts)^2)
  error.test.meae <- median(abs(test.y - pred.rawcounts))
  error.test.eas <- 1 - var(test.y - pred.rawcounts) / var(test.y)
  error.test.r2 <- 1 - sum((test.y - pred.rawcounts)^2) / sum((test.y - mean(test.y))^2)
  
  cat(sprintf("iter:%d, testing: mean AE:%f, mean SE:%f, median AE:%f, exp var:%f, r2:%f\n",
              i, error.test.muae, error.test.mse, error.test.meae, error.test.eas, error.test.r2))

  res <- c(error.train.muae, error.train.mse, error.train.meae, error.train.eas, error.train.r2,
          error.test.muae, error.test.mse, error.test.meae, error.test.eas, error.test.r2)
  
  res.err.score.df <- rbind(res.err.score.df,res)
  res.fit.imp <- cbind(res.fit.imp, fit.imp)
}##//end for

## Write the result to file so we can plot with output from skilearn regression models
write.table(res.err.score.df,"nc.poisson.res.err.score.tsv",sep="\t",row.names=TRUE)
write.table(rowMeans(res.fit.imp), "nc.poisson.res.imp.tsv",sep="\t")
cat("average errors:\n",colMeans(res.err.score.df))
cat("\n")

