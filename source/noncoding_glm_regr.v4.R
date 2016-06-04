library(MASS)
library(boot)
library(bootstrap)
library(DAAG)
library(caret)
require(VGAM)

setwd("/home/tanjinxu/Project/noncoding")
#setwd("/Users/tanjinxu/Documents/project/noncoding")
source("source/noncoding_loaddata.R")

glmRegression <- function(train_set, test_set, model, outstr){
  nrows.train <- nrow(train_set)
  ncols.train <- ncol(train_set)
  nrows.test <- nrow(test_set)
  ncols.test <- ncol(test_set)
  colnames(train_set)[ncols.train] <- "counts"
  
  if (model == "nb"){
    ## apply negative binomial regression
    glm.fit.res <- glm.nb(counts~., data=train_set, link="log")
  }
  else if (model == "poisson"){
    ## apply poisson regression
    glm.fit.res <- glm(counts ~., data = train_set, family=poisson(link="log"))
  }
  else if (model == "quasi"){
    ## apply quasi-poisson regression
    glm.fit.res <- glm(counts ~., data=train_set, family=quasipoisson(link="log"))
  }
  
  ## Check if the model coverged or not
  cat(sprintf("Regression model %s converged: %d\n", model,glm.fit.res$converged))
  
  ## Call varImp from "caret" package to calculate the importance score
  ## so that we can compare this score with the scores from the skilearn
  ## regression models
  fit.imp <- varImp(glm.fit.res, scale=FALSE)
  fit.imp$Overall <- (fit.imp$Overall - min(fit.imp$Overall))/
                     (max(fit.imp$Overall)-min(fit.imp$Overall))
  
  ## Evaluation metrics used for all regression models in this project
  ## Reference: http://scikit-learn.org/stable/modules/model_evaluation.html
  ## 1> Mean Absolute Error (MUAE)
  ## 2> Mean Squared Error (MSE)
  ## 3> Median Absolute Error (MEAE)
  ## 4> Explained Variance Score (EAS)
  ## 5> R-squared Score (R2)
  err.train.muae <- mean(abs(train_set$counts - glm.fit.res$fitted.values))
  err.train.mse  <- mean((train_set$counts - glm.fit.res$fitted.values)^2)
  err.train.meae <- median(abs(train_set$counts - glm.fit.res$fitted.values))
  err.train.eas  <- 1 - var(train_set$counts - glm.fit.res$fitted.values) / var(train_set$counts)
  err.train.r2   <- 1 - (sum((train_set$counts - glm.fit.res$fitted.values)^2) / 
                        sum((train_set$counts - mean(train_set$counts))^2))

  log.counts <- 0.0
  log.fitted <- 0.0
  if (model == "nb"){
    log.counts <- log(train_set$counts + 1)
    log.fitted <- log(glm.fit.res$fitted.values + 1)
  }
  else{
    log.counts <- log(train_set$counts)
    log.fitted <- log(glm.fit.res$fitted.values)
  }
  
  log.err.train.muae <- mean(abs(log.counts - log.fitted))
  log.err.train.mse <- mean((log.counts - log.fitted)^2)
  log.err.train.meae <- median(abs(log.counts - log.fitted))
  log.err.train.eas <- 1 - (var(log.counts - log.fitted) / var(log.counts))
  log.err.train.r2 <- 1 - (sum((log.counts - log.fitted)^2) / 
                           sum((log.counts - mean(log.counts))^2))
  
  cat(sprintf("model:%s,training: MUAE:%f, MUSE:%f, MDAE:%f, \
              expvar:%f, r2:%f\n",
              model,err.train.muae, err.train.mse, err.train.meae, 
              err.train.eas, err.train.r2))
  
  #cat(sprintf("model:%s,training: log.mean AE:%f, log.mean SE:%f, log.median AE:%f, \
  #            log.exp var:%f, log.r2:%f\n",
  #            model,log.err.train.muae, log.err.train.mse, log.err.train.meae, 
  #            log.err.train.eas, log.err.train.r2))
  
  ## Now fitting the testing data to the learned model
  test.x <- test_set[, 1:ncols.test-1]
  test.y <- test_set[, ncols.test]
  pred.counts <- predict.glm(glm.fit.res, newdata=test.x, type="response")
  
  ## Check the ranking
  test.y.rank <- order(test.y,decreasing = TRUE)
  test.y.rank.top100.idx <- which(test.y.rank <= 100)
  test.y.rank.top100 <- test.y.rank[test.y.rank.top100.idx]
  test.yhat.rank.top100 <- (order(pred.counts, decreasing = TRUE))[test.y.rank.top100.idx]
  
  err.test.muae <- mean(abs(test.y - pred.counts))
  err.test.mse <- mean((test.y - pred.counts)^2)
  err.test.meae <- median(abs(test.y - pred.counts))
  err.test.eas <- 1 - (var(test.y - pred.counts) / var(test.y))
  err.test.r2 <- 1 - (sum((test.y - pred.counts)^2) / sum((test.y - mean(test.y))^2))
  err.test.rank <- mean(abs(test.y.rank.top100 - test.yhat.rank.top100))
  
  shuff_times <- 10
  err.test.mse.shuffled <- 0
  err.test.muae.shuffled <- 0
  err.test.meae.shuffled <- 0
  err.test.eas.shuffled <- 0
  err.test.r2.shuffled <- 0
  sum.err.mse.shuffled <- 0
  sum.err.muae.shuffled <- 0
  sum.err.meae.shuffled <- 0
  sum.err.eas.shuffled <- 0
  sum.err.r2.shuffled <- 0
  
  for (i in 1:shuff_times){
    ## Shuffle the features value of testing set and evaluate the random error
    test.x.shuffled <- test.x
    for (col in 1:ncol(test.x)){
      test.x.shuffled[,col] <- sample(test.x[,col])
    }
    
    pred.counts.shuffled <- predict.glm(glm.fit.res, newdata=test.x.shuffled, type="response")
    err.test.mse.shuffled <- mean((test.y - pred.counts.shuffled)^2)
    err.test.muae.shuffled <- mean(abs(test.y - pred.counts.shuffled))
    err.test.meae.shuffled <- median(abs(test.y - pred.counts.shuffled))
    err.test.eas.shuffled <- 1 - (var(test.y - pred.counts.shuffled) / var(test.y))
    err.test.r2.shuffled <- 1- (sum((test.y - pred.counts.shuffled)^2) / 
                                sum((test.y - mean(test.y))^2))
  
    cat(sprintf("model:%s, shuffled testing: mean AE:%f, mean SE:%f, median AE:%f, 
                exp var:%f, r2:%f\n",
                model, err.test.muae.shuffled, err.test.mse.shuffled, err.test.meae.shuffled,
                err.test.eas.shuffled, err.test.r2.shuffled))
    
    sum.err.mse.shuffled <- sum.err.mse.shuffled + err.test.mse.shuffled
    sum.err.muae.shuffled <- sum.err.muae.shuffled + err.test.muae.shuffled
    sum.err.meae.shuffled <- sum.err.meae.shuffled + err.test.meae.shuffled
    sum.err.eas.shuffled <- sum.err.eas.shuffled + err.test.eas.shuffled
    sum.err.r2.shuffled <- sum.err.r2.shuffled + err.test.r2.shuffled
  }
  err.test.mse.shuffled <- sum.err.mse.shuffled / shuff_times
  err.test.muae.shuffled <- sum.err.muae.shuffled / shuff_times
  err.test.meae.shuffled <- sum.err.meae.shuffled / shuff_times
  err.test.eas.shuffled <- sum.err.eas.shuffled / shuff_times
  err.test.r2.shuffled <- sum.err.r2.shuffled / shuff_times
  
  ## Calculate the log-scale testing error
  log.y <- 0.0
  log.preds <- 0.0
  if (model == "nb"){
    log.y <- log(test.y + 1)
    log.preds <- log(pred.counts + 1)
  }
  else{
    log.y <- log(test.y)
    log.preds <- log(pred.counts)
  }
  
  log.err.test.muae <- mean(abs(log.y - log.preds))
  log.err.test.mse <- mean((log.y - log.preds)^2)
  log.err.test.meae <- median(abs(log.y - log.preds))
  log.err.test.eas <- 1 - (var(log.y - log.preds) / var(log.y))
  log.err.test.r2 <- 1 - (sum((log.y - log.preds)^2) / sum((log.y - mean(log.y))^2))
  
  ## Print out the errors
  cat(sprintf("model:%s,testing: MUAE:%f, MUSE:%f, MDAE:%f, \ 
              EXPV:%f, r2:%f, RANK:%f\n",
              model,err.test.muae, err.test.mse, err.test.meae, 
              err.test.eas, err.test.r2, err.test.rank))
  #cat(sprintf("model:%s,testing: log.MUAE:%f, log.MUSE:%f, log.MDAE:%f, 
  #            log.EXPV:%f, log.r2:%f\n",
  #            model,log.err.test.muae, log.err.test.mse, log.err.test.meae, 
  #            log.err.test.eas, log.err.test.r2))
  cat(sprintf("model:%s, shuffled testing: MUAE:%f, MUSE:%f, MDAE:%f, \
             exp var:%f, r2:%f\n",
             model, err.test.muae.shuffled, err.test.mse.shuffled, err.test.meae.shuffled,
             err.test.eas.shuffled, err.test.r2.shuffled))
  
  res <- c(err.train.muae, err.train.mse, err.train.meae, 
           err.train.eas, err.train.r2,
           log.err.train.muae, log.err.train.mse, log.err.train.meae, 
           log.err.train.eas, log.err.train.r2,
           err.test.muae, err.test.mse, err.test.meae, 
           err.test.eas, err.test.r2, err.test.rank,
           log.err.test.muae, log.err.test.mse, log.err.test.meae, 
           log.err.test.eas, log.err.test.r2,
           err.test.muae.shuffled, err.test.mse.shuffled, err.test.meae.shuffled,
           err.test.eas.shuffled, err.test.r2.shuffled)
  res.names <- c("err.train.muae","err.train.mse","err.train.meae",
                 "err.train.eas","err.train.r2",
                 "log.err.train.muae","log.err.train.mse","log.err.train.meae",
                 "log.err.train.eas","log.err.train.r2",
                 "err.test.muae","err.test.mse","err.test.meae",
                 "err.test.eas","err.test.r2","err.test.rank",
                 "log.err.test.muae","log.err.test.mse","log.err.test.meae",
                 "log.err.test.eas","log.err.test.r2",
                 "err.test.muae.shuffled","err.test.mse.shuffled","err.test.meae.shuffled",
                 "err.test.eas.shuffled", "err.test.r2.shuffled")
  
  ## write to tsv file
  errfile <- paste(outstr,model,"err.tsv",sep=".")
  impfile <- paste(outstr,model,"imp.tsv",sep=".")
  write.table(cbind(res.names,res),errfile,sep="\t",row.names=FALSE,quote=FALSE)
  write.table(fit.imp, impfile, sep="\t")
  return(list(err=res,imp=fit.imp))
}

#-----------------------------------------------------------------
# Main function
#-----------------------------------------------------------------
## Load the data file into data frame
## !!! Note: for both training set and testing set
## 1) the last column is only used for re-sampling/split
## 2) the last 2nd and 3rd are response variables and others
##    are predictors/features

## The input whole noncoding mutation dataset, two categories:
## 1) the complete dataset without any filtering
## 2) the dataset with filtered features only
#in.data.file <- "data/for_R/train_test.filtered.420.rda"

## the data file format:
##  "mixed" -- two response variables rawcounts and samplefreq
##             which should be handled separately
##  "331"-- the date when the file was generated
##  "0.3"-- the sampling rate, 30% of the whole data is used 
#load(file="data/train_test.sampled.mixed.0.4.rda")
sampling.rates <- c(1.0)
#sampling.rates <- c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8)
iterations <- 1
resampling <- FALSE #//<= TRUE - sampling from scratch; FALSE - use existing sampled data
#models <- c("poisson","quasi","nb")
models <- c("nb","poisson")
## Uniformly sample data for regression due to huge volume of dataset
for (rate in sampling.rates){
  for (iter in 1:iterations){
    ## For each sampling rate, sample 10 times
    #sampled.data.file <- paste("data/for_R/train_test.sampled.filtered",rate, iter,"421.rda", sep=".") 
    sampled.data.file <- "data/train_test.filtered.trim.504.all.rda"
    
    ## Use fixed sampling rate of 0.6 for test data
    if (resampling){
      ## Here we fixed the sampling rate for test
      sampleData(rate, 0.6, in.data.file, sampled.data.file)
    }
    
    load(file=sampled.data.file)
    train_data <- datlist$train
    test_data <- datlist$test

    train.x <- train_data[,1:(ncol(train_data)-3)]
    train.y.raw <- train_data$rawcounts
    train.y.freq <- train_data$samplefreq
    test.x<- test_data[,1:(ncol(test_data)-3)]
    test.y.raw <- test_data$rawcounts
    test.y.freq <- test_data$samplefreq

    for (model in models){
      cat(sprintf("---start regression analysis with model:%s\n", model))
  
      if (model == "nb"){
        ## Here we substract the rawcounts number by 1 to satisfy
        ## the Negative-Binomial Distribution assumption, which
        ## will be verified by the rawcounts distribution plots
        train.y.raw <- train.y.raw - 1
        train.y.freq <- train.y.freq - 1
        test.y.raw <- test.y.raw - 1
        test.y.freq <- test.y.freq - 1
      }
  
      normalize <- FALSE
      if (normalize){
        train.mean <- apply(train.x, 2, mean)
        train.std <- apply(train.x, 2, sd)
        train.mean.mtx <- matrix(train.mean, nrow=nrow(train.x), ncol=length(train.mean), byrow=TRUE)
        train.std.mtx <- matrix(train.std, nrow=nrow(train.x), ncol=length(train.std), byrow=TRUE)
        train.x.normed <- (train.x - train.mean.mtx) / train.std.mtx
        train.x <- train.x.normed
    
        test.mean <- apply(test.x, 2, mean)
        test.std <- apply(test.x, 2, sd)
        test.mean.mtx <- matrix(test.mean, nrow=nrow(test.x), ncol=length(test.mean), byrow=TRUE)
        test.std.mtx <- matrix(test.std, nrow=nrow(test.x), ncol=length(test.std), byrow=TRUE)
        test.x.normed <- (test.x - test.mean.mtx) / test.std.mtx
        test.x <- test.x.normed
      }
  
      if (FALSE){
        ## Use raw counts as response variable
        reslist.raw <- glmRegression(cbind(train.x, train.y.raw),
                                     cbind(test.x, test.y.raw),
                                     model,
                                     paste("glm_out/nc.res.raw.filtered",rate,iter,sep="."))
      }
  
      ## Use counts with sample frequency as response variable
      reslist.freq <- glmRegression(cbind(train.x, train.y.freq),
                                    cbind(test.x, test.y.freq),
                                    model,
                                    paste("glm_out/nc.res.freq.filtered",rate,iter,sep="."))
    }#//end for
  }#//end for
}#//end for

