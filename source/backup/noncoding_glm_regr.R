#library(gtools)
library(MASS)
library(boot)
#library(ResourceSelection)
#library(car)

setwd("/home/tanjinxu/Project/noncoding")
source("noncoding_loaddata.R")

## Load the data
load(file="train_test.rda")
train_data <- datlist$train
test_data <- datlist$test

## TO-DO: try what we get if the data is normalized/standarized

## Apply Negative-Binomial Regression
#myres <- glm.nb(counts ~ ., data = mydata, control=glm.control(trace=10,maxit=100))
#print(summary(myres))

## Apply Poisson-Regression
## TODO: How to apply cross-validation in GLM model calls?
## Split the dataset into training set and validation set
## Using K-fold cross_validation
kfold <- 2
#iterations <- c(25,50,100,150,200,300)
iterations <- c(25,30)
n_iters <- length(iterations)

avg.train.dev <- rep(0,n_iters)
avg.train.aic <- rep(0,n_iters)
avg.valid.mse <- rep(0,n_iters)

for(i in 1:n_iters){
	mse <- rep(0,kfold)
	pois.fits <- c()

	print("start training...")
  cat(sprintf("start training with %d iters...",iterations[i]))
	for(j in 1:kfold){
		dataset <- splitData(train_data,1-(1/kfold))
		train_subset <- dataset$train
		validation_set <- dataset$test
		X.validate <- validation_set[,1:(ncol(validation_set)-1)]
		y.validate <- validation_set[,ncol(validation_set)]

		## poisson-regression
		fm_pois <- glm(counts ~ .,data=train_subset,family=poisson(link=log),maxit=iterations[i])
		pois.fits <- c(pois.fits,fm_pois)
		if (!fm_pois$converged){
			print("model does not converge in iteration:",iter)
		}
		
		## now apply the model with the validation set and calculate the metrics
		valid.pred.counts <- predict.glm(fm_pois,newdata=X.validate,type="response",se.fit=FALSE)
		mse[j] <- mean((valid.pred.counts - y.validate)^2)
	}

	avg.valid.mse[i] <- mean(mse)
	avg.train.aic[i] <- mean(pois.fits$aic)
	avg.train.dev[i] <- mean(pois.fits$deviance)
}

## Plot the model evaluation metrics
plot(avg.train.dev)
par(new=T)
plot(avg.train.aic)
par(new=T)
plot(avg.valid.mse)
par(new=F)
 
## Get best param (maxiter) for min(MSE),min(deviance),max(AIC)
cat(sprintf("min_MSE: %f, min_dev: %f, max_aic: %f\n", min(avg.valid.mse), min(avg.train.dev), min(avg.train.aic)))

best.iters <- c(which.min(avg.valid.mse), which.min(avg.valid.mse),which.min(avg.train.aic))
cat(sprintf("iters for min MSE/DEV/AIC: %d\n",best.iters))

## Use the parameter for minimum AIC and retrain the model
iter_pred <- best.iters[3]
best.pois.fit <- glm(counts ~ .,data=train_data,family=poisson(link=log),maxit=iter_pred)
print(summary(best.pois.fit))


## Get the predicted counts for each observation in testing set
X.test <- test_data[,1:(ncol(test_data)-1)]
y.test <- test_data[,ncol(test_data)]
test.pred.counts <- predict.glm(best.pois.fit,newdata=X.test,type="response",se.fit=FALSE)
test.mse <- mean((test.pred.counts - y.test)^2)
cof.spearman <- cor(test.pred.counts,y.test,method="spearman")
cof.pearson <- cor(test.pred.counts,y.test,method="pearson")
cat(sprintf("Result: mse=%f,spearman=%f,pearson=%f\n",test.mse,cof.spearman,cof.pearson))

## Evaluate the goodness_of_fit of the trained model
## TODO: What value could be simply used as evaluation of the goodness_of_fit?
#anova(fm_pois,test="Chisq")
#fm_pois$deviance / fm_pois$df.residual
#fm_pois$null.deviance / fm_pois$df.null
#1 - pchisq(fm_pois$deviance,fm_pois$df.residual)
#pchisq(fm_pois$null.deviance - fm_pois$deviance, fm_pois$df.null - fm_pois$df.residual, lower.tail=FALSE)
#hl <- hoslem.test(train_data$counts, fitted(fm_pois))
#cbind(hl$observed,hl$expected)
#pseudo.r.square <- 100*(fm_pois$null.deviance - fm_pois$deviance)/fm_pois$null.deviance

#fm_diag <- glm.diag(fm_pois)
#glm.diag.plots(fm_pois,fm_diag)
#printdat <- data.frame(mydata,pred=fm_pois$fitted.values)
#counts_hat <- predict(fm_pois, type="response")

