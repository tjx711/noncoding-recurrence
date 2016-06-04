
## number of cases to simulate
N <- 100000

## make an example feature matrix
feature1 <- rnorm(N, mean=-10, sd=0.05)
feature2 <- rbeta(N, 1.5, 1.5)
feature3 <- runif(N)
feature4 <- sample.int(5, size=N, replace=TRUE)
feature5 <- rgamma(N, shape=1, scale=0.8)
feature6 <- rnorm(N, mean=-1) + rnorm(N, mean=1)

all_features <- cbind(feature1,
                      feature2,
                      feature3,
                      feature4,
                      feature5,
                      feature6)
                      
coefficients <- c(1, -0.2, -1, -0.1, 1, 0.05)


hidden_covariate_multiplier <- matrix(rep(coefficients, N), ncol=6, byrow=TRUE)

hidden_covariate_real <- apply(hidden_covariate_multiplier * all_features, 1, sum)

library(gtools)

## our hidden underlying "mutation probability" for each 100-bp bin
hidden_covariate <- inv.logit(hidden_covariate_real)

## our fake "mutation counts" for each 100-bp bin
fake_counts <- sapply(hidden_covariate,
                      function(hidden_covariate_val) {
                          rbinom(n=1, size=1000, prob=hidden_covariate_val)
                      })

library(MASS)

mydata <- data.frame(all_features,
                     counts=fake_counts)


myres <- glm.nb(counts ~ .,data = mydata)
summary(myres)
