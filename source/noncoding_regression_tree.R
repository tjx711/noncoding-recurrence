library(rpart)

setwd("/home/tanjinxu/Project/noncoding")

load(file="train_test.rda")
train_data <- datlist$train
test_data <- datlist$test

######################
## Apply Regression Tree
nc.tree <- rpart(counts ~ ., data = train_data, cp = 10^(-6))

names(nc.tree)
min_cp <- nc.tree$cptable[which.min(nc.tree$cptable[,"xerror"]),"CP"]
nc.tree.pruned <- prune(nc.tree,cp=min_cp)
plot(nc.tree.pruned)
png("nctree.pruned.png",width=2400,height=1200)
post(nc.tree.pruned, file="", title.="Pruned Regression Truee for Mutation Counts Frequency")
dev.off()

#TODO: How to prune the tree? Any criteria?
nc.tree$cptable[1:10,]
nc.tree$cptable[dim(nc.tree$cptable)[1] - 9:0,]

cp9 <- which(nc.tree$cptable[,2] == 9)
nc.tree9 <- prune(nc.tree, nc.tree$cptable[cp9, 1])
print(nc.tree9)
summary(nc.tree9)

png("nctree9.png",width=1200,height=800)
post(nc.tree9, file="", title.="Pruned Regression Truee for Mutation Counts Frequency",bp=18)
dev.off()

# Predict output
nc.tree.predicted <- predict(nc.tree,test_data)


########################
## Apply Random Forest
library(randomForest)
nc.rf <- randomForest(counts ~ ., data=train_data, importance=TRUE)
png("ncrf.png",width=1200,height=800)
plot(nc.rf)
dev.off()

png("ncrf_barplots.png")
par(mfrow = c(2, 1))
barplot(nc.rf$importance[, 7], main = "Importance (Dec.Accuracy)")
barplot(nc.rf$importance[, 8], main = "Importance (Gini Index)")
dev.off()

summary(nc.rf)
importance(nc.rf)

# Predict output
nc.rf.predicted <- predict(nc.rf,test_data)
