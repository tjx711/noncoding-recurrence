library(car)
library(ggplot2)

##########################################################################
# This program checks whether the counts data fits some distributions,   #
# i.e., poisson distribution or negative binomial distribution           #
#                                                                        #
# The measurements mainly includs the plots(counts,frequencies) and also #
# the chi-square test, var/mean ratio,etc.                               #
#                                                                        #
# Tanjin Xu, 03/2016                                                     #
##########################################################################

setwd("/home/tanjinxu/Project/noncoding")
load(file="data/train_test.filtered.430.all.rda")
#load(file="data/for_R/train_test.sampled.filtered.0.8.4.421.rda")
train_data <- datlist$train
test_data <- datlist$test
#rawcounts <-c(train_data$rawcounts,test_data$rawcounts)
samplefreq <-c(train_data$samplefreq,test_data$samplefreq)

r <- c(mean(samplefreq),var(samplefreq))
c(mean=r[1],var=r[2],ratio=r[2]/r[1])


mutCountFreq <- data.frame(table(samplefreq))
suffix <- "all.png"

png(paste("plots/noncoding.freq.counts.log.dist",suffix,sep="."), width=400, height=400)
x <- as.matrix(as.numeric(levels(mutCountFreq$samplefreq)))
y <- as.matrix(sapply(mutCountFreq["Freq"],as.numeric))
plot(x,y,log="y",xlab="reccurrence(all)",ylab="log(frequency)",main="(b) Noncoding Mutation Log-frequency")
dev.off()

png(paste("plots/noncoding.freq.counts.dist.100",suffix,sep="."), width=400, height=400)
indices <- which(x<100)
plot(x[indices],y[indices],xlab="reccurrence(<100)",ylab="frequency",main="(c) Noncoding Mutation Frequency Dist")
dev.off()

#png("plots/noncoding.freq.counts.log.dist.png", width=400, height=400)
#plot(log2(x),log2(y),xlab="log(recurrence)",ylab="log(frequency)",main="(c) Noncoding Mutation Log Frequency")
#dev.off()

samplefreq <- samplefreq[which(samplefreq<100)]
k <- unique(samplefreq)
(x_p=table(samplefreq))
png(paste("plots/noncoding.freq.counts.linear.100",suffix,sep="."),width=400, height=400) 
plot(k,log(x_p)+lfactorial(k),xlab="recurrence (< 100)",ylab="log(freq)+lfactorial(recurrence)",yaxt='n',
     main="(d) Noncoding Recurrence (<100) ~ Poisson")
abline(lm((log(x_p)+lfactorial(k))~k), col="red")
lines(lowess(k,(log(x_p)+lfactorial(k))), col="blue")
legend('topleft',c("linear model", "lowess"), lty = 1, col=c('red','blue'),bty='n',cex=.75)
dev.off()


##qq-plot
#mut.subset <- rawcounts[which(rawcounts<100)]
#k.sub <- unique(mut.subset)
#(x.p.sub=table(mut.subset)) 
#plot(k.sub,log(x.p.sub)+lfactorial(k.sub),xlab="counts(freq<100)",ylab="log(freq)+lfactorial(counts)",
#     main="Counts ~ Poisson Distribution")
#abline(lm((log(x.p.sub)+lfactorial(k.sub))~k.sub), col="red")
#lines(lowess(k.sub,(log(x.p.sub)+lfactorial(k.sub))), col="blue")

if (FALSE){
y <- rt(200, df = 5)
qqnorm(y); qqline(y, col = 2)
qqplot(y, rt(300, df = 5))

qqnorm(precip, ylab = "Precipitation [in/yr] for 70 US cities")
}




