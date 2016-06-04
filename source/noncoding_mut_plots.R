library(plyr)
library(ggplot2)

setwd("/home/tanjinxu/Project/noncoding/plots")
metrics <- c("mse","mae","mdae","expvar","r2")
#metrics <- c("mse")

#glm.models <- c("poisson","quasi","nb")
#sampling.rates <- c(0.3,0.4,0.5,0.6,0.7,0.8)

#glm.models <- c("ann","sae")
glm.models <- c("nb","poisson","ann","sae")
#glm.models <- c("poisson","nb","rf","gradient","ann","sae","ada")
sampling.rates <- c(0.3, 0.4, 0.5)

st.err <- function(x) {
  sd(x)/sqrt(length(x))
}

num.rates <- length(sampling.rates)

titles <- list(mse="Mean Squared Error",
               mae="Mean Absolute Error",
               mdae="Median Absolute Error",
               expvar="Explained Variance Score",
               r2="R-squared Score")

for (metric in metrics){
  test.file <- paste("result/nc.res.result.test",metric,sep=".")
  nc.res.mse <- read.table(test.file,
                           sep="\t",
                           header=F,
                           stringsAsFactors=FALSE)
  
  test.colname <- paste("test",metric,sep=".")
  shuff.colname <- paste("shuff",metric,sep=".")
  colnames(nc.res.mse) <- c("rate","model",test.colname, shuff.colname)

  tmp <- nc.res.mse[which(nc.res.mse$rate < 0.6),]
  nc.res.mse <- tmp
  
  res.cols.names <- c("model","rate","metric","sd")
  glm.res <- c()
  for (model in glm.models){
    nc.res.mse.model <- nc.res.mse[nc.res.mse$model == model,]
    nc.num.rows <- nrow(nc.res.mse.model)
      
    test.avg <- aggregate(nc.res.mse.model[,test.colname], by=list(nc.res.mse.model$rate), FUN=mean)[2]
    test.sd <- aggregate(nc.res.mse.model[,test.colname], by=list(nc.res.mse.model$rate), FUN=st.err)[2]
    #test.sd <- test.sd / sqrt(nc.num.rows)
    
    shuff.avg <- aggregate(nc.res.mse.model[,shuff.colname], by=list(nc.res.mse.model$rate), FUN=mean)[2]
    shuff.sd <- aggregate(nc.res.mse.model[,shuff.colname], by=list(nc.res.mse.model$rate), FUN=st.err)[2]
    #shuff.sd <- shuff.sd / sqrt(nc.num.rows)
  
    glm.val.model <- cbind(rep(model,num.rates),sampling.rates, test.avg, test.sd)
    shuffl.val.model <- cbind(rep(paste(model,"shuff",sep="."),num.rates),sampling.rates, shuff.avg, shuff.sd)
  
    names(glm.val.model) <- res.cols.names
    names(shuffl.val.model) <- res.cols.names
    #glm.res <- rbind(glm.val.model, shuffl.val.model,glm.res)
    glm.res <- rbind(glm.val.model,glm.res)
  }

  png(paste("nc.freq.res.test2.all",metric,"png",sep="."),width=6, height=4, units='in',res=350)
  pd <- position_dodge(0.1)
  colnames(glm.res) <- res.cols.names
  print(ggplot(glm.res, aes(x=rate,y=metric,colour=model,group=model,shape=model))+
          theme_bw() +
          theme(axis.line = element_line(colour = "black"),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                #panel.border = element_blank(),
                panel.background = element_blank())+
    ggtitle(paste("Test",titles[[metric]]))+
    ylab(metric)+
    geom_errorbar(aes(ymin=metric-sd, ymax=metric+sd), colour="grey",width=.03, position=pd) +
    geom_line(position=pd)+
    geom_point(position=pd,size=3))
  dev.off()
}

  
