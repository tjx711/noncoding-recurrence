library(plyr)
library(ggplot2)

setwd("/home/tanjinxu/Project/noncoding/plots")
#metrics <- c("mse","mae","mde","expvar",'r2')
#glm.models <- c("ada","gradient","rf")

metrics <- c("mse")
glm.models <- c("ann","sae")

sampling.rates <- c(0.3,0.4,0.5)
num.rates <- length(sampling.rates)

titles <- list(mse="Mean Squared Error",
               mae="Mean Absolute Error",
               mdae="Median Absolute Error",
               expvar="Explained Variance Score",
               r2="R-squared Score")

### plotMetric function
plotMetricByRate <- function(glm.res,metric,rate){
  colnames(glm.res) <- c("rate","model","param","metric","sd")
  png(paste("nc.freq.res.net",metric,rate,"png",sep="."),width=6, height=4, units='in',res=350)
  pd <- position_dodge(0.1)
  gtitle <- sprintf("Training %s (sampling=%.1f)",titles[[metric]],rate)
  print(ggplot(glm.res, aes(x=param,y=metric,colour=model,group=model,shape=model))+
          theme_bw() +
          theme(axis.line = element_line(colour = "black"),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                #panel.border = element_blank(),
                panel.background = element_blank())+
          ggtitle(gtitle)+
          ylab(metric)+
          geom_errorbar(aes(ymin=metric-sd, ymax=metric+sd), colour="grey",width=.1, position=pd) +
          geom_line(position=pd)+
          geom_point(position=pd,size=3))
  dev.off()
}

plotMetricByModel <- function(glm.res,metric,model){
  colnames(glm.res) <- c("rate","model","param","metric","sd")
  glm.res$rate <- as.character(glm.res$rate)
  png(paste("nc.freq.res.net",metric,model,"png",sep="."),width=6, height=4, units='in',res=350)
  pd <- position_dodge(0.1)
  gtitle <- sprintf("Training %s (%s)",titles[[metric]],model)
  print(ggplot(glm.res, aes(x=param,y=metric,colour=rate,group=rate,shape=rate))+
          theme_bw() +
          theme(axis.line = element_line(colour = "black"),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                #panel.border = element_blank(),
                panel.background = element_blank())+
          ggtitle(gtitle)+
          ylab(metric)+
          geom_errorbar(aes(ymin=metric-sd, ymax=metric+sd), colour="grey",width=.1, position=pd) +
          geom_line(position=pd)+
          geom_point(position=pd,size=3))
  dev.off()
}

plotMetric <- function(glm.res,metric){
  colnames(glm.res) <- c("rate","model","metric","sd")
  png(paste("nc.freq.res.net.min",metric,"png",sep="."),width=6, height=4, units='in',res=350)
  pd <- position_dodge(0.03)
  gtitle <- sprintf("Training Min_%s",titles[[metric]])
  print(ggplot(glm.res, aes(x=rate,y=metric,colour=model,group=model,shape=model))+
          theme_bw() +
          theme(axis.line = element_line(colour = "black"),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                #panel.border = element_blank(),
                panel.background = element_blank())+
          ggtitle(gtitle)+
          ylab(metric)+
          geom_errorbar(aes(ymin=metric-sd, ymax=metric+sd), colour="grey",width=.05, position=pd) +
          geom_line(position=pd)+
          geom_point(position=pd,size=3))
  dev.off()
}

### Main function
for (metric in metrics){
  train.file <- ""
  if (metric == "mse"){
    #train.file <- paste("result/nc.res.result.train.ensemble",metric,sep=".")
    train.file <- paste("result/nc.res.result.train.net",metric,sep=".")
  }
  else{
    train.file <- paste("result/nc.res.result.train",metric,sep=".")
  }
  nc.res.mse <- read.table(train.file,
                           sep="\t",
                           header=F,
                           stringsAsFactors=FALSE)
  
  train.colname <- paste("test",metric,sep=".")
  shuff.colname <- paste("shuff",metric,sep=".")
  colnames(nc.res.mse) <- c("rate","model","p1","p2","p3","p4","p5","p6")

  res.cols.names <- c("model","rate","metric","sd")
  param.cols <- c("p1","p2","p3","p4","p5","p6")
  glm.res <- c()
  glm.res.min <- c()
  for (model in glm.models){
    nc.res.mse.model <- nc.res.mse[nc.res.mse$model == model,]
    train.avg <- aggregate(cbind(p1,p2,p3,p4)~rate+model, data = nc.res.mse.model, mean)
    train.sd <- aggregate(cbind(p1,p2,p3,p4)~rate+model, data = nc.res.mse.model, sd)
    
    tmp.col.rates <- rep(sampling.rates, each=length(param.cols))
    tmp.col.model <- rep(model,length(sampling.rates) * length(param.cols))
    tmp.col.avg <- as.vector(t(train.avg[,3:(length(param.cols)-1)]))
    tmp.col.sd <- as.vector(t(train.sd[,3:(length(param.cols)-1)]))
    tmp.col.params <- rep(param.cols,length(sampling.rates))
    glm.val.model <- data.frame(tmp.col.rates, tmp.col.model,tmp.col.params,tmp.col.avg,tmp.col.sd,
                                stringsAsFactors = F)
    
    ## Plot metrics by model
    plotMetricByModel(glm.val.model,metric,model)
    glm.res <- rbind(glm.val.model,glm.res)
    
    if (metric == "mse" || metric == "mae"){
      train.avg.ps <- train.avg[,3:(length(param.cols)-1)]
      train.sd.ps <- train.sd[,3:(length(param.cols)-1)]
      values.min <- apply(train.avg.ps,1,min)
      indices.min <- which(train.avg.ps == values.min, arr.ind=T)
      glm.res.min.model <- data.frame(train.avg[,1:2], train.avg.ps[indices.min], train.sd.ps[indices.min])
      glm.res.min <- rbind(glm.res.min.model, glm.res.min)
    }
  }
  
  if (metric == "mse" || metric == "mae"){
    plotMetric(glm.res.min,metric)
  }
  
  for (ra in sampling.rates){
    colnames(glm.res) <- c("rate","model","param","metric","sd")
    glm.res.rate <- glm.res[which(glm.res$rate==ra),]
    plotMetricByRate(glm.res.rate,metric,ra)
  }
}

  
