library(plyr)
library(ggplot2)

setwd("/home/tanjinxu/Project/noncoding/plots")

glm.models <- c("nb","poisson")
#glm.models <- c("poisson","nb","rf","gradient","ann","sae","ada")
#sampling.rates <- c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8)
#num.rates <- length(sampling.rates)

features <- c("tfbs_count",
              "atf3_count",
              "cebpb_count",
              "cebpd_count",
              "creb1_count",
              "egr1_count",
              "ets1_count",
              "maff_count",
              "dhs_src_count",
              "dhs_max_score",
              "gerp_score",
              "tss_dist",
              "gc_percent")

st.err <- function(x) {
  sd(x)/sqrt(length(x))
}


imp.file <- "result/nc.res.result.glm.imp"
nc.res.imp <- read.table(imp.file,
                         sep=",",
                         header=F,
                         stringsAsFactors=FALSE)
colnames(nc.res.imp) <- c("model","rate",features)

nc.res.imp.avg <- aggregate(. ~ model+rate,data=nc.res.imp, mean)
nc.res.imp.all <- nc.res.imp.avg[,c("model",features)]
nc.res.imp.all.avg <- aggregate(nc.res.imp.all[,features], by=list(nc.res.imp.all$model),FUN=mean)
nc.res.imp.all.sd <- aggregate(nc.res.imp.all[,features], by=list(nc.res.imp.all$model),FUN=st.err)

indices <- order(-nc.res.imp.all.avg[1,features])
features <- features[indices]

tmp.model.cols <- rep(glm.models,each=13)
tmp.feature.cols <- rep(features,2)
tmp.avg.cols <- as.vector(t(nc.res.imp.all.avg[,features]))
tmp.sd.cols <- as.vector(t(nc.res.imp.all.sd[,features]))
nc.res.imp.val <- data.frame(tmp.model.cols,tmp.feature.cols,tmp.avg.cols,tmp.sd.cols,stringsAsFactors = F)
colnames(nc.res.imp.val) <- c("model","feature","importance_score","sd")

png("nc.freq.res.imp.avg.all.png",width=6, height=4, units='in',res=350)
pd <- position_dodge(0.2)
print(ggplot(nc.res.imp.val, aes(x=factor(feature,levels=unique(feature)),y=importance_score,colour=model,group=model,shape=model))+
        theme_bw() +
        theme(axis.line = element_line(colour = "black"),
              axis.text.x = element_text(angle=50,size=9,vjust=0.5),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              #panel.border = element_blank(),
              panel.background = element_blank())+
        ggtitle("Variable Importance Score")+
        xlab("feature")+
        ylab("importance_score")+
        geom_errorbar(aes(ymin=importance_score-sd, ymax=importance_score+sd), colour="grey",width=.03, position=pd) +
        geom_line(position=pd)+
        geom_point(position=pd,size=3))
dev.off()
