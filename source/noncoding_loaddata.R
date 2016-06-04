library(gtools)
library(MASS)
library(boot)

setwd("/home/tanjinxu/Project/noncoding")

## Split the dataset into training set and testing set
splitData <- function(mydata,split_rate){
  mcords <- strsplit(rownames(mydata),"-")
  chr <- lapply(mcords,"[",1)
  pos <- as.numeric(lapply(mcords,"[",2))
  chr_cords_list <- paste(chr,lapply(round(pos/10^6),as.character),sep="-")
  
  uniq_list <- unique(chr_cords_list)
  samplist <- 1:length(uniq_list)
  #train_indices_map <- sample(samplist,0.5*length(uniq_list))
  train_indices_map <- sample(samplist,split_rate*length(uniq_list))
  test_indices_map <- samplist[! samplist %in% train_indices_map]
  
  train_indices <- which(chr_cords_list %in% uniq_list[train_indices_map])
  test_indices <- which(chr_cords_list %in% uniq_list[test_indices_map])
  train_data <- mydata[train_indices,]
  train_data$rname <- rownames(mydata)[train_indices]
  test_data <- mydata[test_indices,]
  test_data$rname <- rownames(mydata)[test_indices]
  
  print(colnames(train_data))
  print(colnames(test_data))
  #print(test_data$rname)
  cat(sprintf("split-train:%s, test:%s\n", dim(train_data), dim(test_data)))
  return(list(train=train_data,test=test_data))
}

## Load the data file
## "result/cosmic.noncoding.mut.reg225.allsets.tsv"
## - the percentage of value of split_rate for training
## - the percentage of valud of (1-split_rate) for testing
loadData <- function(datfile, split_rate, outfl, normalized = FALSE, fromText = FALSE){
  # Read the features data file
  if (fromText){
    noncoding_mutations <- read.table(datfile,
                                      sep="\t",
                                      header=TRUE,
                                      stringsAsFactors=FALSE,
                                      row.names=1)
    save(noncoding_mutations,file="data/noncoding.mixed.331.rda")
  }else{
    # Load the feature data
    load(file="data/noncoding.mixed.331.rda")
  }
  
  # Get the features value
  # the last 2 column are response variable
  print(colnames(noncoding_mutations))
  features <- noncoding_mutations[,2:(ncol(noncoding_mutations)-2)]
  
  # Need to exclude the feautre dhs_avg_score as it has the same value 
  # as dhs_max_score so that the predicted coefficient value is N/A
  features <- features[, !names(features) %in% c("dhs_avg_score")]
  
  # Recurrent counts - dependent variable 
  samplefreq <- noncoding_mutations[,(ncol(noncoding_mutations)-1)]
  rawcounts <- noncoding_mutations[,ncol(noncoding_mutations)]
  
  print(range(samplefreq))
  ## Filter super "outliers"
  if (TRUE){
    used.indices <- which(samplefreq < 600)
    samplefreq <- samplefreq[used.indices]
    rawcounts <- rawcounts[used.indices]
    features <- features[used.indices,]
  }
  
  print(range(samplefreq))
  
  # Log transformation of tss distances, max_score
  #features$tss_dist <- sign(features$tss_dist) * log(abs(features$tss_dist)+0.5)
  features$tss_dist <- log(abs(features$tss_dist) + 1)
  score.index <- grep("max_score|avg_score", names(features))
  for (col in score.index){
    features[,col] <- log10(features[,col] + 1)
  }
  
  ##!!!
  ## Filter out those TF max_score/avg_score features
  filter.index <- grep("tss|gerp|gc|dhs|count", names(features))
  features.filtered <- features[,filter.index]
  features <- features.filtered
  print(colnames(features))
  
  if (normalized){
    feaMean <- apply(features,2,mean)
    feaStd <- apply(features,2,sd)
    meanMtx <- matrix(feaMean, nrow=nrow(features),ncol=length(feaMean),byrow=TRUE)
    stdMtx <- matrix(feaStd,nrow=nrow(features),ncol=length(feaStd),byrow=TRUE)
    feaNormed <- (features - meanMtx) / stdMtx
    features <- feaNormed
  }

  ## Now divide the dataset into training set and testing set
  #split_rate <- 0.80
  mydata <- data.frame(features,
                       samplefreq = samplefreq,
                       rawcounts=rawcounts)
  datlist <- splitData(mydata,split_rate)
  
  ## Save training/testing objects into file
  save(datlist,file=outfl)
}

## Sampling the data setcosmic.noncoding.mut.reg225.allsets.tsv
## Input data file: "train_test.rda"
sampleData <- function(sample_rate_train, sample_rate_test, datfl, outfl){
  # Load the data file
  load(file = datfl)
  train_data <- datlist$train
  test_data <- datlist$test

  cat(sprintf("train:%s, test:%s\n", dim(train_data), dim(test_data)))
  
  ## Sampling subset from the whole data uniformly
  ## The last column contains the hash value string:
  ## <chrom_id> + "_" + "<[coordinate/10^6]>"
  train_subset <- train_data[sample(nrow(train_data),
                                    sample_rate_train*nrow(train_data)),
                             1:ncol(train_data)]
  test_subset <- test_data[sample(nrow(test_data),
                                  sample_rate_test*nrow(test_data)),
                           1:ncol(train_data)]
  
  datlist<- list(train=train_subset, test=test_subset)
  
  # Write the sampled data to disk
  save(datlist, file=outfl)
}

## Load the data file - raw counts
## 80% for training; 20% for testing
#loadData("result/cosmic.noncoding.mut.reg331.allsets.tsv", 0.8, "data/train_test.filtered.trim.504.all.rda", FALSE, FALSE)

## Sample the data
## Sample 30% of the whole training set for training; and
## same 30% for testing
#sampleData(1.0, 1.0,"data/train_test.filtered.trim.504.all.rda", "data/train_test.sampled.filtered.trim.504.all.rda")

