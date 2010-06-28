# commandArgs()
# q()


for (e in commandArgs()) {
  ta = strsplit(e,"=",fixed=TRUE)
  if(! is.na(ta[[1]][2])) {
    temp = ta[[1]][2]
    if(substr(ta[[1]][1],nchar(ta[[1]][1]),nchar(ta[[1]][1])) == "F") {
      FNstem <- temp
      cat("file name stem assigned ",FNstem,"\n")
    }
    if(substr(ta[[1]][1],nchar(ta[[1]][1]),nchar(ta[[1]][1])) == "M") {
      obsEdges <-  1 - as.numeric(temp)
      cat("observed fraction = ",obsEdges, "\n")
    }
  } 
}



library(ROCR)
# read in the file name stem from the command line
# FNstem <- commandArgs()[2]
# obsEdges <- 1-as.numeric(commandArgs()[3])

# get labels

lab <- read.table(paste(FNstem,".labels",sep=""))

scores <- list()
# our scores
sc0 <- read.table(paste(FNstem,".scores0",sep=""))
scores[[length(scores)+1]] <- list("score0",sc0)

sc1 <- read.table(paste(FNstem,".scores1",sep=""))
scores[[length(scores)+1]] <- list("score1",sc1)
sc1$V3 <- sc0$V3 + sc1$V3
scores[[length(scores)+1]] <- list("cum_score1",sc1)

sc2 <- read.table(paste(FNstem,".scores2",sep=""))
scores[[length(scores)+1]] <- list("scores2",sc2)
sc2$V3 <- sc1$V3 + sc2$V3
scores[[length(scores)+1]] <- list("cum_score2",sc2)



sc3 <- read.table(paste(FNstem,".scores3",sep=""))
sc3$V3 <- sc2$V3 + sc3$V3
scores[[length(scores)+1]] <- list("cum_score3",sc3)


sc4 <- read.table(paste(FNstem,".scores4",sep=""))
sc4$V3 <- sc3$V3 + sc4$V3
scores[[length(scores)+1]] <- list("cum_score4",sc4)


sc5 <- read.table(paste(FNstem,".scores5",sep=""))
sc5$V3 <- sc5$V3 + sc4$V3
scores[[length(scores)+1]] <- list("cum_score5",sc5)


sc6 <- read.table(paste(FNstem,".scores6",sep=""))
sc6$V3 <- sc6$V3 + sc5$V3
scores[[length(scores)+1]] <- list("cum_score6",sc6)

# competitor scores
scores[[length(scores)+1]] <- list("jaccard",   read.table(paste(FNstem,".jaccard",sep=""))  )
scores[[length(scores)+1]] <- list("dprod",   read.table(paste(FNstem,".dprod",sep=""))  )
scores[[length(scores)+1]] <- list("cneighb",   read.table(paste(FNstem,".cneighb",sep=""))  )
scores[[length(scores)+1]] <- list("hyperg",   read.table(paste(FNstem,".hyperg",sep=""))  )

fscVec <- c(obsEdges)
aucVec <- c(obsEdges)
# get all predictions

for (s in scores) {
	pred <- prediction(s[[2]]$V3, lab$V3)
	perf <- performance(pred,'auc')
	# aucVec[[length(aucVec)+1]] <- c(s[[1]], perf@y.values[[1]])
	aucVec[[length(aucVec)+1]] <- perf@y.values[[1]]
	perf <- performance(pred,'f')
	Fscore<-max(ifelse(is.nan(perf@y.values[[1]]),NA,perf@y.values[[1]]),na.rm=T);
	#fscVec[[length(fscVec)+1]] <- c(s[[1]], Fscore)
	fscVec[[length(fscVec)+1]] <- Fscore
}

# write performance files

if (! file.exists(paste(FNstem,".auc",sep=""))) {
	heads = c("ObservedFraction")
	for (s in scores) {
		heads <- c(heads,s[[1]])
	}
	write(heads, file = paste(FNstem,".auc",sep=""), ncolumns = length(heads), sep = "\t")
}

write(aucVec, file = paste(FNstem,".auc",sep=""), ncolumns = length(aucVec), sep = "\t", append=T)


if (! file.exists(paste(FNstem,".fsc",sep=""))) {
	heads = c("ObservedFraction")
	for (s in scores) {
		heads <- c(heads,s[[1]])
	}
	write(heads, file = paste(FNstem,".fsc",sep=""), ncolumns = length(heads), sep = "\t")
}
write(fscVec, file = paste(FNstem,".fsc",sep=""), ncolumns = length(fscVec), sep = "\t", append=T)


# quit
q()

