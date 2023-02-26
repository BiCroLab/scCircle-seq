args<-commandArgs(T)
bed<-read.table(args[1],header=FALSE,sep='\t')
bed[,4]<-'CNVkit'
bed <- bed[which(bed[,5] != 0),]
length <- bed[,3] - bed[,2]
bed[,5] <- (bed[,5]/length)*mean(length)
bed[,5] <- bed[,5]/min(bed[,5])
write.table(bed,file=args[2],sep='\t',col.names=FALSE,row.names=FALSE,quote=FALSE)
