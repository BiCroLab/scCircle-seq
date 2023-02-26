args<-commandArgs(T)
table<-read.table(args[1],header=FALSE,sep ="\t",quote="")
write.table(table,file=args[1],sep=',',col.names=FALSE,row.names=FALSE,quote=FALSE)