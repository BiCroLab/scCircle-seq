args<-commandArgs(T)
table<-read.table(args[1],header=FALSE,sep ="\t",quote="")
table_new<-table[,c(1,2,3)]
write.table(table_new,file='for_count.bed',sep='\t',col.names=TRUE,row.names=FALSE,quote=FALSE)
table_new <- paste(table[,1],":",table[,2], "-", table[,3],sep='')
write.table(table_new,file='name.matrix',sep='\t',col.names=TRUE,row.names=FALSE,quote=FALSE)
