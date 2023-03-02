#Spectrum is a fast adaptive spectral clustering method for single or multi-view data
library("Spectrum")
args<-commandArgs(trailingOnly=TRUE)
len<-length(args)

group<-args[1]
filename<-args[2]

omics<-strsplit(tail(unlist(strsplit(filename, "/")), n=1), "_")
omics<-unlist(strsplit(unlist(omics)[4], "-"))

group<-as.matrix(read.csv(group, row.names=1, header=T))
group_n <- length(unique(group))

data<-vector(mode='list', length=(len-2))
for(i in 1:(len-2)){
  data[[i]]<-t(as.matrix(read.csv(args[i+2], row.names=1, header=T)))
}

result <- Spectrum(data,missing=TRUE,method=3, fixk=group_n)

cluster_mat<-as.matrix(result[[1]])
rownames(cluster_mat)<-rownames(group)
colnames(cluster_mat)<-c('cluster')
write.csv(cluster_mat, paste(filename, ".cluster.csv", sep=""),quote=FALSE )