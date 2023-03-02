suppressMessages(library(iClusterPlus))
suppressMessages(library(GenomicRanges))
suppressMessages(library(gplots))
suppressMessages(library(lattice))
suppressMessages(library(caret))
suppressMessages(library(preprocessCore))

args<-commandArgs(trailingOnly=TRUE)
len<-length(args)

group<-args[1]
filename<-args[2]

omics<-strsplit(tail(unlist(strsplit(filename, "/")), n=1), "_")
omics<-unlist(strsplit(unlist(omics)[4], "-"))

group<-as.matrix(read.csv(group, row.names=1, header=T))
group_n <- length(unique(group))

dat<-vector(mode='list', length=(len-2))
omics_type<-vector(mode='list', length=(len-2))
alpha<-vector(mode='list', length=(len-2))

for(i in 3:len){
  dat[[i-2]]<-as.matrix(read.csv(args[i], row.names=1, header=T))
  omics_type[[i-2]]<-"gaussian"
  alpha[[i-1]]<-1
}

k<-group_n-1
if (length(omics_type) == 1){
    cv.fit <- iClusterPlus(dt1=dat[[1]], type=omics_type, K=k, maxiter=100, n.burnin=100, n.draw=200, sdev=0.05, eps=1.0e-4)
}else if (length(omics_type) == 2){
    cv.fit <- iClusterPlus(dt1=dat[[1]], dt2=dat[[2]], type=omics_type, K=k, maxiter=100, n.burnin=100, n.draw=200, sdev=0.05, eps=1.0e-4)
}else if(length(omics_type) == 3){
    cv.fit <- iClusterPlus(dt1=dat[[1]], dt2=dat[[2]], dt3=dat[[3]], type=omics_type, K=k, maxiter=100, n.burnin=100, n.draw=200, sdev=0.05, eps=1.0e-4)
}else if(length(omics_type) == 4){
    cv.fit <- iClusterPlus(dt1=dat[[1]], dt2=dat[[2]], dt3=dat[[3]], dt4=dat[[4]],type=omics_type, K=k, maxiter=100, n.burnin=100, n.draw=200, sdev=0.05, eps=1.0e-4)
}

result<-cv.fit
cluster_mat<-as.matrix(result$clusters)
rownames(cluster_mat)<-rownames(group)
colnames(cluster_mat)<-c('cluster')
write.csv(cluster_mat, paste(filename, ".cluster.csv", sep=""),quote=FALSE )
