library(SNFtool)
K<-10;		# number of neighbors, usually (10~30)
alpha<-0.5;  	# hyperparameter, usually (0.3~0.8)
T<-20; 	# Number of Iterations, usually (10~20)

args<-commandArgs(trailingOnly=TRUE)
len<-length(args)

labelpath<-args[1]
filename<-args[2]

omics<-strsplit(tail(unlist(strsplit(filename, "/")), n=1), "_")
omics<-unlist(strsplit(unlist(omics)[4], "-"))

label<-as.matrix(read.csv(labelpath, row.names=1, header=T))
cluster_n <- length(unique(label))

data<-vector(mode='list', length=(len-2))
for(i in 1:(len-2)){
  omics_data<-as.matrix(read.csv(args[i+2], row.names=1, header=T))
  dist<-(dist2(omics_data,omics_data))^(1/2)
  w<-affinityMatrix(dist, K, alpha)
  #displayClusters(w,label);
  data[[i]]<-w
}

W = SNF(data, K, T)
result = spectralClustering(W,cluster_n);
#displayClusters(W, cluster)

cluster_mat<-as.matrix(result)
rownames(cluster_mat)<-rownames(label)
colnames(cluster_mat)<-c('cluster')
write.csv(cluster_mat, paste(filename, ".cluster.csv", sep=""),quote=FALSE )