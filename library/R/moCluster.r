library("MOVICS")
library(mogsa)

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

for(i in 3:len){
  dat[[i-2]]<-t(as.matrix(read.csv(args[i], row.names=1, header=T)))
  omics_type[[i-2]]<-"gaussian"
}

method_list<-list("MoCluster")
moic.res.list <-getMOIC(data = dat, methodslist = method_list, N.clust = group_n,type = omics_type)
#save(moic.res.list, file = paste(filename, ".rda", sep=""))

movics_n<-strsplit(filename, split = "[.]")[[1]][2]

df<-moic.res.list$clust.res
df$samID<-NULL
names(df)[1]<-'cluster'
write.csv(df, paste(filename, ".cluster.csv", sep=""),quote=FALSE )