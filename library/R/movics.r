library("MOVICS")

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

#"CIMLR"
#"MoCluster"
method_list<-list("PINSPlus", "NEMO", "COCA", "LRAcluster", "ConsensusClustering","IntNMF")
moic.res.list <-getMOIC(data = dat, methodslist = method_list, N.clust = group_n,type = omics_type)

movics_n<-strsplit(filename, split = "[.]")[[1]][2]
nametag<-paste(filename, ".", sep="")
for (i in 1:length(moic.res.list)){
    df<-moic.res.list[[i]]$clust.res
    df$samID<-NULL
    names(df)[1]<-'cluster'
    path<-paste(nametag, moic.res.list[[i]]$mo.method, sep="")
    write.csv(df, paste(path, ".cluster.csv", sep=""),quote=FALSE )
}