suppressMessages(library(ggplot2))
suppressMessages(library(MOFA2))

args<-commandArgs(trailingOnly=TRUE)
len<-length(args)

group<-args[1]
filename<-args[2]

omics<-strsplit(tail(unlist(strsplit(filename, "/")), n=1), "_")
omics<-unlist(strsplit(unlist(omics)[4], "-"))

group<-as.matrix(read.csv(group, row.names=1, header=T))
group_n <- length(unique(group))

dat<-vector(mode='list', length=(len-2))
for(i in 1:(len-2)){
  dat[[i]]<-t(as.matrix(read.csv(args[i+2], row.names=1, header=T)))
}

names(dat)<-omics

# build MOFA object
MOFAobject <- create_mofa(dat)
DataOptions <- get_default_data_options(MOFAobject)
ModelOptions <- get_default_model_options(MOFAobject)
TrainOptions <- get_default_training_options(MOFAobject)
TrainOptions$convergence_mode<-"fast"

# Automatically drop factors that explain less than 2% of variance in all omics

MOFAobject <- prepare_mofa(
  object = MOFAobject,
  data_options = DataOptions,
  model_options = ModelOptions,
  training_options = TrainOptions
)

hdf5<-paste(filename, ".hdf5", sep="")
MOFAobject.trained <- run_mofa(MOFAobject, outfile=hdf5)

# Loading an existing trained model
model <- load_model(hdf5)

clusters <- cluster_samples(model, k=group_n)
predicted_samples <- data.frame(clusters[1])
colnames(predicted_samples)<-c('cluster')
write.csv(predicted_samples, paste(filename, ".cluster.csv", sep=""),quote=FALSE )