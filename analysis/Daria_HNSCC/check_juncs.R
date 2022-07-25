#!/usr/bin/env Rscript

# The purpose of this script is to determine whether any of the junc files have "?" in strand

#------------------------------------------------------------------------------------------------------------------------------------------------#
# loading libraries


#------------------------------------------------------------------------------------------------------------------------------------------------#
# loading junc files

junc_files<-read.table("/media/theron/My_Passport/head_and_neck_DARIA/data/IEDB_benchmarking/HLA-A02-01/simulated_RNAseq_data/junc_files.txt")
qs<-vapply(seq(nrow(junc_files)),function(x){
  junc<-read.table(junc_files[x,])
  any(junc[,6] == "?")
},logical(1))

chr<-lapply(seq(nrow(junc_files)),function(x){
  junc<-read.table(junc_files[x,])
  junc[,1]
})
chr<-unique(unlist(chr))