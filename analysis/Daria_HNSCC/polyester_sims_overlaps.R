#!/usr/bin/env Rscript

# The purpose of this script is to simulate RNAseq for sequences used for validation of splicemutr

#------------------------------------------------------------------------------#
# loading libraries and setting seed

library(polyester)
library(Biostrings)
library(stringr)
library(ensembldb)
library(biomaRt)
library(ggplot2)

#------------------------------------------------------------------------------#
# internal functions

obtain_surrounding<-function(num,surr,max){
  min_num<-num-surr
  max_num<-num+surr
  if (min_num<0){
    min_num<-0
    max_num<-max_num+abs(min_num)
  }
  if (max_num>max){
    max_num<-max
  }
  sequence<-seq(min_num,max_num)
  which_seq<-sequence!=num
  return(sequence[which_seq])
}

#------------------------------------------------------------------------------#
# loading in the fasta files

trancsripts<-"/media/theron/My_Passport/reference_genomes/SEQUENCES/GENCODE/gencode.v36.pc_transcripts.fa"

trancsripts<-readDNAStringSet(trancsripts)

#------------------------------------------------------------------------------#
# creating new fasta file

fasta_dir<-"/media/theron/My_Passport/head_and_neck_DARIA/data/IEDB_benchmarking/overlaps.fa"
writeXStringSet(trancsripts[1:500],fasta_dir,append=FALSE)

overlap_fasta<-readDNAStringSet(fasta_dir)

#------------------------------------------------------------------------------#
# setting up the uniform counts matrix

num_samples = 24 # 6 per normal and not normal
countmat = matrix(200, nrow=length(overlap_fasta), ncol=num_samples)

#------------------------------------------------------------------------------#
# modifying the uniform counts matrix for differential splicing

percentages<-sort(rep(c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8),3))
reads<-500*200
samples<-lapply(percentages,function(x){
  sample(seq(500),x*500)
})
reads_per_tx<-vapply(percentages,function(x){
  return(round(reads/((1-x)*500)))
},numeric(1))

for(x in seq(length(samples))){
  countmat[samples[[x]],x]<-0
  countmat[which(countmat[,x]>0),x]<-round(rnorm(500-length(samples[[x]]),mean=reads_per_tx[x],sd=5))
}

# countmat_ggplot<-data.frame(c(countmat[1:20,1],countmat[1:20,12]))
# 
# countmat_ggplot$isoform<-NA
# countmat_ggplot$isoform[1:20]<-1:20
# countmat_ggplot$isoform[21:40]<-1:20
# countmat_ggplot$group<-NA
# countmat_ggplot$group[1:20]<-"normal"
# countmat_ggplot$group[21:40]<-"tumor"
# 
# colnames(countmat_ggplot)<-c("reads","isoform","group")
# 
# ggplot(data=countmat_ggplot, aes(x=isoform, y=reads, group=group)) +
#   geom_line(aes(color=group))+
#   geom_point(aes(color=group))

#------------------------------------------------------------------------------#
# simulate reads:

simulate_experiment_countmat(fasta_dir, readmat=countmat, paired = T, seed = 27,
                             outdir='/media/theron/My_Passport/head_and_neck_DARIA/data/IEDB_benchmarking/HLA-A02-01/simulated_RNAseq_data_overlaps')
saveRDS(countmat,file="/media/theron/My_Passport/head_and_neck_DARIA/data/IEDB_benchmarking/simulated_RNAseq_data_overlaps/countmat.rds")

