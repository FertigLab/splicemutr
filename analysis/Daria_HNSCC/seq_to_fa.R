#!/usr/bin/env Rscript

# This script takes the validation_set.xlsx file and turns the amino acid sequences and metadata into fasta file. To be used when running polyester simulations.

#------------------------------------------------------------------------------#
# loading libraries

library(seqinr)
library(readxl)

#------------------------------------------------------------------------------#
# pre-processing the validation dataset

validation_data<-read_xlsx("/media/theron/My_Passport/head_and_neck_DARIA/data/IEDB_benchmarking/HLA-A02-01/validation_set.xlsx")
names_seq<-apply(validation_data,1,function(x){
  paste(paste(x[2],x[3],x[5],x[6],sep=";"),x[8],sep=">")
})
names_seq<-unique(names_seq)
names<-unname(vapply(names_seq,function(x){
  strsplit(x,">")[[1]][1]
},character(1)))
sequences<-unname(lapply(names_seq,function(x){
  strsplit(x,">")[[1]][2]
}))


#------------------------------------------------------------------------------#
# writing fasta file

write.fasta(sequences,names,file.out="/media/theron/My_Passport/head_and_neck_DARIA/data/IEDB_benchmarking/HLA-A02-01/transcripts.fa",
            open="w",nbchar=80,as.string=TRUE)
