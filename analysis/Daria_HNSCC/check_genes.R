#!/usr/bin/env Rscript

# purpose of this script is to check which genes have been determined to be differentially used
# The target genes that I am expressing are differentially expressed, but not differentially spliced. This is shown 

#------------------------------------------------------------------------------#
# loading libraries

library(polyester)
library(Biostrings)
library(stringr)
library(readxl)

#------------------------------------------------------------------------------#
# loading in the leafcutter and validation fasta

load("/media/theron/My_Passport/head_and_neck_DARIA/data/IEDB_benchmarking/HLA-A02-01/data.Rdata")
validation_fasta<-readDNAStringSet("/media/theron/My_Passport/head_and_neck_DARIA/data/IEDB_benchmarking/HLA-A02-01/validation.fa")

#------------------------------------------------------------------------------#
# processing the validation_fasta for a list of target and random genes

names<-validation_fasta@ranges@NAMES
genes<-unname(vapply(names,function(x){
  str_split(x,"[|]")[[1]][6]
},character(1)))

genes_target<-unique(genes[seq(119)])
genes_random<-unique(genes[seq(119+1,length(genes))])

#------------------------------------------------------------------------------#
# determining percentages of random vs target gene representation in genes determined to be differentially spliced

diff_genes<-unique(introns$gene)
target_per<-length(which(genes_target %in% diff_genes))/length(genes_target) # 0.04545455, 01_05_2021 run
random_per<-length(which(genes_random %in% diff_genes))/length(genes_random) # 0.4761905, 01_05_2021 run


