#!/usr/bin/env Rscript

# The purpose of this script is to extract data from the splicemutr output

#------------------------------------------------------------------------------#
# loading libraries
print("loading libraries")

library(stringi)
library(stringr)
library(optparse)
library(data.table)
library(dplyr)
library(biomaRt)
library(ensembldb)
library(IRanges)
library(Homo.sapiens)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicFeatures)

#------------------------------------------------------------------------------#
# get introns
txdb<-TxDb.Hsapiens.UCSC.hg38.knownGene
introns<-intronsByTranscript(txdb,use.names=T)

#------------------------------------------------------------------------------#
# internal functions

check_junc<-function(junc,gene){
  #inputs:
  # junc: 3 element vector c(chr, start, end)
  # gene: 3 elements vector c(chr, start, end)
  #output:
  # junc_in: logical True if junction within gene coord
  
  same_chr<- junc[1] == gene[1]
  junc<-as.numeric(junc[2:3])
  gene<-as.numeric(gene[2:3])
  start_in<- junc[1] <= gene[2] & junc[1] >= gene[1]
  end_in<- junc[2]<= gene[2] & junc[2] >= gene[1]
  junc_in <- start_in & end_in
  return(junc_in)
}

same_junc<-function(junc1,junc2){ # don't really need to use this
  return(all(junc1==junc2))
}

split_block_sizes<-function(block_sizes){
  return(str_split(block_sizes,"[,]"))
}