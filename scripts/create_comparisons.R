#!/usr/bin/env Rscript

library(dplyr)
library(stringr)
library(optparse)

#------------------------------------------------------------------------------#
# handling command line input

arguments <- parse_args(OptionParser(usage = "",
                 description="",
                 option_list=list(
                   make_option(c("-d","--splice_dat_file"),
                               default = sprintf("%s",getwd()),
                               help="splicemutr file"),
                   make_option(c("-s","--sample_files"),
                               default = sprintf("%s",getwd()),
                               help="the sample files"),
                   make_option(c("-o","--out_prefix"),
                               default = sprintf("%s",getwd()),
                               help="the output prefix"))))
opt=arguments
splice_dat_file <- opt$splice_dat_file
sample_files<-read.table(opt$sample_files)
out_prefix <- opt$out_prefix

#------------------------------------------------------------------------------#
# internal functions

parse_kmers <- function(kmers){
  a<-unlist(lapply(strsplit(as.character(kmers[seq(3,length(kmers))]),":"),function(kmer_split){
    kmer_split <- unique(kmer_split)
    length(kmer_split)-length(which(str_detect(kmer_split,"NAZZZZZZZ")))
  }))
  a<-c(as.character(kmers[seq(1,2)]),a)
  return(a)
}

#------------------------------------------------------------------------------#
# reading in the data

splice_dat <- readRDS(splice_dat_file)
splice_dat$juncs <- sprintf("%s:%s:%s:%s",splice_dat$chr,splice_dat$start,splice_dat$end,splice_dat$strand)

#------------------------------------------------------------------------------#
# creating kmer files

kmers <- data.frame(seq(nrow(splice_dat)))
colnames(kmers)<-"rows"
kmers$juncs <- splice_dat$juncs

for (file in sample_files$V1){
  juncs <- as.data.frame(splice_dat$juncs)
  if (file.size(file)!=0){
    sample_dat <- readRDS(file)
    sample <- str_remove(basename(file),"_splicemutr_kmers.rds")
    kmers[,sample]<-sample_dat$kmers
  } else {
    kmers[,sample]<-""
  }
}

splice_dat$rows <- seq(nrow(splice_dat))
kmers_parsed<-as.data.frame(t(apply(kmers,1,parse_kmers)))
colnames(kmers_parsed) <- colnames(kmers)

#------------------------------------------------------------------------------#
# saving comparison info

splice_dat$deltapsi<-as.numeric(splice_dat$deltapsi)
splice_dat_norm <- splice_dat %>% dplyr::filter(deltapsi<0)
splice_dat_tumor <- splice_dat %>% dplyr::filter(deltapsi>0)
saveRDS(splice_dat_norm,file=sprintf("%s_splice_dat_norm.rds",out_prefix))
saveRDS(splice_dat_tumor,file=sprintf("%s_splice_dat_tumor.rds",out_prefix))
saveRDS(kmers_parsed,file=sprintf("%s_kmers.rds",out_prefix))
