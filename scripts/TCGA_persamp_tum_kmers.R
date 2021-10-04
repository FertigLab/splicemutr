#!/usr/bin/env Rscript

#author: "Theron Palmer"
#date: "08/06/2021"

library(stringr)
library(dplyr)
library(optparse)

#------------------------------------------------------------------------------#
# handling command line input

arguments <- parse_args(OptionParser(usage = "",
                                     description="",
                                     option_list=list(
                                       make_option(c("-t","--tumor_dir"),
                                                   default = sprintf("%s",getwd()),
                                                   help="tumor_dir"),
                                       make_option(c("-c","--cancer"),
                                                   default = sprintf("%s",getwd()),
                                                   help="cancer"),
                                       make_option(c("-n","--tum_norm"),
                                                   default = sprintf("%s",getwd()),
                                                   help="tum_norm"))))
opt=arguments
tumor_dir <- opt$tumor_dir
cancer <- opt$cancer
tum_norm_file <- opt$tum_norm

#------------------------------------------------------------------------------#
# internal functions

filter_tumor_kmers <- function(row){
  normal_kmers<-strsplit(as.character(row[length(row)]),":")[[1]]
  normal_kmers<-normal_kmers[normal_kmers!=""]
  kmers_split <- lapply(seq(ncol(row)),function(col_val){
    a<-strsplit(as.character(row[col_val]),":")[[1]]
    a<-vapply(a,function(val){
      str_replace_all(val,"Z","")
    },character(1))
    a<-a[a !="" & a != "NA"]
    paste(a[!(a %in% normal_kmers)],collapse=":")
  })
  return(unlist(kmers_split))
}

#------------------------------------------------------------------------------#
# reading in the necessary files

genotypes_file <- sprintf("%s/%s/%s_genotypes.txt",tumor_dir,cancer,cancer)
genotypes <- read.table(genotypes_file,header=T)
tumor_samples <- genotypes$type == "T"

kmer_file <- sprintf("%s/kmer_files.txt",tumor_dir)
kmer_files <- read.table(kmer_file)

splice_dat_file <- sprintf("%s/%s/%s_splicemutr_dat.txt",tumor_dir,cancer,cancer)
splice_dat <- read.table(splice_dat_file,header=T,sep="\t",)

tum_norm_kmers <-read.table(tum_norm_file,header=T,sep="\t")

for (i in nrow(kmer_files)){
  single_kmer_file <- as.character(kmer_files[i,])
  kmer_dat <- read.table(single_kmer_file,header=T,sep="\t")
  if (i == 1){
    kmer_dat_whole<-kmer_dat
  }
  kmer_dat_whole <- cbind(kmer_dat_whole,kmer_dat)
}

#------------------------------------------------------------------------------#
# splitting the kmer_dat_whole into tumor and normal samples

tumor_kmers <- kmer_dat_whole[,tumor_samples]
normal_kmers <- tum_norm_kmers$normal
tumor_kmers <- cbind(tumor_kmers,normal_kmers)

#------------------------------------------------------------------------------#
# filtering tumor_kmers

start<-Sys.time()
tumor_kmers_fill <- tumor_kmers
for (i in seq(nrow(tumor_kmers))[seq(1000)]){
  row<-tumor_kmers[i,]
  tumor_kmers_fill[i,]<-filter_tumor_kmers(row)
}
end<-Sys.time()
duration<-end-start
duration

a<-apply(tumor_kmers,1,filter_tumor_kmers)




