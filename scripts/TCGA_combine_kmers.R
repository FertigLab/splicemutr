#!/usr/bin/env Rscript

library(dplyr)
library(stringr)
library(optparse)

#------------------------------------------------------------------------------#
# handling command line input

arguments <- parse_args(OptionParser(usage = "",
                                     description="",
                                     option_list=list(
                                       make_option(c("-j","--junc_dir"),
                                                   default = "",
                                                   help="junction_dir"))))
opt=arguments
junc_dir <- opt$junc_dir

#------------------------------------------------------------------------------#
# internal functions

count_kmers <- function(vals){
  vapply(strsplit(vals,":"),function(val){
    length(which(val != "NAZZZZZZZ"))
  },numeric(1))
}

#------------------------------------------------------------------------------#
# creating junc expression

cancer<-basename(junc_dir)
kmer_files <- read.table(sprintf("%s/kmer_files.txt",junc_dir))
genotypes <- read.table(sprintf("%s/%s_genotypes.txt",junc_dir,cancer),header=T)

for (i in seq(nrow(kmer_files))){
  file <- kmer_files[i,1]
  if (i == 1){
    kmer_counts_all <- read.table(file,sep="\t",header=T)
    last_bit <- kmer_counts_all[,1,drop=F]
    kmer_counts_all<-kmer_counts_all[,seq(2,ncol(kmer_counts_all))]
    kmer_counts_all <- data.frame(apply(kmer_counts_all,2,count_kmers))
  } else {
    kmer_counts_fill <- read.table(file,sep="\t")
    kmer_counts_fill<-kmer_counts_all[,seq(2,ncol(kmer_counts_fill))]
    kmer_counts_fill <- data.frame(apply(kmer_counts_fill,2,count_kmers))
    kmer_counts_all <- cbind(kmer_counts_all,kmer_counts_fill)
  }
}
kmer_counts_all <- cbind(last_bit,kmer_counts_all)

colnames(kmer_counts_all)<-c("row",genotypes$external_id)

#------------------------------------------------------------------------------#
# saving junction expression file

saveRDS(kmer_counts_all,file=sprintf("%s/kmer_counts_all.rds",junc_dir))

