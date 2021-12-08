#!/usr/bin/env Rscript

library(dplyr)
library(stringr)
library(optparse)
library(DESeq2)

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
# creating junc expression

cancer<-basename(junc_dir)
kmer_files <- read.table(sprintf("%s/kmer_files.txt",junc_dir))
genotypes <- read.table(sprintf("%s/%s_genotypes.txt",junc_dir,cancer),header=T)

for (i in seq(nrow(kmer_files))){
  file <- kmer_files[i,1]
  if (i == 1){
    kmer_counts_all <- read.table(str_replace(file,".txt","_filt.txt"),sep="\t")
    last_bit <- kmer_counts_all[,seq(ncol(kmer_counts_all)-1,ncol(kmer_counts_all))]
    kmer_counts_all<-kmer_counts_all[,seq(ncol(kmer_counts_all)-2)]
  } else {
    kmer_counts_fill <- read.table(file,sep="\t")
    kmer_counts_all <- cbind(kmer_counts_all,kmer_counts_fill[,seq(ncol(kmer_counts_fill)-2)])
  }
}
kmer_counts_all <- cbind(last_bit,kmer_counts_all)

colnames(kmer_counts_all)<-c("row","cluster",genotypes$external_id)

#------------------------------------------------------------------------------#
# saving junction expression file

saveRDS(kmer_counts_all,file=sprintf("%s/kmer_counts_all.rds",junc_dir))

