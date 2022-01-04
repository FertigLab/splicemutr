#!/usr/bin/env Rscript

library(dplyr)
library(stringr)
library(optparse)
library(DESeq2)

#------------------------------------------------------------------------------#
# handling command line input

arguments <- parse_args(OptionParser(usage = "",
                                     description="",
                                     option_list=list(make_option(c("-f","--featurecount_files"),
                                                                  default = sprintf("%s",getwd()),
                                                                  help="featurecount_files"))))
opt=arguments
featurecount_files <- opt$featurecount_files
rep_str <- "Aligned.out.bam.sorted"

#------------------------------------------------------------------------------#
# reading in the featurecount_files

featurecount_files_file <- read.table(featurecount_files,header=F)

#------------------------------------------------------------------------------#
# compiling featurecount files


for (i in seq(nrow(featurecount_files_file))){
  print(i)
  if (i == 1){
    featurecounts_all <- read.table(featurecount_files_file[i,],header=T,sep="\t",check.names=F)
  } else {
    featurecounts <- read.table(featurecount_files_file[i,],header=T,sep="\t",check.names=F)
    featurecounts_all[,colnames(featurecounts)[ncol(featurecounts)]] <- featurecounts[,ncol(featurecounts)]
  }
}
featurecounts_all <- featurecounts_all[,c(1,seq(7,ncol(featurecounts_all)))]
colnames(featurecounts_all)[seq(2,ncol(featurecounts_all))] <- vapply(colnames(featurecounts_all)[seq(2,ncol(featurecounts_all))],function(fname){str_replace(basename(fname),rep_str,"")},character(1))
rownames(featurecounts_all)<-featurecounts_all$Geneid
featurecounts_all <- featurecounts_all[,seq(2,ncol(featurecounts_all))]
saveRDS(featurecounts_all,file=sprintf("%s/%s",dirname(featurecount_files_file[i,]),"featurecounts_all.rds"))

#------------------------------------------------------------------------------#
# variance stabilizing featurecounts

featurecounts_all_vst <- as.data.frame(varianceStabilizingTransformation(as.matrix(featurecounts_all),fitType="mean"))
colnames(featurecounts_all_vst) <- colnames(featurecounts_all)
rownames(featurecounts_all_vst) <- rownames(featurecounts_all)

saveRDS(featurecounts_all_vst,file=sprintf("%s/%s",dirname(featurecount_files_file[i,]),"featurecounts_all_vst.rds"))
