#!/usr/bin/env Rscript

# The purpose of this script is to play with the coef files in order to determine the nature of the data, in ZAPPA

#------------------------------------------------------------------------------#
# loading in the libraries

library(dplyr)
library(monocle3)

#------------------------------------------------------------------------------#
# loading in a coefs files

coefs_file<-'/home/tpalme15/LUCIANE_LM_LMT_model_sc/analysis/run_12_11_2020/diff_dat/gsea_files.txt'
coefs<-read.table(coefs_file,header=F)
coefs<-readRDS(coefs[1,])

cols<-colnames(coesf)
for (col in cols){
  print(col)
  print(coef[1:10,col])
}