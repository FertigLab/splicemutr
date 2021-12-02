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
junc_rse <- readRDS(sprintf("%s/%s_junc_rse.rds",junc_dir,cancer))
print("extracting juncs")
junc_expr_comb <- junc_rse@assays@data@listData[["counts"]]
rm(junc_rse)

coldata <- data.frame(factor(rep("a",ncol(junc_expr_comb))))
colnames(coldata)<-"samp"
print("create_deseq")
junc_expr_comb_dds <- DESeqDataSetFromMatrix(junc_expr_comb,colData = coldata,design=~1)
print("size_fact")
junc_expr_comb_dds <- estimateSizeFactors(junc_expr_comb_dds)
print("dispersions")
junc_expr_comb_dds <- estimateDispersions(junc_expr_comb_dds,fitType="glmGamPoi")

#------------------------------------------------------------------------------#
# creating junction expression file

total <- ncol(junc_expr_comb)
iter<-1
for (i in seq(1,total,100)){
  start <- i
  end <- i + 99
  if (end > total){
    end <- total
  }
  print("sub")
  junc_expr_comb_sub<-junc_expr_comb_dds[,seq(start,end)]
  print("vst")
  junc_expr_comb_vst <- vst(junc_expr_comb_sub)
  # junc_expr_comb_vst <- as.data.frame(junc_expr_comb_vst@assays@data@listData[[1]])
  # junc_expr_comb_sub <- as.data.frame(junc_expr_comb_sub@assays@data@listData[[1]])
  saveRDS(junc_expr_comb_sub,file=sprintf("%s/junc_expr_combined_%d.rds",junc_dir,iter))
  saveRDS(junc_expr_comb_vst,file=sprintf("%s/junc_expr_combined_vst_%d.rds",junc_dir,iter))
  iter<-iter+1
}

