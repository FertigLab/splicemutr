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
junc_expr_comb <- junc_rse@assays@data@listData[["counts"]]

#------------------------------------------------------------------------------#
# creating junction expression file

junc_expr_comb_vst <- varianceStabilizingTransformation(as.matrix(junc_expr_comb))

#------------------------------------------------------------------------------#
# saving junction expression file

saveRDS(junc_expr_comb,file=sprintf("%s/junc_expr_combined.rds",junc_dir))
saveRDS(junc_expr_comb_vst,file=sprintf("%s/junc_expr_combined_vst.rds",junc_dir))

