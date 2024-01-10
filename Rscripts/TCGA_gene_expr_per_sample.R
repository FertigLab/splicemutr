#!/usr/bin/env Rscript

#------------------------------------------------------------------------------#
# loading libraries

library(recount3)
library(stringr)
library(optparse)
library(DESeq2)

#------------------------------------------------------------------------------#
# handling command line input

arguments <- parse_args(OptionParser(usage = "",
             description="create gene expression data per sample",
             option_list=list(
               make_option(c("-o","--output_directory"),
                           default = sprintf("%s",getwd()),
                           help="The output directory for the tcga directory"))))
opt=arguments

tcga_junc_dir <- opt$output_directory

#------------------------------------------------------------------------------#
# reading in the appropriate cancer

gene_expression <- readRDS(sprintf("%s/gene_expression.rds",tcga_junc_dir))

#------------------------------------------------------------------------------#
# normalizing and saving the gene expression data

gene_expression_norm <- as.data.frame(varianceStabilizingTransformation(gene_expression))

saveRDS(gene_expression_norm,file=sprintf("%s/gene_expression_vst.rds",tcga_junc_dir))
