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

cancer <- basename(tcga_junc_dir)

human_projects <- available_projects()
proj_info <- subset(
  human_projects,
  file_source %in% c("tcga")
)
rse_dat<-create_rse(proj_info[which(proj_info$project == cancer),],type="gene")
gene_expression <- rse_dat@assays@data@listData$raw_counts

#------------------------------------------------------------------------------#
# normalizing and saving the gene expression data

gene_expression_norm <- as.data.frame(varianceStabilizingTransformation(gene_expression))

saveRDS(gene_expression,file=sprintf("%s/gene_expression.rds",tcga_junc_dir))
