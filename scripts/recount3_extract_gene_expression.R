#!/usr/bin/env Rscript

# created: 07/03/2021

# Taking the recount3 junctions and building .junc files per cancer type and tumor

#------------------------------------------------------------------------------#
# loading libraries

library(recount3)
library(stringr)
library(optparse)
library(DESeq2)

#------------------------------------------------------------------------------#
# handling command line input

arguments <- parse_args(OptionParser(usage = "%prog [options] counts_file groups_file",
                 description="extract out the TCGA cancer gene expression data",
                 option_list=list(
                   make_option(c("-c","--cancer_type"),
                               default = sprintf("%s",getwd()),
                               help="The TCGA cancer subtype"),
                   make_option(c("-o","--output_directory"),
                               default = sprintf("%s",getwd()),
                               help="Output directory"))))

opt=arguments

cancer_type <- opt$cancer_type
output_directory <- dirname(output_directory)

#------------------------------------------------------------------------------#
# downloading the rse file

gene_rse <- recount3::create_rse_manual(
  project = cancer_type,
  project_home = "data_sources/tcga",
  organism = "human",
  annotation = "gencode_v26",
  type = "gene"
)

#------------------------------------------------------------------------------#
# extracting the rse

gene_counts <- cpm(data.frame(gene_rse@assays@data@listData[["raw_counts"]],check.names = F))

saveRDS(gene_counts,file=sprintf("%s/gene_expression.rds",output_directory))

#------------------------------------------------------------------------------#
# normalizing the gene expression

gene_expression_norm <- as.data.frame(varianceStabilizingTransformation(gene_counts))

saveRDS(gene_expression_norm,file=sprintf("%s/gene_expression_vst.rds",output_directory))

