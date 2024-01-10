#!/usr/bin/env Rscript

# created: 07/03/2021

# Taking the recount3 junctions and building .junc files per cancer type and tumor

#------------------------------------------------------------------------------#
# loading libraries

library(recount3)
library(stringr)
library(optparse)

#------------------------------------------------------------------------------#
# handling command line input

arguments <- parse_args(OptionParser(usage = "%prog [options] counts_file groups_file",
                 description="form transcripts per junction for the given input junction file",
                 option_list=list(
                   make_option(c("-g","--gene_expression"),
                               default = sprintf("%s",getwd()),
                               help="The raw gene expression file"))))

opt=arguments

gene_expression_file <- opt$gene_expression
output_directory <- dirname(gene_expression_file)

#------------------------------------------------------------------------------#
# reading in the raw gene expression data

gene_expression <- log2(data.frame(readRDS(gene_expression_file),check.names=F,check.rows=F)+1)

#------------------------------------------------------------------------------#
# reading data in from each tcga project



saveRDS(junc_rse,file=sprintf("%s/%s/junc_rse.rds",tcga_junc_dir,cancer))
saveRDS(junc_metadata,file=sprintf("%s/%s/junc_metadata.rds",tcga_junc_dir,cancer))


all_cancer_dirs <-data.frame(all_cancer_dirs)
colnames(all_cancer_dirs)<-"cancer_dir"
write.table(all_cancer_dirs,file=sprintf("%s/cancer_dirs.txt",tcga_junc_dir),col.names = F,row.names = F,quote = F)


