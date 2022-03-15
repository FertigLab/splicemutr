#!/usr/bin/env Rscript

#author: Theron Palmer
#date created: 08/11/2021


#------------------------------------------------------------------------------#
# loading libraries

library(optparse)
library(maftools)


#------------------------------------------------------------------------------#
# command line input

arguments <- parse_args(OptionParser(usage = "",
                             description="",
                             option_list=list(
                               make_option(c("-m","--maf_file"),
                                           default = sprintf("%s",getwd()),
                                           help=".maf file"),
                               make_option(c("-o","--output_dir"),
                                           default = sprintf("%s",getwd()),
                                           help="the output directory"),
                               make_option(c("-s","--target_genes_file"),
                                           default = sprintf("%s",getwd()),
                                           help="the target genes file"))))
opt=arguments

maf_file <- opt$maf_file
output_dir <- opt$output_dir
target_genes_file<- opt$target_genes_file

#------------------------------------------------------------------------------#
# summarizing missense mutations

mc3_maf = read.maf(maf_file)
target_genes <- read.table(target_genes_file,header=F)
target_genes <- as.character(target_genes[,1])

target_gene_maf = subsetMaf(mc3_maf,genes=target_genes)

#------------------------------------------------------------------------------#
# saving maf summary

write.table(target_gene_maf,
            file=sprintf("%s/maf_filt.maf",output_dir,
            sep="\t",
            quote=F,
            col.names=T,
            row.names=F)


