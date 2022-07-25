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
                               make_option(c("-m","--naf_file"),
                                           default = sprintf("%s",getwd()),
                                           help=".maf file"))))
opt=arguments

maf_file <- opt$maf_file

#------------------------------------------------------------------------------#
# summarizing missense mutations

mc3_maf = read.maf(maf_file)
mc3_perSampMut = getSampleSummary(mc3_maf)


#------------------------------------------------------------------------------#
# saving maf summary

write.table(mc3_perSampMut,
            file=sprintf("%s/maf_summary.txt",dirname(maf_file)),
            sep="\t",
            quote=F,
            col.names=T,
            row.names=F)


