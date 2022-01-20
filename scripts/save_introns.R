#!/usr/bin/env Rscript

#author: "Theron Palmer"
#date: "08/13/2021"

#------------------------------------------------------------------------------#
# loading in the necessray libraries

library(optparse)

#------------------------------------------------------------------------------#
# handling command line input

arguments <- parse_args(OptionParser(usage = "",
               description="",
               option_list=list(
                 make_option(c("-i","--intron_file"),
                             default = sprintf("%s",getwd()),
                             help="The file containing the leafcutter intron data"),
                 make_option(c("-o","--out_prefix"),
                             default = sprintf("%s",getwd()),
                             help="The output prefix"))))
opt=arguments

intron_file <- opt$intron_file
out_prefix <- opt$out_prefix

#------------------------------------------------------------------------------#
# reading in the data.Rdata and saving the introns

load(intron_file)

write.table(introns,
            file=sprintf("%s/%s_introns.txt",dirname(intron_file),out_prefix),
            sep="\t",
            quote=F,
            col.names=T,
            row.names=F)
