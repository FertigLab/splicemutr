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
                 make_option(c("-i","--intron_dir"),
                             default = sprintf("%s",getwd()),
                             help="The directory containing the leafcutter intron data"),
                 make_option(c("-o","--out_prefix"),
                             default = sprintf("%s",getwd()),
                             help="The output prefix"))))
opt=arguments

intron_dir <- opt$intron_dir
out_prefix <- opt$out_prefix

#------------------------------------------------------------------------------#
# reading in the data.Rdata and saving the introns

data_file <- sprintf("%s/data.Rdata",intron_dir)
load(data_file)

write.table(introns,
            file=sprintf("%s/%s_introns.txt",intron_dir,output_prefix),
            sep="\t",
            quote=F,
            col.names=T,
            row.names=F)
