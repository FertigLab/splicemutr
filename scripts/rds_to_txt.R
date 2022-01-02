#!/usr/bin/env Rscript

library(optparse)
library(stringr)

#------------------------------------------------------------------------------#
# handling command line input

arguments <- parse_args(OptionParser(usage = "",
                                     description="",
                                     option_list=list(
                                       make_option(c("-r","--rds"),
                                                   default = "",
                                                   help="the rds file"))))
opt=arguments
rds <- opt$rds
out<-str_replace(rds,"rds","txt")

#------------------------------------------------------------------------------#
# reading in and saving the data

rds_dat <- readRDS(rds)
write.table(rds_dat,file=out,col.names=T,row.names=F,quote=F,sep="\t")