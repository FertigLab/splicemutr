d#!/usr/bin/env Rscript

# created: 09/28/2021

# converting the

#------------------------------------------------------------------------------#
# loading libraries

library(optparse)
library(dplyr)

#------------------------------------------------------------------------------#
# command line input

arguments <- parse_args(OptionParser(usage = "",
             description="filtering SJ.out.tab file for biologically relevant junctions",
             option_list=list(
               make_option(c("-o","--output_directory"), default = sprintf("%s",getwd()), help="The output directory for the filtered SJ.out.tab file"),
               make_option(c("-s","--sj_file"), default=NULL, help="The star file"))))

opt=arguments

out_dir<-opt$output_directory
sj_file<-opt$sj_file

#------------------------------------------------------------------------------#
# reading in the SJ.out.tab file

splice_juncs <- read.table(sj_file)
splice_juncs <- splice_juncs %>% dplyr::filter(V5 != 0)
splice_juncs <- splice_juncs %>% dplyr::filter(V7 >= 10)

#------------------------------------------------------------------------------#
# saving the biologically relevant splice junctions

out <- sprintf("%s/%s.filt",out_dir,basename(sj_file))
write.table(splice_juncs,
            file=out,
            col.names=F,
            quote = F,
            row.names = F,
            sep="\t")
