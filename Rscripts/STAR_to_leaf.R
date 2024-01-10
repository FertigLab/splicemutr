#!/usr/bin/env Rscript

# created: 05/28/2021

# converting the

#------------------------------------------------------------------------------#
# loading libraries

library(stringr)
library(optparse)
# library(splicemute)

#------------------------------------------------------------------------------#
# command line input

arguments <- parse_args(OptionParser(usage = "%prog [options] counts_file groups_file",
                 description="form transcripts per junction for the given input junction file",
                 option_list=list(
                   make_option(c("-o","--output_directory"), default = sprintf("%s",getwd()), help="The output directory for the leafcutter junc file data"),
                   make_option(c("-s","--star_file"), default=NULL, help="The star file"),
                   make_option(c("-f","--funcs"), default=NULL, help="The functions file"))))

opt=arguments

out_dir<-opt$output_directory
star_file<-opt$star_file
funcs <- opt$funcs
source(funcs)

#------------------------------------------------------------------------------#
# converting to junc file

junc_dat <- star_to_leaf(star_file)

#------------------------------------------------------------------------------#
# saving junc file

star_file_name <- basename(star_file)
star_file_name <- str_remove_all(star_file_name,"SJ.out.tab")

out<-sprintf("%s/%s.junc",out_dir,star_file_name)
write.table(junc_dat,file=out,col.names=F,quote = F,row.names = F,sep="\t")
