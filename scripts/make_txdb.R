#!/usr/bin/env Rscript

# created: 06/02/2021

# The purpose of this script is to create and save the txdb object

#------------------------------------------------------------------------------#
# loading libraries

library(ensembldb)

#------------------------------------------------------------------------------#
# handling command line input

arguments <- parse_args(OptionParser(usage = "%prog [options] counts_file groups_file",
               description="form transcripts per junction for the given input junction file",
               option_list=list(
                 make_option(c("-o","--output_file"), default = sprintf("%s",getwd()), help="The output path and file to save the txdb as"),
                 make_option(c("-g","--gtf"), default=NULL, help="The gtf file used to create the txdb object"))))

opt=arguments

out_file<-opt$output_file
gtf_file<-opt$gtf

#------------------------------------------------------------------------------#
# creating the txdb

txdb<-makeTxDbFromGFF(gtf_file) # making the txdb from gtf

#------------------------------------------------------------------------------#
# saving the txdb

saveRDS(txdb,file=out_file)