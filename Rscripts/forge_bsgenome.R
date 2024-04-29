#!/usr/bin/env Rscript

# created: 05/15/2021
# updated: 09/29/2021

# The purpose of this script is to modify reference transcripts using
# reference transcript information

# build in exons next to one another function checking if the target exons are directly next to one another in the reference transcript
# this implies that the intron modification is only taking place in a single reference transcript
# look to see if there are different

#------------------------------------------------------------------------------#
# loading libraries

library(BSgenome)
library(optparse)

#------------------------------------------------------------------------------#
# handling command line input

arguments <- parse_args(OptionParser(usage = "%prog [options] counts_file groups_file",
  description="form transcripts per junction for the given input junction file",
  option_list=list(
   make_option(c("-s","--seed_file"), default = sprintf("%s",getwd()), help="The BSgenome seed file"))))

opt=arguments

seed_file <- opt$seed_file

forgeBSgenomeDataPkg(seed_file,replace=TRUE)
