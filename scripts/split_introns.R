#!/usr/bin/env Rscript

# splitting the introns input file into a set of introns files

#------------------------------------------------------------------------------#
# loading libraries

library(optparse)

#------------------------------------------------------------------------------#
# command line input

arguments <- parse_args(OptionParser(usage = "",
               description="split introns into separate intron files",
               option_list=list(
                 make_option(c("-o","--output_directory"), default = sprintf("%s",getwd()), help="The output directory for the kmer data"),
                 make_option(c("-i","--introns"), default=NULL, help="The overall intron file"),
                 make_option(c("-s","--split_num"), default=100, help="The number of introns per split"))))

opt=arguments

out_dir<-opt$output_directory
introns_file<-opt$introns
split_num<-as.numeric(opt$split_num)

#------------------------------------------------------------------------------#
# reading in the intron file

load(introns_file)

#------------------------------------------------------------------------------#
# splitting the introns file
total_rows<-nrow(introns)
breakdown<-seq(1,total_rows,by=split_num)
logs<-lapply(seq(length(breakdown)),function(x){
  if ((breakdown[x]+split_num-1) > total_rows){
    saveRDS(introns[seq(breakdown[x],total_rows),],sprintf("%s/intron%s.rds",out_dir,x))
  } else {
    saveRDS(introns[seq(breakdown[x],breakdown[x]+split_num-1),],sprintf("%s/intron%s.rds",out_dir,x))
  }
})
