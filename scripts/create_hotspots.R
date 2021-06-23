#!/usr/bin/env Rscript

# created: 06/23/2021

# The purpose of this script is to turn the HLA binding affinity scores into peptide hostpots

#------------------------------------------------------------------------------#
# loading libraries

library(GenomicRanges)
library(stringr)
library(optparse)
library(dplyr)

#------------------------------------------------------------------------------#
# handling command line input

arguments <- parse_args(OptionParser(usage = "",
                 description="turning peptide binding affinity scores into peptide hotspot scores",
                 option_list=list(
                   make_option(c("-o","--output_directory"), default = sprintf("%s",getwd()), help="The output directory for the kmer data"),
                   make_option(c("-m","--hla_file"), default=NULL, help="The hla binding summary file"),
                   make_option(c("-s","--splicemutr_file"), default=NULL, help="The splicemutr file containing the peptides"))))

opt=arguments

out_dir<-opt$output_directory
hla_file<-opt$hla_file
splicemutr_file <- opt$splicemutr_file

#------------------------------------------------------------------------------#
# reading in the hla file and subsetting by binding rows

hla_dat <-read.table(hla_file)
splicemutr_dat <- read.table(splicemutr_file,sep="\t")
has_binders <- as.numeric(hla_dat[,1])+1
proteins <- splicemutr_dat[has_binders,14]

# start_time <- Sys.time()
hotspots<-lapply(seq(length(proteins)),function(row_val){
  row_val<<-row_val
  filler <<- rep(0,nchar(proteins[row_val]))
  peps <- unique(str_replace_all(str_split(hla_dat[row_val,2],":")[[1]],"Z",""))
  locations<<-str_locate_all(proteins[row_val],peps)
  locs<<-unlist(locations)
  locations_df <<- data.frame(matrix(locs, nrow=length(locs)/2, byrow=TRUE))
  vapply(seq(nrow(locations_df)),function(row){
    filler[seq(locations_df[row,1],locations_df[row,2])] <<- filler[seq(locations_df[row,1],locations_df[row,2])] + 1
    return(T)
  },logical(1))
  return(filler)
})
names(hotspots) <- hla_dat[seq(1000),1]
# end_time <- Sys.time()
# duration <- end_time - start_time

#------------------------------------------------------------------------------#
# saving the hotspots rds

out_file <- sprintf("%s/%s",out_dir,str_replace_all(basename(hla_file),"_tx_dict_summary.txt","_hotspots.rds"))
saveRDS(hotspots,file=out_file)
