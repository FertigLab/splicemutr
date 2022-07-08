#!/usr/bin/env Rscript

#------------------------------------------------------------------------------#
# loading libraries

library(splicemute)
library(Biostrings)
library(stringr)

#------------------------------------------------------------------------------#
# handling command line input

arguments <- parse_args(OptionParser(usage = "%prog [options] counts_file groups_file",
               description="calculate the coding potential of splicemutr output transcripts",
               option_list=list(
                 make_option(c("-s","--splicemutr_data"), default = sprintf("%s",getwd()), help="The splicemutr data file"),
                 make_option(c("-f","--transcript_fasta"), default=NULL, help="The transcript fasta"))))



opt=arguments

splicemutr_data_file <- opt$splicemutr_data
transcript_data_file <- opt$transcript_data

#------------------------------------------------------------------------------#
# playing with internal data

splicemutr_data_file <- "F:/head_and_neck_DARIA/data/splicemutr_02_13_2022/formed_transcripts/intron1_data_splicemutr.rds"
transcript_data_file <- "F:/head_and_neck_DARIA/data/splicemutr_02_13_2022/formed_transcripts/intron1_sequences.fa"

#------------------------------------------------------------------------------#
# reading in the splicemutr data and transcript fasta

splicemutr_data <- readRDS(splicemutr_data_file)
transcript_data <- as.character(readDNAStringSet(transcript_data_file))

#------------------------------------------------------------------------------#
# calculating the coding potential

sense_codons <- unname(vapply(splicemutr_data$peptide,function(pep){return(nchar(pep)-1)},numeric(1)))
coding_potential_LGC <- vapply(seq(nrow(splicemutr_data)),function(row_val){
  calc_coding_potential_LGC(transcript_data[row_val],sense_codons[row_val])
},numeric(1))
coding_potential <- vapply(seq(nrow(splicemutr_data)),function(row_val){
  calc_coding_potential(transcript_data[row_val],sense_codons[row_val])
},numeric(1))
splicemutr_data$coding_potential <- coding_potential
splicemutr_data$coding_potential_LGC <- coding_potential_LGC

