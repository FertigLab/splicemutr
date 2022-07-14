#!/usr/bin/env Rscript

#------------------------------------------------------------------------------#
# loading libraries

library(Biostrings)
library(stringr)
library(optparse)
# library(splicemute)

#------------------------------------------------------------------------------#
# handling command line input

arguments <- parse_args(OptionParser(usage = "%prog [options] counts_file groups_file",
               description="calculate the coding potential of splicemutr output transcripts",
               option_list=list(
                 make_option(c("-s","--splicemutr_data"), default = sprintf("%s",getwd()), help="The splicemutr data file"),
                 make_option(c("-t","--transcript_fasta"), default=NULL, help="The transcript fasta"),
                 make_option(c("-o","--out_dir"), default = sprintf("%s",getwd()), help="The output directory for the splicemutr data"),
                 make_option(c("-f","--funcs"), default = sprintf("%s",getwd()), help="The output directory for the splicemutr data"))))



opt=arguments

splicemutr_data_file <- opt$splicemutr_data
transcript_fasta_file <- opt$transcript_fasta
out_dir <- opt$out_dir
source(opt$funcs)

#------------------------------------------------------------------------------#
# playing with internal data

# splicemutr_data_file <- "F:/head_and_neck_DARIA/data/splicemutr_02_13_2022/formed_transcripts/intron1_data_splicemutr.rds"
# transcript_fasta_file <- "F:/head_and_neck_DARIA/data/splicemutr_02_13_2022/formed_transcripts/intron1_sequences.fa"

#------------------------------------------------------------------------------#
# reading in the splicemutr data and transcript fasta

if (str_detect(splicemutr_data_file,".rds")){
  splicemutr_data <- readRDS(splicemutr_data_file)
} else if (str_detect(splicemutr_data_file,".txt")){
  splicemutr_data <- read.table(splicemutr_data_file,header=T,sep="\t")
}
transcript_fasta <- as.character(readDNAStringSet(transcript_fasta_file))

#------------------------------------------------------------------------------#
# calculating the coding potential

sense_codons <- unname(vapply(splicemutr_data$peptide,function(pep){return(nchar(pep)-1)},numeric(1)))
coding_potential_LGC <- vapply(seq(nrow(splicemutr_data)),function(row_val){
  calc_coding_potential_LGC(transcript_fasta[row_val],sense_codons[row_val])
},numeric(1))
coding_potential <- vapply(seq(nrow(splicemutr_data)),function(row_val){
  calc_coding_potential(transcript_fasta[row_val],sense_codons[row_val])
},numeric(1))
splicemutr_data$coding_potential <- coding_potential
splicemutr_data$coding_potential_LGC <- coding_potential_LGC

out_txt<-sprintf("%s/%s_%s%s",out_dir,str_remove(basename(splicemutr_data_file),".rds"),"_cp_corrected",".txt")
out_rds <- sprintf("%s/%s_%s%s",out_dir,str_remove(basename(splicemutr_data_file),".rds"),"_cp_corrected",".rds")
saveRDS(splicemutr_data,file=out_rds)
write.table(splicemutr_data,file=out_txt,col.names=T,row.names=F,quote=F,sep="\t")

