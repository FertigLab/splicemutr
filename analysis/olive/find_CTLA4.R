#!/usr/bin/env Rscript

# created: 04/19/2022

# The purpose of this script is to search STAR SJ.out.tab files for annotated CTLA4 introns and retain them

#------------------------------------------------------------------------------#
# loading libraries

library(GenomicRanges)
library(rtracklayer)
library(ensembldb)
library(optparse)

#------------------------------------------------------------------------------#
# loading libraries

arguments <- parse_args(OptionParser(usage = "%prog [options] counts_file groups_file",
         description="form transcripts per junction for the given input junction file",
         option_list=list(
           make_option(c("-o","--output_directory"), default = sprintf("%s",getwd()), help="The output directory for the CTLA4 introns and counts"),
           make_option(c("-t","--txdb"), default=NULL, help="The txdb object for the annotation"),
           make_option(c("-j","--juncs"), default=NULL, help="The SJ.out.tab file"))))

opt=arguments

out_dir<-opt$output_directory
txdb_file<-opt$txdb

#------------------------------------------------------------------------------#
#internal functions

calc_strand <- function(strands){
  a<-vapply(strands,function(val){
    if (val == 0){
      return("*")
    } else if (val == 1){
      return("+")
    } else {
      return("-")
    }
  },character(1))
  return(a)
}

#------------------------------------------------------------------------------#
# preparing the reference

txdb<-loadDb(txdb_file) # making the txdb from gtf
tx_by_gene<-transcriptsBy(txdb,by="gene")
introns_by_tx<-intronsByTranscript(txdb,use.names=T)


#------------------------------------------------------------------------------#
# forming the CTLA4 intron reference

CTLA4_ensembl <- "ENSG00000163599"
CTLA4_transcripts <- tx_by_gene[[CTLA4_ensembl]]$tx_name
CTLA_4_introns_transcripts <- names(unlist(introns_by_tx[CTLA4_transcripts]))
CTLA_introns <- data.frame(unlist(introns_by_tx[CTLA4_transcripts]))
CTLA_introns$transcript <- CTLA_4_introns_transcripts
reference_juncs <- sprintf("%s:%s:%s:%s",CTLA_introns$seqnames,
                           as.numeric(CTLA_introns$start)-1,
                           as.numeric(CTLA_introns$end)+1,
                           CTLA_introns$strand)

#------------------------------------------------------------------------------#
# finding CTLA4 introns

#------------------------------------------------------------------------------#
# preparing the input junction file

juncs_file <- read.table(opt$juncs) # loading in the introns data

for (file in juncs_file[,1]){
  introns <- read.table(file)
  strands <- introns[,4]
  strands_interpreted <- calc_strand(strands) # internal function
  target_juncs <- sprintf("%s:%s:%s:%s",introns[,1],introns[,2],
                          introns[,3],strands_interpreted)
  print(which(target_juncs %in% reference_juncs))
}
