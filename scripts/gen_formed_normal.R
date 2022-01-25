#!/usr/bin/env Rscript

# created: 05/15/2021
# updated: 09/29/2021

# The purpose of this script is to modify reference transcripts using
# reference transcript information

#------------------------------------------------------------------------------#
# loading libraries

library(GenomicRanges)
library(GenomicFeatures)
library(ensembldb)

library(rtracklayer)
library(annotate)
library(org.Hs.eg.db)
library(biomaRt)
library(DataCombine)
library(stringi)
library(stringr)
library(optparse)
library(dplyr)
library(rlist)
# library(splicemute)
library(AnnotationDbi)

#------------------------------------------------------------------------------#
# handling command line input

arguments <- parse_args(OptionParser(usage = "%prog [options] counts_file groups_file",
  description="form transcripts per junction for the given input junction file",
  option_list=list(
   make_option(c("-s","--splice_dat"), default = sprintf("%s",getwd()), help="splicemutr data file"),
   make_option(c("-b","--bsgenome_name"), default = sprintf("%s",getwd()), help="the bsgenome object"),
   make_option(c("-t","--txdb_file"), default=NULL, help="The txdb object"),
   make_option(c("-f","--funcs"), default=NULL, help="The splicemute functions to source"),
   make_option(c("-m","--chr_map"),default=NA,help="the chromosome map"))))

opt=arguments

splice_dat_file<-opt$splice_dat
bsgenome_name <- opt$bsgenome_name
txdb_file <- opt$txdb_file
funcs<-opt$funcs
source(funcs)

#------------------------------------------------------------------------------#
# assigning bsgenome object to "bsgenome" variable

library(bsgenome_name,character.only = T)
assign("bsgenome",get(bsgenome_name))

#------------------------------------------------------------------------------#
# preparing the references for transcript formation and kmerization

print("reading in txdb")
txdb<-loadDb(txdb_file) # making the txdb from gtf
if (typeof(chr_map)!="logical"){
  cds_by_tx <- map_chroms(cdsBy(txdb,by="tx",use.names=T),chr_map)
}
cds_by_tx <- cdsBy(txdb,by="tx",use.names=T)

#------------------------------------------------------------------------------#
# reading in splice dat filtering for transcripts and proteins

if (str_detect(splice_dat_file,"txt")){
  splice_dat <- read.table(splice_dat_file,header=F,sep="\t")
  splice_dat <- splice_dat[,c(6,14)]
  colnames(splice_dat) <- c("transcript","protein")
} else {
  splice_dat <- readRDS(splice_dat_file)
  splice_dat <- splice_dat[,c(6,14)]
  colnames(splice_dat) <- c("transcript","protein")
}

unique_ref_txs <- unique(splice_dat$transcript)
ref_transcripts <- unlist(lapply(seq(length(unique_ref_txs)),function(tx_num){
  if (tx_num %% 1000 == 0){
    print(tx_num)
  }
  tx <- unique_ref_txs[tx_num]
  if (str_detect(tx,"-")){
    txs <- unlist(lapply(unlist(str_split(tx,"-")),function(in_tx){
      sequence<-getSeq(bsgenome,cds_by_tx[[in_tx]])
      sequ<-paste(as.character(sequence[1:length(sequence)]), collapse="")
      return(as.character(Biostrings::translate(Biostrings::DNAStringSet(sequ),no.init.codon = F)))
    }))
    paste(txs,collapse="-")
  } else {
    sequence<-getSeq(bsgenome,cds_by_tx[[tx]])
    sequ<-paste(as.character(sequence[1:length(sequence)]), collapse="")
    return(as.character(Biostrings::translate(Biostrings::DNAStringSet(sequ),no.init.codon = F)))
  }
}))

ref_txs <- data.frame(ref_transcripts)
rownames(ref_txs)<-unique_ref_txs

splice_dat$ref_transcript <- ref_txs[splice_dat$transcript,"ref_transcripts"]

#------------------------------------------------------------------------------------------------------------------------------------------------#
# saving data

out<-sprintf("%s%s%s",out_prefix,"mod_ref_transcripts",".txt")
write.table(splice_dat,file=out,row.names=F,col.names=F,quote=F,sep="\t")
