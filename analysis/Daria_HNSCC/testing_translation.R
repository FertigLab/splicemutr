#!/usr/bin/env Rscript

#------------------------------------------------------------------------------#
# loading libraries

library(Biostrings)
library(GenomicRanges)
library(rtracklayer)
library(annotate)
library(org.Hs.eg.db)
library(ensembldb)
library(biomaRt)
library(DataCombine)
library(stringi)
library(stringr)
library(optparse)
library(dplyr)
library(BSgenome.Hsapiens.GENCODE.GRCh38.p13)
library(ggplot2)

#------------------------------------------------------------------------------#
# command line input

arguments <- parse_args(OptionParser(usage = "%prog [options] counts_file groups_file", 
                   description="form transcripts per junction for the given input junction file",
                   option_list=list(
                     make_option(c("-o","--output_directory"), default = sprintf("%s",getwd()), help="The output directory for the kmer data"),
                     make_option(c("-n","--num"), default=NULL, help="the segment to calculate"),
                     make_option(c("-g","--gtf"), default=NULL, help="The gtf file used to create the txdb object"))))

opt=arguments

out_dir<-opt$output_directory
gtf_file<-opt$gtf
num<-as.numeric(opt$num)

#------------------------------------------------------------------------------#
# forming txdb stuff

txdb<-makeTxDbFromGFF(gtf_file) # making the txdb from gtf
all_genes<-genes(txdb)
exons_by_gene<-exonsBy(txdb,by="gene")
exons_by_tx<-exonsBy(txdb,by=c("tx"),use.names=T)
tx_by_gene<-transcriptsBy(txdb,by="gene")
cds_by_tx <- cdsBy(txdb,by="tx",use.names=T)
bsgenome<-BSgenome.Hsapiens.GENCODE.GRCh38.p13

#------------------------------------------------------------------------------#
# testing translation methods

start_end <- data.frame(rep(NA,length(cds_by_tx)),rep(NA,length(cds_by_tx)),
                        rep(NA,length(cds_by_tx)),rep(NA,length(cds_by_tx)),
                        rep(NA,length(cds_by_tx)),rep(NA,length(cds_by_tx)))
colnames(start_end) <- c("tx","start","end","length_whole_tx","cds_start","cds_stop")

txs <- names(cds_by_tx)
seq_vals <- seq(1,length(txs),by=1000)
curr_txs <- txs[seq(seq_vals[num],(seq_vals[num]+999))]
curr_txs <- curr_txs[!is.na(curr_txs)]
a<-vapply(seq(length(curr_txs)),function(i){
  tx<-curr_txs[i]
  sequence<-getSeq(bsgenome,makeGRangesFromDataFrame(exons_by_tx[[tx]],keep.extra.columns=TRUE))
  sequ<-paste(as.character(sequence[1:length(sequence)]), collapse="")
  sequence_cds<-getSeq(bsgenome,makeGRangesFromDataFrame(cds_by_tx[[tx]],keep.extra.columns=TRUE))
  sequ_cds<-paste(as.character(sequence_cds[1:length(sequence_cds)]), collapse="")
  sequ_cds_len <- nchar(sequ_cds)
  start_end[i, c("tx","start","end","length_whole_tx","cds_start","cds_stop")] <<- c(tx,str_locate(sequ,sequ_cds),
                                                                       nchar(sequ),
                                                                       substr(sequ_cds,1,3),
                                                                       substr(sequ_cds,sequ_cds_len-2,sequ_cds_len))
  return(T)
},logical(1))

write.table(start_end,file=sprintf("%s/start_end_%d.txt",out_dir,num),quote=F,row.names=F,col.names=F)

#------------------------------------------------------------------------------#
# playing with reference data

start_end_file <-"/media/theron/My_Passport/testing_translation/start_end.txt"
start_end_all <- read.table(start_end_file)
colnames(start_end_all) <- c("tx","start","end","length_whole_tx","cds_start","cds_stop")
start_end_all <- start_end_all[!is.na(start_end_all$tx),]

start_codons <- "ATG"
stop_codons <- c("TAG","TAA","TGA")
valid_start <- vapply(start_end_all$cds_start,function(start_val){
  return(start_val == start_codons)
},logical(1))
valid_stop <- vapply(start_end_all$cds_stop,function(stop_val){
  return(stop_val %in% stop_codons)
},logical(1))

start_end_all$valid_start <- valid_start
start_end_all$valid_stop <- valid_stop

base_shift <- vapply(start_end_all$start,function(start_val){start_val %% 3},numeric(1))
start_end_all$base_shift <- base_shift

#------------------------------------------------------------------------------#
# looking at valid proteins

valid_start_end_all <- start_end_all %>% dplyr::filter((valid_start & valid_stop) == T)

#------------------------------------------------------------------------------#
# plotting distributions of shifts from 5' UTR

