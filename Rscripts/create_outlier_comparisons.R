#!/usr/bin/env Rscript

library(dplyr)
library(stringr)
library(optparse)

#------------------------------------------------------------------------------#
# handling command line input

arguments <- parse_args(OptionParser(usage = "",
                 description="",
                 option_list=list(
                   make_option(c("-c","--comparison_juncs_file"),
                               default = sprintf("%s",getwd()),
                               help="comparison junction file"),
                   make_option(c("-d","--splice_dat_file"),
                               default = sprintf("%s",getwd()),
                               help="splicemutr file"),
                   make_option(c("-e","--comparisons_file"),
                               default = sprintf("%s",getwd()),
                               help="the comparisons file"),
                   make_option(c("-n","--comp_num"),
                               default = sprintf("%s",getwd()),
                               help="the comparisons number"),
                   make_option(c("-j","--junc_dir"),
                               default = sprintf("%s",getwd()),
                               help="the junctions directory"),
                   make_option(c("-o","--out_dir"),
                               default = sprintf("%s",getwd()),
                               help="the output directory"),
                   make_option(c("-l","--leaf_dir"),
                               default = sprintf("%s",getwd()),
                               help="the leafcutter dir"))))
opt=arguments
comparison_juncs_file <- opt$comparison_juncs_file
splice_dat_file <- opt$splice_dat_file
comparisons_file <- opt$comparisons_file
comp_num <- as.numeric(opt$comp_num)
junc_dir <- opt$junc_dir
out_dir <- opt$out_dir
leaf_dir <- opt$leaf_dir

#------------------------------------------------------------------------------#
# internal functions

parse_fraction <- function(frac){
  a<-as.data.frame(matrix(unlist(strsplit(frac,"/")),byrow=T,nrow=length(frac)))
  a$num <- as.numeric(a$V1)
  a$denom <- as.numeric(a$V2)
  return(a$num/a$denom)
}

parse_kmers <- function(kmers){
  a<-unlist(lapply(strsplit(as.character(kmers[seq(3,length(kmers))]),":"),function(kmer_split){
    kmer_split <- unique(kmer_split)
    length(kmer_split)-length(which(str_detect(kmer_split,"NAZZZZZZZ")))
  }))
  a<-c(as.character(kmers[seq(1,2)]),a)
  return(a)
}

#------------------------------------------------------------------------------#
# local play

# comparison_juncs_file <- "/media/theron/My_Passport/Valsamo/analysis/leafcutterMD/run_12072021/comparison_juncs.rds"
# comparisons_file <- "/media/theron/My_Passport/Valsamo/analysis/leafcutterMD/run_12072021/comparisons.rds"
# splice_dat_file <- "/media/theron/My_Passport/Valsamo/analysis/splicemutr_output/run_02052022/data_splicemutr_all_pep.txt"
# comp_num <- 1
# junc_dir <- "/media/theron/My_Passport/Valsamo/juncs"
# leaf_dir <- "/media/theron/My_Passport/Valsamo/analysis/leafcutterMD/run_12072021"

#------------------------------------------------------------------------------#
# reading in the data

comparison_juncs <- readRDS(comparison_juncs_file)
names(comparison_juncs) <- str_remove_all(names(comparison_juncs),"_pVals.txt")

splice_dat <- readRDS(splice_dat_file)
comparisons <- readRDS(comparisons_file)

splice_dat$juncs <- sprintf("%s:%s:%s:%s",splice_dat$chr,splice_dat$start,splice_dat$end,splice_dat$strand)

#------------------------------------------------------------------------------#
# creating comparison-specific junc files

comparison <- names(comparisons)[comp_num]

targets <- comparisons[[comparison]]$targets
comparators <- comparisons[[comparison]]$comparators
psi <- data.frame(seq(nrow(splice_dat)))
colnames(psi)<-"rows"
psi$juncs <- splice_dat$juncs

kmers <- data.frame(seq(nrow(splice_dat)))
colnames(kmers)<-"rows"
kmers$juncs <- splice_dat$juncs

comp_juncs <- unname(unlist(comparison_juncs))
comp_juncs_parsed <- as.data.frame(matrix(unlist(strsplit(comp_juncs,":")),
                                          byrow=T,nrow=length(comp_juncs)))
clu_split <- as.data.frame(matrix(unlist(strsplit(comp_juncs_parsed$V4,"_")),
                                  byrow=T,nrow=nrow(comp_juncs_parsed)))
strand <- clu_split$V3
comp_juncs_parsed <- sprintf("%s:%s:%s:%s",
                             comp_juncs_parsed$V1,
                             comp_juncs_parsed$V2,
                             comp_juncs_parsed$V3,
                             strand)

psi_all <- data.frame(unique(comp_juncs_parsed))
colnames(psi_all)<-"comp_juncs_parsed"
rownames(psi_all)<-psi_all$comp_juncs_parsed

for (target in targets){
  psi_file <- sprintf("%s/%s_%s_juncs_perind.counts.gz",leaf_dir,target,"baseline_POST_treatment")
  psi_dat <- read.table(psi_file,header=T)
  juncs <- as.data.frame(matrix(unlist(strsplit(psi_dat$chrom,":")),byrow=T,nrow=nrow(psi_dat)))
  strand <- as.data.frame(matrix(unlist(strsplit(juncs$V4,"_")),byrow=T,nrow=nrow(psi_dat)))
  strand <- strand$V3
  juncs <- sprintf("%s:%s:%s:%s",juncs$V1,juncs$V2,juncs$V3,strand)
  psi_dat$juncs <- juncs
  psi_dat <- psi_dat %>% dplyr::filter(juncs %in% psi_all$comp_juncs_parsed)
  psi_all[psi_dat$chrom,target] <- parse_fraction(psi_dat[,sprintf("%s.filt",str_replace_all(str_replace_all(target,"-","."),"[+]","."))])
  if (file.size(sprintf("%s/%s_splicemutr_kmers.txt",junc_dir,target))!=0){
    target_dat <- read.table(sprintf("%s/%s_splicemutr_kmers.txt",junc_dir,target),sep="\t",header=T)
    kmers[,target]<-target_dat$kmers
  } else {
    kmers[,target]<-""
  }
}
for (comp in comparators){
  if (file.size(sprintf("%s/%s_splicemutr_kmers.txt",junc_dir,comp))!=0){
    comp_dat <- read.table(sprintf("%s/%s_splicemutr_kmers.txt",junc_dir,comp),sep="\t",header=T)
    psi_all[psi_dat$chrom,sprintf("%s_comp",comp)] <- parse_fraction(psi_dat[,sprintf("%s.filt",str_replace_all(str_replace_all(comp,"-","."),"[+]","."))])
    kmers[,comp]<-comp_dat$kmers
  } else {
    kmers[,comp]<-""
  }
}
psi_all[is.na(psi_all)]<-0
splice_dat$rows <- seq(nrow(splice_dat))
splice_dat_specific <- splice_dat %>% dplyr::filter(juncs %in% comp_juncs_parsed)
splice_dat_specific <- splice_dat_specific %>% dplyr::filter(!(verdict == "annotated" & modified == "changed"))
kmers_specific <- kmers[splice_dat_specific$rows,]
kmers_specific_parsed <- kmers_specific
kmers_specific_parsed<-as.data.frame(t(apply(kmers_specific_parsed,1,parse_kmers)))
colnames(kmers_specific_parsed) <- colnames(kmers_specific)
psi_all <- psi_all[splice_dat_specific$juncs,]

#------------------------------------------------------------------------------#
# saving comparison info

saveRDS(splice_dat_specific,file=sprintf("%s/splice_dat_%s.rds",out_dir,comparison))
saveRDS(kmers_specific_parsed,file=sprintf("%s/kmers_specific_%s.rds",out_dir,comparison))
saveRDS(psi_all,file=sprintf("%s/psi_all_%s.rds",out_dir,comparison))

