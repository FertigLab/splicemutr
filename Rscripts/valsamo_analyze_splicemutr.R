#!/usr/bin/env Rscript

library(dplyr)
library(ggplot2)
library(stringr)
library(optparse)

#------------------------------------------------------------------------------#
# handling command line input

arguments <- parse_args(OptionParser(usage = "",
                 description="",
                 option_list=list(
                   make_option(c("-g","--genotypes_file"),
                               default = sprintf("%s",getwd()),
                               help="the genotypes file"),
                   make_option(c("-s","--summary_dir"),
                               default = sprintf("%s",getwd()),
                               help="summary file dir"),
                   make_option(c("-d","--splice_dat_file"),
                               default = sprintf("%s",getwd()),
                               help="splicemutr file"),
                   make_option(c("-c","--counts_file"),
                               default = sprintf("%s",getwd()),
                               help="sample junc file"),
                   make_option(c("-o","--out_dir"),
                               default = sprintf("%s",getwd()),
                               help="output directory"),
                   make_option(c("-t","--summary_type"),
                               default = "",
                               help="summary type"))))
opt=arguments
genotypes_file <- opt$genotypes_file
summary_dir <- opt$summary_dir
splice_dat_file <- opt$splice_dat_file
counts_file <- opt$counts_file
out_dir <- opt$out_dir
summary_type <- opt$summary_type

#------------------------------------------------------------------------------#
# local directories and file inputs for testing

# genotypes_file <- "/Volumes/One_Touch/Valsamo/genotypes/genotypes.rds"
# summary_dir <- "/Volumes/One_Touch/Valsamo/mhcnuggets_out/run_12072021/predictions_1"
# splice_dat_file <- "/Volumes/One_Touch/Valsamo/analysis/splicemutr_output/run_12072021/data_splicemutr.txt"
# counts_file <- "/Volumes/One_Touch/Valsamo/juncs/Q21777-Plate-1-A01_L15.filt.junc"

#------------------------------------------------------------------------------#
# reading in the data necessary for creating specific splicemutr data

genotypes <- readRDS(genotypes_file)
splice_dat <- readRDS(splice_dat_file)
counts <- read.table(counts_file)

#------------------------------------------------------------------------------#
# formatting counts file

counts$V2 <- as.numeric(counts$V2)
counts$V3 <- as.numeric(counts$V3)
counts$V5 <- as.numeric(counts$V5)
counts$juncs <- sprintf("%s:%s:%s:%s",counts$V1,counts$V2,counts$V3,counts$V6)

#------------------------------------------------------------------------------#
# reading in the genotype for the specific sample

sample <- basename(str_replace(counts_file,".filt.junc",""))
sample_geno <- unname(genotypes[[which(names(genotypes) == sample)]])
class_1 <- c("HLA-A","HLA-B","HLA-C")
sample_geno <- str_replace(sample_geno[which(str_detect(sample_geno,class_1[1]) |
                                               str_detect(sample_geno,class_1[2]) |
                                               str_detect(sample_geno,class_1[3]))], ":","-")

sample_kmers <- data.frame(seq(nrow(splice_dat)))
colnames(sample_kmers) <- "rows"
rownames(sample_kmers) <- sample_kmers$rows
sample_kmers$kmers <- NA

sample_kmers_ret <- vapply(seq(length(sample_geno)),function(geno_val){
  HLA<-sample_geno[geno_val]
  if (summary_type == "perc"){
    HLA_summ_file <- sprintf("%s/%s_tx_dict_summary_perc.txt",summary_dir,HLA,summary_type)
  } else {
    HLA_summ_file <- sprintf("%s/%s_tx_dict_summary.txt",summary_dir,HLA,summary_type)
  }
  if (file.size(HLA_summ_file)!=0){
    HLA_summ <- read.table(HLA_summ_file,header=F,sep="\t")
    HLA_summ$V1 <- as.numeric(HLA_summ$V1)+1
    sample_kmers[HLA_summ$V1,"kmers"] <<- vapply(seq(nrow(HLA_summ)),function(row_val){
      dat_row_val <- HLA_summ[row_val,"V1"]
      if (is.na(sample_kmers$kmers[dat_row_val])){
        HLA_summ[row_val,"V2"]
      } else {
        paste(c(sample_kmers$kmers[dat_row_val],HLA_summ[row_val,"V2"]),collapse=":")
      }
    },character(1))
  }
  return(TRUE)
},logical(1))

#------------------------------------------------------------------------------#
# annotating sample_kmers with junctions

juncs <- sprintf("%s:%s:%s:%s",splice_dat$chr,splice_dat$start,splice_dat$end,splice_dat$strand)
splice_dat$juncs <- juncs
sample_kmers$juncs <- juncs
rownames(counts)<-counts$juncs
sample_kmers$counts <- counts[sample_kmers$juncs,"V5"]

#------------------------------------------------------------------------------#
# saving samples_kmers file

sample_kmers_file <- sprintf("%s/%s_splicemutr_kmers.rds",out_dir,sample)
saveRDS(sample_kmers,file=sample_kmers_file)
