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
                   make_option(c("-o","--outlier_file"),
                               default = sprintf("%s",getwd()),
                               help="the genotypes file"),
                   make_option(c("-s","--significance_level"),
                               default = sprintf("%s",getwd()),
                               help="the significance level"))))
opt=arguments
outlier_file <- opt$outlier_file
significance_level <- as.numeric(opt$significance_level)

#------------------------------------------------------------------------------#
# creating the leafcutter comparison juncs list file

comparison_junctions <- list()
outlier_juncs <- read.table(outlier_file,check.names = F)
outlier_junc_cols <- str_replace_all(colnames(outlier_juncs),".filt","")
colnames(outlier_juncs) <- str_replace_all(outlier_junc_cols,"[-_+.]","_")
file_split <- strsplit(basename(outlier_file),"_")[[1]]
sample <- basename(file_split[1])
sig_juncs <- rownames(outlier_juncs)[as.numeric(outlier_juncs[,1]) <= significance_level]
comparison_junctions[[sample]] <- unique(c(comparison_junctions[[sample]],sig_juncs))

#------------------------------------------------------------------------------#
# creating the introns object from the pvalue data

comparison_juncs_linear <- data.frame(t(vapply(unname(unlist(comparison_junctions)),function(junc){
  junc_vals<-strsplit(junc,":")[[1]]
  strand<-strsplit(junc_vals[4],"_")[[1]][3]
  c(junc_vals[1],junc_vals[2],junc_vals[3],strand)
},character(4))))
rownames(comparison_juncs_linear)<-seq(nrow(comparison_juncs_linear))
colnames(comparison_juncs_linear)<-c("chr","start","end","strand")
comparison_juncs_linear<-unique(comparison_juncs_linear)

#------------------------------------------------------------------------------#
# saving the introns object for input into form_transcripts.R and the comparison juncs file

saveRDS(comparison_juncs_linear,
        sprintf("%s/comparison_juncs_linear.rds",dirname(outlier_file)))
saveRDS(comparison_junctions,file=sprintf("%s/comparison_juncs.rds",dirname(outlier_file)))
