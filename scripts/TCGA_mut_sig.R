#!/usr/bin/env Rscript

library(dplyr)
library(ggplot2)
library(TCGAutils)
library(maftools)
library(stringr)
library(optparse)

#------------------------------------------------------------------------------#
# handling command line input

arguments <- parse_args(OptionParser(usage = "",
                 description="",
                 option_list=list(
                   make_option(c("-m","--maf_file"),
                               default = sprintf("%s",getwd()),
                               help="maf file"),
                   make_option(c("-s","--sbs_file"),
                               default = sprintf("%s",getwd()),
                               help="sbs file"))))
opt=arguments
maf_file <- opt$maf_file
sbs_file<-opt$sbs_file


#------------------------------------------------------------------------------#
# internal functions

calc_l1 <- function(comp,act){
  sum(abs(act-comp))
}
calc_l2<-function(comp,act){
  sum((act-comp)**2)
}



#------------------------------------------------------------------------------#
# reading in the maf file

maf_dat <- read.maf(maf_file)
sbs_dat <- read.table(sbs_file,header=T)

#------------------------------------------------------------------------------#
# prepping contents to describe mutation type

maf_mut_dat<-maf_dat@data
maf_mut_dat <- maf_mut_dat %>% dplyr::filter(Variant_Type == "SNP")
maf_mut_dat$type <- vapply(maf_mut_dat$HGVSc,function(type){
  substr(type,nchar(type)-2,nchar(type))
},character(1))
maf_mut_dat_5p_base <- vapply(maf_mut_dat$CONTEXT,function(CON){substr(CON,5,5)},character(1))
maf_mut_dat_3p_base <- vapply(maf_mut_dat$CONTEXT,function(CON){substr(CON,7,7)},character(1))
maf_mut_dat$type <- sprintf("%s[%s]%s",maf_mut_dat_5p_base,maf_mut_dat$type,maf_mut_dat_3p_base)

#------------------------------------------------------------------------------#
# annotating each mutation as apobec or not

apobec_sig <- c("T[C>G]T","T[C>G]A","T[C>T]A","T[C>T]T")
maf_mut_dat$apobec <- vapply(maf_mut_dat$type,function(type){
  if (type %in% apobec_sig){
    return("apobec")
  } else {
    return("not")
  }
},character(1))

#------------------------------------------------------------------------------#
# annotating each tumor sample with apobec percentage

mutation_types <- sbs_dat$Type
tumor_samples <- unique(as.character(maf_mut_dat$Tumor_Sample_Barcode))
percentages <- data.frame(t(vapply(tumor_samples,function(samp){
  maf_mut_dat_small <- maf_mut_dat %>% dplyr::filter(Tumor_Sample_Barcode == samp)
  vapply(mutation_types,function(type){
    length(which(maf_mut_dat_small$type == type))/nrow(maf_mut_dat_small)
  },numeric(1))
},numeric(96))))
colnames(percentages)<-mutation_types

L1_loss <- data.frame(t(vapply(seq(nrow(percentages)),function(row_val){
  tumor_vals <- as.numeric(percentages[row_val,])
  vapply(seq(2,ncol(sbs_dat)),function(col_val){
    calc_l1(tumor_vals,as.numeric(sbs_dat[,col_val]))
  },numeric(1))
},numeric(ncol(sbs_dat)-1))))
colnames(L1_loss)<-colnames(sbs_dat)[seq(2,ncol(sbs_dat))]
rownames(L1_loss)<-rownames(percentages)

L2_loss <- data.frame(t(vapply(seq(nrow(percentages)),function(row_val){
  tumor_vals <- as.numeric(percentages[row_val,])
  vapply(seq(2,ncol(sbs_dat)),function(col_val){
    calc_l2(tumor_vals,as.numeric(sbs_dat[,col_val]))
  },numeric(1))
},numeric(ncol(sbs_dat)-1))))
colnames(L2_loss)<-colnames(sbs_dat)[seq(2,ncol(sbs_dat))]
rownames(L2_loss)<-rownames(percentages)

mut_sig_type <- colnames(sbs_dat)[seq(2,ncol(sbs_dat))]
mut_sig_types_L1 <- data.frame(vapply(seq(nrow(L1_loss)),function(row_val){
  row_vals <- as.numeric(L1_loss[row_val,])
  mut_type <- mut_sig_type[which(row_vals == min(row_vals))]
  mut_type <- paste(mut_type,collapse=":")
},character(1)))
mut_sig_types$tumor_sample <- tumor_samples
colnames(mut_sig_types) <- c("SBS_TYPE","TUMOR_SAMPLE")

mut_sig_types_L2 <- data.frame(vapply(seq(nrow(L1_loss)),function(row_val){
  row_vals <- as.numeric(L2_loss[row_val,])
  mut_type <- mut_sig_type[which(row_vals == min(row_vals))]
  mut_type <- paste(mut_type,collapse=":")
},character(1)))
mut_sig_types$tumor_sample <- tumor_samples
colnames(mut_sig_types) <- c("SBS_TYPE","TUMOR_SAMPLE")

#------------------------------------------------------------------------------#
# saving the mutation type information

saveRDS(mut_sig_types_L1,file=sprintf("%s/mut_sig_types_L1.rds",dirname(maf_file)))
saveRDS(mut_sig_types_L2,file=sprintf("%s/mut_sig_types_L2.rds",dirname(maf_file)))
