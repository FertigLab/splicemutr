#!/usr/bin/env Rscript

library(dplyr)
library(stringr)
library(optparse)
library(DESeq2)

#------------------------------------------------------------------------------#
# handling command line input

arguments <- parse_args(OptionParser(usage = "",
                                     description="",
                                     option_list=list(
                                       make_option(c("-j","--junc_rse"),
                                                   default = "",
                                                   help="junction rse file"),
                                       make_option(c("-s","--splice_dat_file"),
                                                   default = "",
                                                   help="the splicemutr data file"),
                                       make_option(c("-o","--output_dir"),
                                                   default = "",
                                                   help="output_dir"))))
opt=arguments
junc_rse <- opt$junc_rse
splice_dat_file <- opt$splice_dat_file
output_dir <- opt$output_dir

#------------------------------------------------------------------------------#
# creating and filtering junc expression

junc_rse <- readRDS(junc_rse)
juncs_for_rows <- junc_rse@rowRanges@ranges@NAMES
junc_expr_comb <- junc_rse@assays@data@listData[["counts"]]
rownames(junc_expr_comb)<-juncs_for_rows
rm(junc_rse)
junc_mat_rows <- as.data.frame(matrix(unlist(strsplit(rownames(junc_expr_comb),"[:]")),
                                      byrow=T,nrow=nrow(junc_expr_comb)))
coords <- as.data.frame(matrix(unlist(strsplit(junc_mat_rows$V2,"[-]")),nrow=nrow(junc_mat_rows),byrow=T))
junc_mat_rows$start<-as.numeric(coords$V1)
junc_mat_rows$end<-as.numeric(coords$V2)
junc_mat_rows<-junc_mat_rows[,c(1,4,5,3)]
colnames(junc_mat_rows)<-c("chr","start","end","strand")
junc_mat_rows$start <- as.numeric(junc_mat_rows$start)-1
junc_mat_rows$end <- as.numeric(junc_mat_rows$end)+1
junc_linear <- sprintf("%s:%d-%d:%s",junc_mat_rows$chr,as.numeric(junc_mat_rows$start),
                       as.numeric(junc_mat_rows$end),junc_mat_rows$strand)
rm(junc_mat_rows)
rownames(junc_expr_comb) <- junc_linear

#------------------------------------------------------------------------------#
# processing splice dat

splice_dat <- readRDS(splice_dat_file)
splice_dat$deltapsi <- as.numeric(splice_dat$deltapsi)
splice_dat <- splice_dat %>% dplyr::filter(deltapsi > 0)
strands <- as.data.frame(matrix(unlist(strsplit(splice_dat$cluster,"_")),byrow=T,nrow=nrow(splice_dat)))
strands <- as.vector(strands[,3])
juncs <- sprintf("%s:%d-%d:%s",splice_dat$chr,as.numeric(splice_dat$start),as.numeric(splice_dat$end),strands)
splice_dat$juncs <- juncs
juncs_unique <- unique(juncs)

junc_expr_comb<-junc_expr_comb[juncs_unique,]

coldata <- data.frame(factor(rep("a",ncol(junc_expr_comb))))
colnames(coldata)<-"samp"
print("create_deseq")
junc_expr_comb_dds <- DESeqDataSetFromMatrix(junc_expr_comb,colData = coldata,design=~1)
print("size_fact")
junc_expr_comb_dds <- estimateSizeFactors(junc_expr_comb_dds)
print("dispersions")
junc_expr_comb_dds <- estimateDispersions(junc_expr_comb_dds,fitType="parametric")

#------------------------------------------------------------------------------#
# creating junction expression file

total <- ncol(junc_expr_comb)
iter<-1
for (i in seq(1,total,100)){
  start <- i
  end <- i + 99
  if (end > total){
    end <- total
  }
  print("sub")
  junc_expr_comb_sub<-junc_expr_comb_dds[,seq(start,end)]
  print("vst")
  junc_expr_comb_vst <- varianceStabilizingTransformation(junc_expr_comb_sub,blind=F,fitType="parametric")
  junc_expr_comb_vst <- as.data.frame(junc_expr_comb_vst@assays@data@listData[[1]])
  junc_expr_comb_sub <- as.data.frame(junc_expr_comb_sub@assays@data@listData[[1]])
  saveRDS(junc_expr_comb_sub,file=sprintf("%s/junc_expr_combined_%d.rds",output_dir,iter))
  saveRDS(junc_expr_comb_vst,file=sprintf("%s/junc_expr_combined_vst_%d.rds",output_dir,iter))
  iter<-iter+1
}
