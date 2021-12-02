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
                                       make_option(c("-j","--junc_dir"),
                                                   default = "",
                                                   help="junction_dir"))))
opt=arguments
junc_dir <- opt$junc_dir
splice_dat_file <- opt$splice_dat

#------------------------------------------------------------------------------#
# format integers

# format_int <- function(vec){
#   gsub(" ","",format(as.numeric(vec),scientific=F))
# }

#------------------------------------------------------------------------------#
# creating and filtering junc expression

cancer<-basename(junc_dir)
junc_rse <- readRDS(sprintf("%s/%s_junc_rse.rds",junc_dir,cancer))
print("extracting juncs")
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
# rm(junc_mat_rows)
rownames(junc_expr_comb) <- junc_linear

splice_dat <- read.table(sprintf("%s/%s_splicemutr_dat.txt",junc_dir,cancer),header=T,sep="\t",quote="")
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
junc_expr_comb_dds <- estimateDispersions(junc_expr_comb_dds,fitType="glmGamPoi")


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
  junc_expr_comb_vst <- vst(junc_expr_comb_sub)
  junc_expr_comb_vst <- as.data.frame(junc_expr_comb_vst@assays@data@listData[[1]])
  junc_expr_comb_sub <- as.data.frame(junc_expr_comb_sub@assays@data@listData[[1]])
  saveRDS(junc_expr_comb_sub,file=sprintf("%s/junc_expr_combined_%d.rds",junc_dir,iter))
  saveRDS(junc_expr_comb_vst,file=sprintf("%s/junc_expr_combined_vst_%d.rds",junc_dir,iter))
  iter<-iter+1
}

