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
                                      make_option(c("-f","--junc_files"),
                                                   default = "",
                                                   help="junction_files"),
                                       make_option(c("-o","--out_dir"),
                                                   default = "",
                                                   help="output directory"))))
opt=arguments
junc_files <- opt$junc_files
out_dir <- opt$out_dir

#------------------------------------------------------------------------------#
# creating junction expression file

junc_expr <- list()
juncs_all<-c()

junc_files <- read.table(junc_files,header=F)

for (file in junc_files$V1){
  print(file)
  sample_juncs_counts <- read.table(file,header=F,sep="\t")
  juncs <- sprintf("%s:%s:%s:%s",sample_juncs_counts$V1,
                   sample_juncs_counts$V2,
                   sample_juncs_counts$V3,
                   sample_juncs_counts$V6)
  juncs_all <- unique(c(juncs_all,juncs))
  counts <- sample_juncs_counts$V5
  junc_expr[[str_replace(basename(file),".filt|.junc","")]]<-data.frame(juncs,counts)
}
junc_expr_comb <- data.frame(juncs_all)
rownames(junc_expr_comb)<-junc_expr_comb$juncs_all
junc_ret <- vapply(names(junc_expr),function(sample){
  sample_juncs <- junc_expr[[sample]]
  junc_expr_comb[sample_juncs$juncs,sample] <<- sample_juncs$counts
  return(T)
},logical(1))
rm(junc_expr)
junc_expr_comb[is.na(junc_expr_comb)]<-0
junc_expr_comb<-junc_expr_comb[,seq(2,ncol(junc_expr_comb))]
junc_expr_comb_vst <- as.data.frame(varianceStabilizingTransformation(as.matrix(junc_expr_comb)))

#------------------------------------------------------------------------------#
# saving junction expression file

saveRDS(junc_expr_comb,file=sprintf("%s/junc_expr_combined.rds",out_dir))
saveRDS(junc_expr_comb_vst,file=sprintf("%s/junc_expr_combined_vst.rds",out_dir))

