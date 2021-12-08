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
                                       make_option(c("-g","--gene_expression"),
                                                   default = "",
                                                   help="the gene expression file"),
                                       make_option(c("-s","--splice_dat_file"),
                                                   default = "",
                                                   help="splicemutr file"),
                                       make_option(c("-k","--kmer_counts"),
                                                   default = "",
                                                   help="the kmer counts file"),
                                       make_option(c("-j","--junc_expr_file"),
                                                   default = "",
                                                   help="junction expression file"),
                                       make_option(c("-o","--out"),
                                                   default = "",
                                                   help="output_file_prefix"))))
opt=arguments
gene_expression_file <- opt$gene_expression
splice_dat_file <- opt$splice_dat_file
kmer_counts_file <- opt$kmer_counts
junc_expr_file <- opt$junc_expr_file
out<-opt$out

#------------------------------------------------------------------------------#
# internal functions

calc_gene_expression <- function(gene_tar,gene_expression){
  gene_expr_target <- gene_expression[gene_tar,]
  if (length(gene_tar)>1){
    gene_expr_tar_min <- apply(gene_expr_target,2,min)
  } else {
    gene_expr_tar_min <- gene_expr_target
  }
  gene_expr_tar_min[is.na(gene_expr_tar_min)]<-0
  return(gene_expr_tar_min)
}
count_kmers <- function(vals){
  vapply(strsplit(vals,":"),function(val){
    length(which(val != "NAZZZZZZZ"))
  },numeric(1))
}

#------------------------------------------------------------------------------#
# local play

gene_expression_file <- "/media/theron/My_Passport/TCGA_junctions/TCGA_cancers/CHOL/gene_expression.rds"
splice_dat_file <- "/media/theron/My_Passport/TCGA_junctions/TCGA_cancers/CHOL/CHOL_splicemutr_dat.txt"
kmer_counts_file <- "/media/theron/My_Passport/TCGA_junctions/TCGA_cancers/CHOL/kmer_counts_all.rds"
vst <-1
junc_expr_file <- "/media/theron/My_Passport/TCGA_junctions/TCGA_cancers/CHOL/junc_expr_combined_vst_1.rds"
tcga<-T

#------------------------------------------------------------------------------#
# reading in the files

gene_expression <- readRDS(gene_expression_file)

if (str_detect(splice_dat_file,".txt")){
  splice_dat <- read.table(splice_dat_file,header=T,sep="\t")
} else {
  splice_dat <- readRDS(splice_dat_file)
}
if (str_detect(kmer_counts_file,".txt")){
  kmer_counts <- read.table(kmer_counts_file)
} else {
  kmer_counts <- readRDS(kmer_counts_file)
}
kmer_counts[,1]<-as.numeric(kmer_counts[,1])

kmer_counts[,seq(3,ncol(kmer_counts))] <- apply(kmer_counts[,seq(3,ncol(kmer_counts))],2,count_kmers)
splice_dat$deltapsi <- as.numeric(splice_dat$deltapsi)
splice_dat <- splice_dat %>% dplyr::filter(deltapsi>0)
splice_dat_strands <- matrix(unlist(strsplit(splice_dat$cluster,"_")),byrow=T,nrow=nrow(splice_dat))[,3]
splice_dat$juncs <- sprintf("%s:%s",splice_dat$juncs,splice_dat_strands)
splice_dat_juncs <- as.data.frame(matrix(unlist(strsplit(splice_dat$juncs,":")),byrow=T,nrow=nrow(splice_dat)))
splice_dat$juncs <- sprintf("%s:%s-%s:%s",splice_dat_juncs$V1,
                            splice_dat_juncs$V2,splice_dat_juncs$V3,
                            splice_dat_juncs$V4)
kmer_counts[,seq(3,ncol(kmer_counts))] <- mutate_all(kmer_counts[,seq(3,ncol(kmer_counts))], function(x) as.numeric(x))
splice_dat$X <- as.numeric(splice_dat$X)

junc_expr_comb <- readRDS(junc_expr_file)
junc_expr_comb <- mutate_all(junc_expr_comb, function(x) as.numeric(x))

splice_dat_filt <- splice_dat[!duplicated(splice_dat[,seq(1,ncol(splice_dat)-1)]),]
kmer_counts_filt <- kmer_counts %>% dplyr::filter(row %in% splice_dat_filt$X)
samples <- colnames(kmer_counts_filt)[seq(3,ncol(kmer_counts_filt))]
gene_expression_filt <- gene_expression[,samples]
junc_expr_comb_filt <- unique(junc_expr_comb[splice_dat_filt$juncs,samples])

rm(gene_expression)
rm(kmer_counts)
rm(splice_dat)

#------------------------------------------------------------------------------#
# calculating per-gene metric

genes <- unique(splice_dat_filt$gene)

gene_metric_mean <- as.data.frame(t(vapply(genes,function(gene_tar){
  splice_dat_small <- splice_dat_filt %>% dplyr::filter(gene==gene_tar)
  kmer_counts_small <- kmer_counts_filt %>% dplyr::filter(as.numeric(kmer_counts_filt[,1]) %in% splice_dat_small$X)
  kmer_counts <- kmer_counts_small[,samples]
  gene_split <- strsplit(gene_tar,"-")[[1]]
  gene_expr <- calc_gene_expression(gene_split,gene_expression_filt)
  junc_expr <- junc_expr_comb_filt[splice_dat_small$juncs,]
  dup_num <-nrow(splice_dat_small)
  gene_expr_dup <- as.data.frame(matrix(rep(as.numeric(gene_expr),dup_num),byrow=T,nrow=dup_num))
  a<-as.numeric(apply((kmer_counts*junc_expr)/gene_expr_dup,2,mean))
},numeric(length(samples)))))
colnames(gene_metric_mean)<-samples

gene_metric_max <- as.data.frame(t(vapply(genes,function(gene_tar){
  splice_dat_small <- splice_dat_filt %>% dplyr::filter(gene==gene_tar)
  kmer_counts_small <- kmer_counts_filt %>% dplyr::filter(as.numeric(kmer_counts_filt[,1]) %in% splice_dat_small$rows)
  kmer_counts <- kmer_counts_small[,samples]
  gene_split <- strsplit(gene_tar,"-")[[1]]
  gene_expr <- calc_gene_expression(gene_split,gene_expression_filt)
  junc_expr <- junc_expr_comb_filt[splice_dat_small$juncs,]
  dup_num <-nrow(splice_dat_small)
  gene_expr_dup <- as.data.frame(matrix(rep(as.numeric(gene_expr),dup_num),byrow=T,nrow=dup_num))
  a<-as.numeric(apply((kmer_counts*junc_expr)/gene_expr_dup,2,max))
},numeric(length(samples)))))
colnames(gene_metric_max)<-samples

#------------------------------------------------------------------------------#
# saving gene metric data

saveRDS(gene_metric_mean,file=sprintf("%s_gene_metric_mean.rds",out))
write.table(gene_metric_mean,
            file=sprintf("%s_gene_metric_mean.txt",out),quote=F, col.names = T, row.names = T, sep = "\t")

saveRDS(gene_metric_max,file=sprintf("%s_gene_metric_max.rds",out))
write.table(gene_metric_max,
            file=sprintf("%s_gene_metric_max.txt",out),quote=F, col.names = T, row.names = T, sep = "\t")
