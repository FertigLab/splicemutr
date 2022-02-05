#!/usr/bin/env Rscript
# testing clones

library(dplyr)
library(stringr)
library(optparse)

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
                               help="output_file_prefix"),
                   make_option(c("-t","--tcga"),
                               default = F,
                               help="tcga or not"))))
opt=arguments
gene_expression_file <- opt$gene_expression
splice_dat_file <- opt$splice_dat_file
kmer_counts_file <- opt$kmer_counts
junc_expr_file <- opt$junc_expr_file
tcga <- as.logical(opt$tcga)
out<-opt$out

#------------------------------------------------------------------------------#
# internal functions

calc_gene_expression <- function(gene_tar,gene_expression){
  gene_expr_target <- gene_expression[gene_tar,]
  gene_expr_tar_min <- as.data.frame(t(apply(gene_expr_target,2,min)))
  gene_expr_tar_min[is.na(gene_expr_tar_min)]<-0
  return(gene_expr_tar_min[1,])
}

#------------------------------------------------------------------------------#
# local play

gene_expression_file <- "/media/theron/My_Passport/head_and_neck_DARIA/data/featurecounts_out/featurecounts_all_vst.rds"
# splice_dat_file <- "/media/theron/My_Passport/head_and_neck_DARIA/data/splicemutr_05_26_2021/GENE_METRIC_01032022/full_splicemutr_dat.rds"
splice_dat_file <- "/media/theron/My_Passport/head_and_neck_DARIA/data/splicemutr_05_26_2021/GENE_METRIC_01032022/full_splicemutr_dat_norm.rds"
# kmer_counts_file <- "/media/theron/My_Passport/head_and_neck_DARIA/data/splicemutr_05_26_2021/GENE_METRIC_01032022/full_kmers_no_junc.rds"
kmer_counts_file <- "/media/theron/My_Passport/head_and_neck_DARIA/data/splicemutr_05_26_2021/GENE_METRIC_01032022/full_kmers_norm.rds"
tcga<-F
junc_expr_file <- "/media/theron/My_Passport/head_and_neck_DARIA/data/junc_expr_combined_vst.rds"
out<- "/media/theron/My_Passport/head_and_neck_DARIA/data/splicemutr_05_26_2021/GENE_METRIC_01032022/HNSCC_filt_norm"

#------------------------------------------------------------------------------#
# reading in the files

gene_expression <- readRDS(gene_expression_file)

if (str_detect(splice_dat_file,".txt")){
  splice_dat <- read.table(splice_dat_file,header=T,sep="\t")
} else {
  splice_dat <- readRDS(splice_dat_file)
}
splice_dat$juncs <- sprintf("%s:%s-%s:%s",splice_dat$chr,splice_dat$start,splice_dat$end,
                            as.data.frame(matrix(unlist(str_split(splice_dat$cluster,"_")),byrow=T,nrow=nrow(splice_dat)))[,3])
if (str_detect(kmer_counts_file,".txt")){
  kmer_counts <- read.table(kmer_counts_file)
  colnames(kmer_counts) <- c("rows",colnames(kmer_counts)[seq(2,ncol(kmer_counts))])
} else {
  kmer_counts <- readRDS(kmer_counts_file)
  colnames(kmer_counts) <- c("rows",colnames(kmer_counts)[seq(2,ncol(kmer_counts))])
}
kmer_counts$rows<-as.numeric(kmer_counts$rows)
kmer_counts[,seq(2,ncol(kmer_counts))] <- mutate_all(kmer_counts[,seq(2,ncol(kmer_counts))], function(x) as.numeric(x))
if (tcga){
  splice_dat$X <- as.numeric(splice_dat$X)
} else {
  splice_dat$rows <- as.numeric(splice_dat$rows)
}

junc_expr_comb <- readRDS(junc_expr_file)
junc_expr_comb <- mutate_all(junc_expr_comb, function(x) as.numeric(x))

splice_dat_filt <- splice_dat[!duplicated(splice_dat[,seq(1,ncol(splice_dat)-1)]),]
splicemutr_data_ann <- splice_dat_filt %>% dplyr::filter(error == "tx" & verdict == "annotated")
splicemutr_data <- splice_dat_filt %>% dplyr::filter(verdict != "annotated")
splice_dat_filt <- rbind(splicemutr_data,splicemutr_data_ann)

if (tcga){
  kmer_counts<-kmer_counts[,!duplicated(colnames(kmer_counts))]
  kmer_counts_filt <- kmer_counts %>% dplyr::filter(rows %in% splice_dat_filt$X)
} else {
  kmer_counts_filt <- kmer_counts %>% dplyr::filter(rows %in% splice_dat_filt$rows)
}
samples <- colnames(kmer_counts_filt)[seq(2,ncol(kmer_counts_filt))]
samples <- samples[which(samples %in% colnames(junc_expr_comb))]
gene_expression_filt <- gene_expression[,samples,drop=F]

strands <- as.data.frame(matrix(unlist(str_split(splice_dat_filt$cluster,"_")),byrow=T,nrow=nrow(splice_dat_filt)))[,3]
splice_dat_filt$juncs <- sprintf("%s:%s:%s:%s",splice_dat_filt$chr,splice_dat_filt$start,splice_dat_filt$end,strands)
junc_expr_comb_filt <- unique(junc_expr_comb[splice_dat_filt$juncs,samples,drop=F])

rm(gene_expression)
rm(kmer_counts)
rm(splice_dat)

#------------------------------------------------------------------------------#
# calculating per-gene metric

genes <- unique(splice_dat_filt$gene)

gene_metric_mean <- as.data.frame(t(vapply(genes,function(gene_tar){
    splice_dat_small <- splice_dat_filt %>% dplyr::filter(gene==gene_tar)
    if(tcga){
      kmer_counts_small <- kmer_counts_filt %>% dplyr::filter(rows %in% splice_dat_small$X)
    } else {
      kmer_counts_small <- kmer_counts_filt %>% dplyr::filter(rows %in% splice_dat_small$rows)
    }
    kmer_counts <- kmer_counts_small[,samples,drop=F]
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
  kmer_counts_small <- kmer_counts_filt %>% dplyr::filter(rows %in% splice_dat_small$rows)
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
