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
               make_option(c("-t","--tumor_dir"),
                           default = sprintf("%s",getwd()),
                           help="tumor type dir"),
               make_option(c("-c","--cancer"),
                           default = sprintf("%s",getwd()),
                           help="cancer"))))
opt=arguments
tumor_dir <- opt$tumor_dir
cancer <- opt$cancer

#------------------------------------------------------------------------------#
# internal functions

split_num <- function(vals){
  a<-data.frame(matrix(as.numeric(unlist(str_split(vals,"/"))),byrow=T,nrow=length(vals)))[,1]
  return(a)
}
split_dom <- function(vals){
  a<-data.frame(matrix(as.numeric(unlist(str_split(vals,"/"))),byrow=T,nrow=length(vals)))[,2]
  return(a)
}

#------------------------------------------------------------------------------#
# creating per_cluster score

nogos <- c("ESCA","MESO","PAAD","KIRC","GBM")

print(cancer)

splice_dat_file <- sprintf("%s/%s/%s_splicemutr_dat.txt",tumor_dir,cancer,cancer)
splice_dat <- read.table(splice_dat_file,header=T,sep="\t")
tumor_geno_file <- sprintf("%s/%s/%s_genotypes.txt",tumor_dir,cancer,cancer)
tumor_geno <- read.table(tumor_geno_file,header=T)
tumor_geno$rows <- seq(nrow(tumor_geno))
summary_file <- sprintf("%s/%s/summaries.txt",tumor_dir,cancer)
summaries <- read.table(summary_file)
summaries <- summaries$V1
summaries<-unname(vapply(summaries,function(summ){
  str_replace(summ,"kmers_summary","persamp_line")
},character(1)))
meta_file <- sprintf("%s/%s/%s_metadata.rds",tumor_dir,cancer,cancer)
meta_dat <- readRDS(meta_file)
rownames(meta_dat) <- meta_dat$external_id
#psi_file <- sprintf("%s/%s/leafcutter_run_1/data_perind.counts",tumor_dir,cancer)
psi_file <- sprintf("%s/%s/data_perind.counts",tumor_dir,cancer)

psi_dat <- read.table(psi_file,header=T,check.names=F)
if (!(any(tumor_geno$external_id %in% colnames(psi_dat)))){
  colnames(psi_dat)<-c("chrom",meta_dat[vapply(meta_dat$tcga.tcga_barcode,function(code){which(code == colnames(psi_dat)[seq(2,ncol(psi_dat))])},numeric(1)),"external_id"])
  psi_dat <- psi_dat[,c("chrom",tumor_geno$external_id)]
}


for (summ in seq(length(summaries))){
  if (summ == 1){
    summaries_combined <- read.table(summaries[summ],header=F,sep="\t")
    if (length(summaries)>1){
      summaries_combined <- summaries_combined[,seq(ncol(summaries_combined)-2)]
    }
  } else if (summ == length(summaries)){
    summaries_fill <- read.table(summaries[summ],header=F,sep="\t")
    summaries_combined <- cbind(summaries_combined,summaries_fill)
  } else {
    summaries_fill <- read.table(summaries[summ],header=F,sep="\t")
    summaries_fill <- summaries_fill[,seq(ncol(summaries_fill)-2)]
    summaries_combined <- cbind(summaries_combined,summaries_fill)
  }
}

cols_to_keep <- unname(vapply(colnames(psi_dat)[seq(2,ncol(psi_dat))],function(col_name){
  which(tumor_geno$external_id == col_name)[1]
},numeric(1)))
summaries_combined <- summaries_combined[,c(cols_to_keep,ncol(summaries_combined)-1,ncol(summaries_combined))]
tumor_geno <- tumor_geno[cols_to_keep,]
tumor_geno$rows <- seq(nrow(tumor_geno))
tumor_geno <- tumor_geno %>% dplyr::filter(type=="T")
cols_to_keep <- tumor_geno$rows

summaries_combined <- summaries_combined[,c(cols_to_keep,ncol(summaries_combined)-1,ncol(summaries_combined))]
psi_dat <- psi_dat[,c(1,cols_to_keep+1)]
colnames(summaries_combined) <- c(colnames(psi_dat)[seq(2,ncol(psi_dat))],"row","cluster")

# At this point summaries_combined, psi_dat, and tumor_geno only include tumor sample information

num <- data.frame(apply(psi_dat[,seq(2,ncol(psi_dat))],2,split_num))
denom <- data.frame(apply(psi_dat[,seq(2,ncol(psi_dat))],2,split_dom))
psi <- num/denom
is.nan.data.frame <- function(x){do.call(cbind, lapply(x, is.nan))}
psi[is.nan(psi)]<-0
psi$chrom <- psi_dat$chrom
colnames(psi)[seq(1,ncol(psi)-1)] <- colnames(psi_dat)[seq(2,ncol(psi_dat))]
rownames(psi)<-psi$chrom
splice_dat$chrom <- sprintf("%s:%s:%s:%s",splice_dat$chr,splice_dat$start,splice_dat$end,splice_dat$cluster)
psi_splice_dat <- psi[splice_dat$chrom,]
psi_summary <- psi_splice_dat[summaries_combined$row+1,]

summaries_combined_psi <- psi_summary[,seq(1,ncol(psi_summary)-1)]*summaries_combined[,seq(1,ncol(summaries_combined)-2)]

splice_dat_specific <- splice_dat[summaries_combined$row+1,]
rows_to_keep <- which(!is.na(splice_dat_specific$peptide))

summaries_combined_psi <- summaries_combined_psi[rows_to_keep,]
summaries_combined_psi[is.na(summaries_combined_psi)]<-0
splice_dat_specific <- splice_dat_specific[rows_to_keep,]
rows_to_keep <- !(splice_dat_specific$verdict == "annotated" & splice_dat_specific$modified == "changed")
summaries_combined_psi <- summaries_combined_psi[rows_to_keep,]
summaries_combined <- summaries_combined[rows_to_keep,]
splice_dat_specific <- splice_dat_specific[rows_to_keep,]
rows_to_keep <- which(splice_dat_specific$deltapsi > 0)
summaries_combined_psi <- summaries_combined_psi[rows_to_keep,]
summaries_combined <- summaries_combined[rows_to_keep,]
splice_dat_specific <- splice_dat_specific[rows_to_keep,]

clusters <- data.frame(table(splice_dat_specific$cluster))
rownames(clusters)<-clusters$Var1
clusters$Var1<-as.character(clusters$Var1)

clusters$genes <- vapply(clusters$Var1,function(clu){
  splice_dat_small <- splice_dat_specific %>% dplyr::filter(cluster == clu)
  gene<-paste(unique(splice_dat_small$gene),collapse=":")
},character(1))

splice_dat_clusters <- data.frame(t(vapply(clusters$Var1,function(clu){
  cluster_rows <- which(splice_dat_specific$cluster == clu)
  summaries_combined_small <- summaries_combined_psi[cluster_rows,]
  apply(summaries_combined_small,2,sum)/clusters[clu,"Freq"]
},numeric(ncol(summaries_combined_psi)))))

saveRDS(splice_dat_clusters,file=sprintf("%s/%s/%s_splice_dat_clusters.rds",tumor_dir,cancer,cancer))
saveRDS(clusters,file=sprintf("%s/%s/%s_clusters.rds",tumor_dir,cancer,cancer))


splice_dat_clusters_filt<-splice_dat_clusters[which(!(apply(splice_dat_clusters,1,sd)==0 & apply(splice_dat_clusters,1,sum)==0)),
                                              unname(which(apply(splice_dat_clusters,2,sum)>0))]
saveRDS(splice_dat_clusters_filt,file=sprintf("%s/%s/%s_splice_dat_clusters_filt.rds",tumor_dir,cancer,cancer))

#------------------------------------------------------------------------------#
# printing per_cluster heatmaps

pdf(file=sprintf("%s/%s/perclust_heatmaps.pdf",tumor_dir,cancer),width=10, height=10)

print(pheatmap::pheatmap(log10(splice_dat_clusters_filt+1),show_rownames=F,show_colnames=F,main="log10"))
print(pheatmap::pheatmap(splice_dat_clusters_filt,show_rownames=F,show_colnames=F,scale="row",main="zscore across row"))

dev.off()


