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
                   # make_option(c("-g","--gene_expression"),
                   #             default = "",
                   #             help="the gene expression file"),
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
# gene_expression_file <- opt$gene_expression
splice_dat_file <- opt$splice_dat_file
kmer_counts_file <- opt$kmer_counts
junc_expr_file <- opt$junc_expr_file
out<-opt$out

#------------------------------------------------------------------------------#
# internal functions

# calc_gene_expression <- function(gene_tar,gene_expression){
#   gene_expr_target <- gene_expression[gene_tar,,drop=F]
#   gene_expr_tar_min <- as.data.frame(t(apply(gene_expr_target,2,min)))
#   gene_expr_tar_min[is.na(gene_expr_tar_min)]<-0
#   return(gene_expr_tar_min[1,])
# }

calc_kmers <- function(peptides){
  K<-9
  peptides <- vapply(peptides,function(pep){substr(pep,1,nchar(pep)-1)},character(1))
  peptides<-lapply(peptides,function(pep){
    if (nchar(pep)>K){
      stringi::stri_sub(str = pep,from = seq(1, nchar(pep) - K + 1,by = 1),length = K)
    } else {
      return(1)
    }})
  return(length(unique(unlist(peptides))))
}

format_juncs <- function(juncs){
  junc_split <- strsplit(juncs[1],":")[[1]]
  if (length(junc_split)==3){
    return(unlist(lapply(juncs,function(junc){
      a<-unlist(strsplit(junc,":"))
      a[2]<-str_replace(a[2],"-",":")
      paste(a,collapse=":")
    })))
  } else {
    return(juncs)
  }
}

#------------------------------------------------------------------------------#
# local testing
# # 
# splice_dat_file <- "/Users/tpalme15/Desktop/splicemutr_paper/simulated_reads/data_splicemutr_all_pep.txt"
# kmer_counts_file <- "/Users/tpalme15/Desktop/splicemutr_paper/simulated_reads/all_kmers_counts.txt"
# junc_expr_file <- "/Users/tpalme15/Desktop/splicemutr_paper/simulated_reads/junc_expr_combined_vst.rds"

# splice_dat_file <- "/Users/tpalme15/Desktop/splicemutr_paper/melanoma_cohort/input_data/splice_dat_baseline_vs_post_treatment.rds"
# kmer_counts_file <- "/Users/tpalme15/Desktop/splicemutr_paper/melanoma_cohort/input_data/kmers_specific_baseline_vs_post_treatment.rds"
# junc_expr_file <- "/Users/tpalme15/Desktop/splicemutr_paper/simulated_reads/junc_expr_combined_vst.rds"

# splice_dat_file <- "/Users/tpalme15/Desktop/splicemutr_paper/melanoma_cohort/input_data/splice_dat_data.rds"
# kmer_counts_file <- "/Users/tpalme15/Desktop/splicemutr_paper/melanoma_cohort/input_data/kmers_specific_data.rds"
# junc_expr_file <- "/Users/tpalme15/Desktop/splicemutr_paper/melanoma_cohort/input_data/junc_expr_combined_vst.rds"

#------------------------------------------------------------------------------#
# reading in the files

# if (str_detect(gene_expression_file,"vst")){
#   gene_expression <- readRDS(gene_expression_file)
# } else {
#   gene_expression <- data.frame(readRDS(gene_expression_file),check.names=F,check.rows=F)
# }

if (str_detect(splice_dat_file,".txt")){
  splice_dat <- read.table(splice_dat_file,header=T,sep="\t")
  splice_dat$juncs <- format_juncs(splice_dat$juncs)
} else {
  splice_dat <- readRDS(splice_dat_file)
  splice_dat$juncs <- format_juncs(splice_dat$juncs)
}
if (str_detect(kmer_counts_file,".txt")){
  kmer_counts <- read.table(kmer_counts_file,header=T,check.names=F)
  kmer_counts$rows <- seq(nrow(kmer_counts))
} else {
  kmer_counts <- readRDS(kmer_counts_file)
  kmer_counts$rows <- seq(nrow(kmer_counts))
  samples <- colnames(kmer_counts)[seq(3,ncol(kmer_counts))]
  kmer_counts <- kmer_counts[,c(samples,"rows")]
}
kmer_counts<- mutate_all(kmer_counts, function(x) as.numeric(x))

splice_dat$rows <- as.numeric(seq(nrow(splice_dat)))


junc_expr_comb <- readRDS(junc_expr_file)
junc_expr_comb <- mutate_all(junc_expr_comb, function(x) as.numeric(x))
colnames(junc_expr_comb) <- str_remove(colnames(junc_expr_comb),".junc")
rownames(junc_expr_comb) <- format_juncs(rownames(junc_expr_comb))

splice_dat_filt <- splice_dat[!duplicated(splice_dat[,seq(1,ncol(splice_dat)-5)]),]
splice_dat_filt <- splice_dat_filt %>% dplyr::filter(nchar(peptide)-1>50)

kmer_counts_filt <- kmer_counts[splice_dat_filt$rows,]

samples <- colnames(kmer_counts_filt)

if (str_detect(samples[1],".filt")){
  samples<-str_remove_all(samples,".filt")
}
samples <- samples[which(samples %in% colnames(junc_expr_comb))]
colnames(kmer_counts_filt) <- c(samples,"rows")
# gene_expression_filt <- gene_expression[,samples,drop=F]

splice_dat_filt <- splice_dat_filt[splice_dat_filt$juncs %in% rownames(junc_expr_comb),]
junc_expr_comb_filt <- unique(junc_expr_comb[splice_dat_filt$juncs,samples,drop=F])

# rm(gene_expression)
rm(kmer_counts)
rm(splice_dat)

#------------------------------------------------------------------------------#
# calculating per-gene metric

if (all(is.na(splice_dat_filt$deltapsi))){
  genes <- unique(splice_dat_filt$gene)
  # gene_metric_mean <- as.data.frame(t(vapply(genes,function(gene_tar){
  #   g<<-gene_tar
  #   splice_dat_small <- splice_dat_filt %>% dplyr::filter(gene==gene_tar)
  #   splice_dat_small$junc_kmers <- vapply(splice_dat_small$peptide,calc_kmers,numeric(1))
  #   splice_dat_small$gene_kmers <- calc_kmers(splice_dat_small$peptide)
  #   kmer_counts_small <- kmer_counts_filt[vapply(splice_dat_small$rows,function(val){which(kmer_counts_filt$rows == val)},numeric(1)),]
  #   kmer_counts <- kmer_counts_small[,samples,drop=F]/splice_dat_small$junc_kmers
  #   gene_split <- strsplit(gene_tar,"-")[[1]]
  #   gene_expr <- calc_gene_expression(gene_split,gene_expression_filt)
  #   junc_expr <- junc_expr_comb_filt[splice_dat_small$juncs,]
  #   dup_num <-nrow(splice_dat_small)
  #   gene_expr_dup <- as.data.frame(matrix(rep(as.numeric(gene_expr),dup_num),byrow=T,nrow=dup_num))
  #   junc_expr[junc_expr<0]<-NA
  #   gene_expr_dup[gene_expr_dup<0]<-NA
  #   a<-as.numeric(apply((kmer_counts*junc_expr)/gene_expr_dup,2,mean,na.rm=T))
  # },numeric(length(samples)))))
  # colnames(gene_metric_mean)<-samples
  #
  # saveRDS(gene_metric_mean,file=sprintf("%s_gene_metric_mean_len_norm.rds",out))
  # write.table(gene_metric_mean,
  #             file=sprintf("%s_gene_metric_mean_len_norm.txt",out),quote=F, col.names = T, row.names = T, sep = "\t")

  if (length(samples)>1){
    gene_metric_mean_no_gene_norm <- as.data.frame(t(vapply(genes,function(gene_tar){
      g<<-gene_tar
      splice_dat_small <- splice_dat_filt %>% dplyr::filter(gene==gene_tar)
      splice_dat_small$junc_kmers <- vapply(splice_dat_small$peptide,calc_kmers,numeric(1))
      splice_dat_small$gene_kmers <- calc_kmers(splice_dat_small$peptide)
      kmer_counts_small <- kmer_counts_filt[vapply(splice_dat_small$rows,function(val){which(kmer_counts_filt$rows == val)},numeric(1)),]
      kmer_counts <- kmer_counts_small[,samples,drop=F]/splice_dat_small$junc_kmers
      gene_split <- strsplit(gene_tar,"-")[[1]]
      # gene_expr <- calc_gene_expression(gene_split,gene_expression_filt)
      junc_expr <- junc_expr_comb_filt[splice_dat_small$juncs,]
      dup_num <-nrow(splice_dat_small)
      # gene_expr_dup <- as.data.frame(matrix(rep(as.numeric(gene_expr),dup_num),byrow=T,nrow=dup_num))
      junc_expr[junc_expr<0]<-NA
      # gene_expr_dup[gene_expr_dup<0]<-NA
      a<-as.numeric(apply((kmer_counts*junc_expr),2,mean,na.rm=T))
    },numeric(length(samples)))))
    colnames(gene_metric_mean_no_gene_norm)<-samples
  } else {
    gene_metric_mean_no_gene_norm <- as.data.frame((vapply(genes,function(gene_tar){
      g<<-gene_tar
      splice_dat_small <- splice_dat_filt %>% dplyr::filter(gene==gene_tar)
      splice_dat_small$junc_kmers <- vapply(splice_dat_small$peptide,calc_kmers,numeric(1))
      splice_dat_small$gene_kmers <- calc_kmers(splice_dat_small$peptide)
      kmer_counts_small <- kmer_counts_filt[vapply(splice_dat_small$rows,function(val){which(kmer_counts_filt$rows == val)},numeric(1)),]
      kmer_counts <- kmer_counts_small[,samples,drop=F]/splice_dat_small$junc_kmers
      gene_split <- strsplit(gene_tar,"-")[[1]]
      # gene_expr <- calc_gene_expression(gene_split,gene_expression_filt)
      junc_expr <- junc_expr_comb_filt[splice_dat_small$juncs,]
      dup_num <-nrow(splice_dat_small)
      # gene_expr_dup <- as.data.frame(matrix(rep(as.numeric(gene_expr),dup_num),byrow=T,nrow=dup_num))
      junc_expr[junc_expr<0]<-NA
      # gene_expr_dup[gene_expr_dup<0]<-NA
      a<-as.numeric(apply((kmer_counts*junc_expr),2,mean,na.rm=T))
    },numeric(length(samples)))))
    colnames(gene_metric_mean_no_gene_norm)<-samples
  }

  saveRDS(gene_metric_mean_no_gene_norm,file=sprintf("%s_gene_metric_mean_len_norm_no_gene_norm.rds",out))
  write.table(gene_metric_mean_no_gene_norm,
              file=sprintf("%s_gene_metric_mean_len_norm_no_gene_norm.txt",out),quote=F, col.names = T, row.names = T, sep = "\t")

  coding_potential_LGC <- as.data.frame(vapply(genes,function(gene_tar){
    g<<-gene_tar
    splice_dat_small <- splice_dat_filt %>% dplyr::filter(gene==gene_tar)
    return(mean(splice_dat_small$coding_potential_LGC,na.rm=T))
  },numeric(1)))
  colnames(coding_potential_LGC)<-"coding_potential"

  saveRDS(coding_potential_LGC,file=sprintf("%s_coding_potential_LGC.rds",out))
  write.table(coding_potential_LGC,
              file=sprintf("%s_coding_potential_LGC.txt",out),quote=F,col.names=T,row.names=T,sep="\t")

  coding_potential <- as.data.frame(vapply(genes,function(gene_tar){
    g<<-gene_tar
    splice_dat_small <- splice_dat_filt %>% dplyr::filter(gene==gene_tar)
    return(mean(splice_dat_small$coding_potential,na.rm=T))
  },numeric(1)))
  colnames(coding_potential)<-"coding_potential"

  saveRDS(coding_potential,file=sprintf("%s_coding_potential.rds",out))
  write.table(coding_potential,
              file=sprintf("%s_coding_potential.txt",out),quote=F,col.names=T,row.names=T,sep="\t")

} else {
  splice_dat_filt_normal <- splice_dat_filt[as.numeric(splice_dat_filt$deltapsi)<0,]
  splice_dat_filt_tumor <- splice_dat_filt[as.numeric(splice_dat_filt$deltapsi)>0,]
  genes_normal <- unique(splice_dat_filt_normal$gene)
  genes_tumor <- unique(splice_dat_filt_tumor$gene)

  # gene_metric_mean_normal <- as.data.frame(t(vapply(genes_normal,function(gene_tar){
  #     g<<-gene_tar
  #     splice_dat_small <- splice_dat_filt_normal %>% dplyr::filter(gene==gene_tar)
  #     splice_dat_small$junc_kmers <- vapply(splice_dat_small$peptide,calc_kmers,numeric(1))
  #     splice_dat_small$gene_kmers <- calc_kmers(splice_dat_small$peptide)
  #     kmer_counts_small <- kmer_counts_filt[vapply(splice_dat_small$rows,function(val){which(kmer_counts_filt$rows == val)},numeric(1)),]
  #     kmer_counts <- kmer_counts_small[,samples,drop=F]/splice_dat_small$junc_kmers
  #     gene_split <- strsplit(gene_tar,"-")[[1]]
  #     gene_expr <- calc_gene_expression(gene_split,gene_expression_filt)
  #     junc_expr <- junc_expr_comb_filt[splice_dat_small$juncs,]
  #     dup_num <-nrow(splice_dat_small)
  #     gene_expr_dup <- as.data.frame(matrix(rep(as.numeric(gene_expr),dup_num),byrow=T,nrow=dup_num))
  #     junc_expr[junc_expr<0]<-NA
  #     gene_expr_dup[gene_expr_dup<0]<-NA
  #     a<-as.numeric(apply((kmer_counts*junc_expr)/gene_expr_dup,2,mean,na.rm=T))
  # },numeric(length(samples)))))
  # colnames(gene_metric_mean_normal)<-samples
  #
  # saveRDS(gene_metric_mean_normal,file=sprintf("%s_gene_metric_mean_len_norm_normal.rds",out))
  # write.table(gene_metric_mean_normal,
  #             file=sprintf("%s_gene_metric_mean_len_norm_normal.txt",out),quote=F, col.names = T, row.names = T, sep = "\t")
  # rm(gene_metric_mean_normal)
  
  if (length(samples)>1){
    gene_metric_mean_normal_no_gene_norm <- as.data.frame(t(vapply(genes_normal,function(gene_tar){
      g<<-gene_tar
      splice_dat_small <- splice_dat_filt_normal %>% dplyr::filter(gene==gene_tar)
      splice_dat_small$junc_kmers <- vapply(splice_dat_small$peptide,calc_kmers,numeric(1))
      splice_dat_small$gene_kmers <- calc_kmers(splice_dat_small$peptide)
      kmer_counts_small <- kmer_counts_filt[vapply(splice_dat_small$rows,function(val){which(kmer_counts_filt$rows == val)},numeric(1)),]
      kmer_counts <- kmer_counts_small[,samples,drop=F]/splice_dat_small$junc_kmers
      gene_split <- strsplit(gene_tar,"-")[[1]]
      # gene_expr <- calc_gene_expression(gene_split,gene_expression_filt)
      junc_expr <- junc_expr_comb_filt[splice_dat_small$juncs,]
      dup_num <-nrow(splice_dat_small)
      junc_expr[junc_expr<0]<-NA
      a<-as.numeric(apply((kmer_counts*junc_expr),2,mean,na.rm=T))
    },numeric(length(samples)))))
    colnames(gene_metric_mean_normal_no_gene_norm)<-samples
  } else {
    gene_metric_mean_normal_no_gene_norm <- as.data.frame((vapply(genes_normal,function(gene_tar){
      g<<-gene_tar
      splice_dat_small <- splice_dat_filt_normal %>% dplyr::filter(gene==gene_tar)
      splice_dat_small$junc_kmers <- vapply(splice_dat_small$peptide,calc_kmers,numeric(1))
      splice_dat_small$gene_kmers <- calc_kmers(splice_dat_small$peptide)
      kmer_counts_small <- kmer_counts_filt[vapply(splice_dat_small$rows,function(val){which(kmer_counts_filt$rows == val)},numeric(1)),]
      kmer_counts <- kmer_counts_small[,samples,drop=F]/splice_dat_small$junc_kmers
      gene_split <- strsplit(gene_tar,"-")[[1]]
      # gene_expr <- calc_gene_expression(gene_split,gene_expression_filt)
      junc_expr <- junc_expr_comb_filt[splice_dat_small$juncs,]
      dup_num <-nrow(splice_dat_small)
      junc_expr[junc_expr<0]<-NA
      a<-as.numeric(apply((kmer_counts*junc_expr),2,mean,na.rm=T))
    },numeric(length(samples)))))
    colnames(gene_metric_mean_normal_no_gene_norm)<-samples    
  }

  saveRDS(gene_metric_mean_normal_no_gene_norm,file=sprintf("%s_gene_metric_mean_len_norm_no_gene_norm_normal.rds",out))
  write.table(gene_metric_mean_normal_no_gene_norm,
              file=sprintf("%s_gene_metric_mean_len_norm_no_gene_norm_normal.txt",out),quote=F, col.names = T, row.names = T, sep = "\t")
  rm(gene_metric_mean_normal_no_gene_norm)

  coding_potential_LGC_normal <- as.data.frame(vapply(genes_normal,function(gene_tar){
    g<<-gene_tar
    splice_dat_small <- splice_dat_filt_normal %>% dplyr::filter(gene==gene_tar)
    return(mean(splice_dat_small$coding_potential_LGC,na.rm=T))
  },numeric(1)))
  colnames(coding_potential_LGC_normal)<-"coding_potential_LGC"

  saveRDS(coding_potential_LGC_normal,file=sprintf("%s_coding_potential_LGC_normal.rds",out))
  write.table(coding_potential_LGC_normal,
              file=sprintf("%s_coding_potential_LGC_normal.txt",out),quote=F,col.names=T,row.names=T,sep="\t")

  coding_potential_normal <- as.data.frame(vapply(genes_normal,function(gene_tar){
    g<<-gene_tar
    splice_dat_small <- splice_dat_filt_normal %>% dplyr::filter(gene==gene_tar)
    return(mean(splice_dat_small$coding_potential,na.rm=T))
  },numeric(1)))
  colnames(coding_potential_normal)<-"coding_potential"

  saveRDS(coding_potential_normal,file=sprintf("%s_coding_potential_normal.rds",out))
  write.table(coding_potential_normal,
              file=sprintf("%s_coding_potential_normal.txt",out),quote=F,col.names=T,row.names=T,sep="\t")

  # gene_metric_mean_tumor <- as.data.frame(t(vapply(genes_tumor,function(gene_tar){
  #   g<<-gene_tar
  #   splice_dat_small <- splice_dat_filt_tumor %>% dplyr::filter(gene==gene_tar)
  #   splice_dat_small$junc_kmers <- vapply(splice_dat_small$peptide,calc_kmers,numeric(1))
  #   splice_dat_small$gene_kmers <- calc_kmers(splice_dat_small$peptide)
  #   kmer_counts_small <- kmer_counts_filt[vapply(splice_dat_small$rows,function(val){which(kmer_counts_filt$rows == val)},numeric(1)),]
  #   kmer_counts <- kmer_counts_small[,samples,drop=F]/splice_dat_small$junc_kmers
  #   gene_split <- strsplit(gene_tar,"-")[[1]]
  #   gene_expr <- calc_gene_expression(gene_split,gene_expression_filt)
  #   junc_expr <- junc_expr_comb_filt[splice_dat_small$juncs,]
  #   dup_num <-nrow(splice_dat_small)
  #   gene_expr_dup <- as.data.frame(matrix(rep(as.numeric(gene_expr),dup_num),byrow=T,nrow=dup_num))
  #   junc_expr[junc_expr<0]<-NA
  #   gene_expr_dup[gene_expr_dup<0]<-NA
  #   a<-as.numeric(apply((kmer_counts*junc_expr)/gene_expr_dup,2,mean,na.rm=T))
  # },numeric(length(samples)))))
  # colnames(gene_metric_mean_tumor)<-samples
  #
  # saveRDS(gene_metric_mean_tumor,file=sprintf("%s_gene_metric_mean_len_norm_tumor.rds",out))
  # write.table(gene_metric_mean_tumor,
  #             file=sprintf("%s_gene_metric_mean_len_norm_tumor.txt",out),quote=F, col.names = T, row.names = T, sep = "\t")
  if (length(samples)>1){
    gene_metric_mean_tumor_no_gene_norm <- as.data.frame((vapply(genes_tumor,function(gene_tar){
      g<<-gene_tar
      splice_dat_small <- splice_dat_filt_tumor %>% dplyr::filter(gene==gene_tar)
      splice_dat_small$junc_kmers <- vapply(splice_dat_small$peptide,calc_kmers,numeric(1))
      splice_dat_small$gene_kmers <- calc_kmers(splice_dat_small$peptide)
      kmer_counts_small <- kmer_counts_filt[vapply(splice_dat_small$rows,function(val){which(kmer_counts_filt$rows == val)},numeric(1)),]
      kmer_counts <- kmer_counts_small[,samples,drop=F]/splice_dat_small$junc_kmers
      gene_split <- strsplit(gene_tar,"-")[[1]]
      # gene_expr <- calc_gene_expression(gene_split,gene_expression_filt)
      junc_expr <- junc_expr_comb_filt[splice_dat_small$juncs,]
      dup_num <-nrow(splice_dat_small)
      junc_expr[junc_expr<0]<-NA
      a<-as.numeric(apply((kmer_counts*junc_expr),2,mean,na.rm=T))
    },numeric(length(samples)))))
    colnames(gene_metric_mean_tumor_no_gene_norm)<-samples
  } else {
    gene_metric_mean_tumor_no_gene_norm <- as.data.frame(t(vapply(genes_tumor,function(gene_tar){
      g<<-gene_tar
      splice_dat_small <- splice_dat_filt_tumor %>% dplyr::filter(gene==gene_tar)
      splice_dat_small$junc_kmers <- vapply(splice_dat_small$peptide,calc_kmers,numeric(1))
      splice_dat_small$gene_kmers <- calc_kmers(splice_dat_small$peptide)
      kmer_counts_small <- kmer_counts_filt[vapply(splice_dat_small$rows,function(val){which(kmer_counts_filt$rows == val)},numeric(1)),]
      kmer_counts <- kmer_counts_small[,samples,drop=F]/splice_dat_small$junc_kmers
      gene_split <- strsplit(gene_tar,"-")[[1]]
      # gene_expr <- calc_gene_expression(gene_split,gene_expression_filt)
      junc_expr <- junc_expr_comb_filt[splice_dat_small$juncs,]
      dup_num <-nrow(splice_dat_small)
      junc_expr[junc_expr<0]<-NA
      a<-as.numeric(apply((kmer_counts*junc_expr),2,mean,na.rm=T))
    },numeric(length(samples)))))
    colnames(gene_metric_mean_tumor_no_gene_norm)<-samples    
  }
  
  saveRDS(gene_metric_mean_tumor_no_gene_norm,file=sprintf("%s_gene_metric_mean_len_norm_no_gene_norm_tumor.rds",out))
  write.table(gene_metric_mean_tumor_no_gene_norm,
              file=sprintf("%s_gene_metric_mean_len_norm_no_gene_norm_tumor.txt",out),quote=F, col.names = T, row.names = T, sep = "\t")

  coding_potential_LGC_tumor <- as.data.frame(vapply(genes_tumor,function(gene_tar){
    g<<-gene_tar
    splice_dat_small <- splice_dat_filt_tumor %>% dplyr::filter(gene==gene_tar)
    return(mean(splice_dat_small$coding_potential_LGC,na.rm=T))
  },numeric(1)))
  colnames(coding_potential_LGC_tumor)<-"coding_potential_LGC"

  saveRDS(coding_potential_LGC_tumor,file=sprintf("%s_coding_potential_LGC_tumor.rds",out))
  write.table(coding_potential_LGC_tumor,
              file=sprintf("%s_coding_potential_LGC_tumor.txt",out),quote=F,col.names=T,row.names=T,sep="\t")

  coding_potential_tumor <- as.data.frame(vapply(genes_tumor,function(gene_tar){
    g<<-gene_tar
    splice_dat_small <- splice_dat_filt_normal %>% dplyr::filter(gene==gene_tar)
    return(mean(splice_dat_small$coding_potential,na.rm=T))
  },numeric(1)))
  colnames(coding_potential_tumor)<-"coding_potential"

  saveRDS(coding_potential_tumor,file=sprintf("%s_coding_potential_tumor.rds",out))
  write.table(coding_potential_tumor,
              file=sprintf("%s_coding_potential_tumor.txt",out),quote=F,col.names=T,row.names=T,sep="\t")
  }
