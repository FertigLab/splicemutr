#!/usr/bin/env Rscript

library(dplyr)
library(stringr)
library(optparse)

#------------------------------------------------------------------------------#
# handling command line input

arguments <- parse_args(OptionParser(usage = "",
                   description="",
                   option_list=list(
                     make_option(c("-g","--genotype_file"),
                                 default = "",
                                 help="the genotype file"),
                     make_option(c("-s","--sample_file"),
                                 default = "",
                                 help="the sample file"),
                     make_option(c("-a","--allele_dir"),
                                 default = "",
                                 help="the allele directory"),
                     make_option(c("-f","--full_splice_file"),
                                 default = "",
                                 help="the full splice dat"),
                     make_option(c("-d","--small_splice_file"),
                                 default = "",
                                 help="the protein-coding splice data"))))
opt=arguments
genotype_file <- opt$genotype_file
sample_file <- opt$sample_file
allele_dir <- opt$allele_dir
full_splice_file <- opt$full_splice
small_splice_file <- opt$small_splice

#------------------------------------------------------------------------------#
# internal functions

create_filler <- function(splice_dat,samples){
  filler_num <- nrow(splice_dat)
  filler_mat <- matrix(rep(0,filler_num*length(samples)),nrow=filler_num,ncol=length(samples))
  colnames(filler_mat)<-samples
  return(as.data.frame(filler_mat))
}

create_full_kmers <- function(splice_dat,samples){
  total_rows <- nrow(splice_dat)
  num_samples <- length(samples)+2
  strand <- as.data.frame(matrix(unlist(str_split(splice_dat$cluster,"_")),byrow=T,nrow=nrow(splice_dat)))[,3]
  juncs <- sprintf("%s:%s:%s:%s",splice_dat$chr,splice_dat$start,splice_dat$end,strand)
  filler_mat <- matrix(rep(0,total_rows*(length(samples)+2)),nrow=total_rows,ncol=length(samples)+2)
  colnames(filler_mat)<-c("rows","junc",samples)
  filler_mat[,"rows"] <- seq(total_rows)
  filler_mat[,"junc"] <- juncs
  return(as.data.frame(filler_mat))
}

format_alleles <- function(alleles){
  alleles <- str_replace(alleles,":","-")
  alleles <- sprintf("%s_tx_dict_summary_perc.txt",alleles)
}

# #------------------------------------------------------------------------------#
# # local files
#
# genotype_file <- "/media/theron/My_Passport/head_and_neck_DARIA/data/splicemutr_05_26_2021/processed_data/genotype_vecs.rds"
# sample_file <- "/media/theron/My_Passport/head_and_neck_DARIA/data/splicemutr_05_26_2021/samples.txt"
# allele_dir <- "/media/theron/My_Passport/head_and_neck_DARIA/data/splicemutr_05_26_2021/processed_data/class1"
# full_splice_file <- "/media/theron/My_Passport/head_and_neck_DARIA/data/splicemutr_05_26_2021/data_splicemutr.txt"
# small_splice_file <- "/media/theron/My_Passport/head_and_neck_DARIA/data/splicemutr_05_26_2021/data_splicemutr_protein_coding.txt"

#------------------------------------------------------------------------------#
# reading in the genotype and samples file

genotypes <- readRDS(genotype_file)
samples <- read.table(sample_file,header=F)
names(genotypes) <- samples$V1
rm(samples)
full_splice <- read.table(full_splice_file, sep=" ",header=T)
full_splice$deltapsi <- as.numeric(full_splice$deltapsi)
small_splice <- read.table(small_splice_file, sep="\t",header=F)
colnames(small_splice) <- colnames(full_splice)

#------------------------------------------------------------------------------#
# reading in the genotypes per file

filler_mat <- create_filler(small_splice,names(genotypes))
a<-vapply(names(genotypes),function(sample){
  alleles <- genotypes[[sample]]
  alleles <- alleles[which(str_detect(alleles,"HLA-A") | str_detect(alleles,"HLA-B") | str_detect(alleles,"HLA-C"))]
  alleles <- format_alleles(alleles)
  vapply(alleles[1],function(allele){
    allele_dat <- read.table(sprintf("%s/%s",allele_dir,allele),header=F,sep="\t")
    filler_mat[as.numeric(allele_dat$V1)+1,sample] <<- as.numeric(filler_mat[as.numeric(allele_dat$V1)+1,sample]) +
      vapply(str_split(allele_dat$V2,":"),function(dat){length(dat)},numeric(1))
    return(T)
  },logical(1))
  return(T)
},logical(1))

#------------------------------------------------------------------------------#
# creating full kmer file

full_splice$rows <- seq(nrow(full_splice))
small_splice_proto <- full_splice %>% dplyr::filter(protein_coding=="Yes")
full_kmers <- create_full_kmers(full_splice,names(genotypes))
full_kmers[small_splice_proto$rows,colnames(genotypes)] <- filler_mat
full_splice <- full_splice %>% dplyr::filter(!(verdict == "annotated" & modified == "changed") & deltapsi > 0)
full_kmers <- full_kmers[full_splice$rows,]

rm(small_splice_proto,small_splice,genotypes)

saveRDS(full_kmers,file=sprintf("%s/full_kmers.rds",dirname(full_splice_file)))
write.table(full_kmers,
            file=sprintf("%s/full_kmers.txt",dirname(full_splice_file)),quote=F, col.names = T, row.names = F, sep = "\t")

saveRDS(full_splice,file=sprintf("%s/full_splicemutr_dat.rds",dirname(full_splice_file)))
write.table(full_splice,
            file=sprintf("%s/full_splicemutr_dat.txt",dirname(full_splice_file)),quote=F, col.names = T, row.names = F, sep = "\t")
