#!/usr/bin/env Rscript

# The purpose of this script is to simulate RNAseq for sequences used for validation of splicemutr. 

#------------------------------------------------------------------------------#
# loading libraries and setting seed

dir<-"/media/theron/My_Passport/splicemutr"
source(sprintf("%s/%s",dir,"libraries.R"))
source(sprintf("%s/%s",dir,"functions.R"))
seed=27

#------------------------------------------------------------------------------#
# loading in the fasta files

target_tx<-"/media/theron/My_Passport/head_and_neck_DARIA/data/IEDB_benchmarking/HLA-A02-01/tx_target.fa"
filler_tx<-"/media/theron/My_Passport/reference_genomes/SEQUENCES/GENCODE/gencode.v36.pc_transcripts.fa"

target_tx<-readDNAStringSet(target_tx)
filler_tx<-readDNAStringSet(filler_tx)

#------------------------------------------------------------------------------#
# obtaining gene chromosomes from filler fasta file

names<-filler_tx@ranges@NAMES
genes<-unname(vapply(names,function(x){
  gene<-str_split(x,"[|]")[[1]][2]
  str_split(gene,"[.]")[[1]][1]
},character(1)))

mouse = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
chromosomes<-getBM(mart = mouse, attributes = c("chromosome_name","ensembl_gene_id"), filter = "ensembl_gene_id", values = genes, uniqueRows=FALSE)
chr<-unique(chromosomes$chromosome_name)
chr<-chr[seq(length(chr)-1)]

#------------------------------------------------------------------------------#
# randomly selecting 30 genes from each chromosome

gene_locations<-unlist(lapply(chr,function(x){
  gene_loc<-which(chromosomes$chromosome_name == x)
  set.seed(27)
  sample(gene_loc,30)
}))

#------------------------------------------------------------------------------#
# creating new fasta file

fasta_dir<-"/media/theron/My_Passport/head_and_neck_DARIA/data/IEDB_benchmarking/HLA-A02-01/validation.fa"
writeXStringSet(target_tx,fasta_dir,append=FALSE)
writeXStringSet(filler_tx[gene_locations],fasta_dir,append=TRUE)

validation_fasta<-readDNAStringSet(fasta_dir)

#------------------------------------------------------------------------------#
# setting up the uniform counts matrix

num_samples = 12 # 6 per normal and not normal
countmat = matrix(200, nrow=length(validation_fasta), ncol=num_samples)

#------------------------------------------------------------------------------#
# modifying the uniform counts matrix for differential splicing

tumor<-seq(7,12)
normal<-seq(1,6)

target<-seq(length(target_tx))

filler<-sample(seq(length(target_tx)+1,length(validation_fasta)),100)

for (x in seq(length(validation_fasta))){
  test <- (x %% 2) == 0
  if (test){
    countmat[x,tumor]<-round(rnorm(length(tumor),mean=2^(log2(200)+1),sd=5))*(countmat[x,tumor]/200)
  } else {
    countmat[x,normal]<-round(rnorm(length(normal),mean=2^(log2(200)+1),sd=5))*(countmat[x,normal]/200)
  }
}

countmat_ggplot<-data.frame(c(countmat[1:20,1],countmat[1:20,12]))

countmat_ggplot$isoform<-NA
countmat_ggplot$isoform[1:20]<-1:20
countmat_ggplot$isoform[21:40]<-1:20
countmat_ggplot$group<-NA
countmat_ggplot$group[1:20]<-"normal"
countmat_ggplot$group[21:40]<-"tumor"

colnames(countmat_ggplot)<-c("reads","isoform","group")

ggplot(data=countmat_ggplot, aes(x=isoform, y=reads, group=group)) +
  geom_line(aes(color=group))+
  geom_point(aes(color=group))

#------------------------------------------------------------------------------#
# simulate reads:

simulate_experiment_countmat(fasta_dir, readmat=countmat, paired = T, seed = 27,
                             outdir='/media/theron/My_Passport/head_and_neck_DARIA/data/IEDB_benchmarking/HLA-A02-01/simulated_RNAseq_data')
saveRDS(countmat,file="/media/theron/My_Passport/head_and_neck_DARIA/data/IEDB_benchmarking/HLA-A02-01/simulated_RNAseq_data/countmat.rds")

