#!/usr/bin/env Rscript

# converting splicemutr summary ensembl gene names to gene symbols using gtf file

# created: 05/16/2022

#------------------------------------------------------------------------------#
# loading libraries


#------------------------------------------------------------------------------#
# reading in summary file

splicemutr_summary_file <- "F:/head_and_neck_DARIA/data/splicemutr_05_26_2021/GENE_METRIC_01032022/tumor_gene_metric_all.csv"
splicemutr_summary <- read.table(splicemutr_summary_file,sep=",")
ense_gene_names <- rownames(splicemutr_summary)

splicemutr_dat_file <- "F:/head_and_neck_DARIA/data/splicemutr_05_26_2021/GENE_METRIC_01032022/full_splicemutr_dat.rds"
splicemutr_dat <- readRDS(splicemutr_dat_file)
splicemutr_dat_genes <- splicemutr_dat$gene

#------------------------------------------------------------------------------#
# reading in associated gtf file

gtf_file <- "F:/reference_genomes/SEQUENCES/GENCODE/GRCh38_Ensembl99_sparseD3_sjdbOverhang99/Homo_sapiens.GRCh38.99.gtf"
gtf <- read.table(gtf_file,sep="\t")

gtf_genes <- gtf[gtf$V3=="gene",]
rm(gtf)

gtf_map <- data.frame(t(vapply(gtf_genes$V9,function(val){
  a<-strsplit(val,";")
  vapply(a[[1]],function(a_val){
    a_vals<-unlist(strsplit(a_val," "))
    return(a_vals[length(a_vals)])
  },character(1))
},character(5))))

gtf_map <- gtf_map[,c(1,3)]
colnames(gtf_map)<-c("ensembl","gene_name")
rownames(gtf_map)<-gtf_map$ensembl
gtf_map <- unique(gtf_map)

splicemutr_summary_gene_names <- vapply(ense_gene_names,function(ens){
  a<-strsplit(ens,"-")[[1]]
  if (length(a)>1){
    return(paste(c(gtf_map[a[1],2]),gtf_map[a[2],2],sep="-"))
  } else {
    return(gtf_map[a[1],2])
  }
},character(1))

splicemutr_summary$gene_name<-splicemutr_summary_gene_names

write.table(splicemutr_summary,
            file="F:/head_and_neck_DARIA/data/splicemutr_05_26_2021/GENE_METRIC_01032022/tumor_gene_metric_all_gene_names.csv",
            quote=F,
            sep=",",
            col.names = T,
            row.names = F)

splicemutr_dat_gene_names <- vapply(splicemutr_dat_genes,function(ens){
  a<-strsplit(ens,"-")[[1]]
  if (length(a)>1){
    return(paste(c(gtf_map[a[1],2]),gtf_map[a[2],2],sep="-"))
  } else {
    return(gtf_map[a[1],2])
  }
},character(1))

splicemutr_dat$gene_name <- splicemutr_dat_gene_names
write.table(splicemutr_dat,
            file="F:/head_and_neck_DARIA/data/splicemutr_05_26_2021/GENE_METRIC_01032022/full_splicemutr_dat_gene_names.txt",
            quote=F,
            sep="\t",
            col.names = T,
            row.names = F)

