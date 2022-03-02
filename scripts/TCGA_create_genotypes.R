#!/usr/bin/env Rscript

#author: "Theron Palmer"
#date: "08/06/2021"

library(stringr)
library(dplyr)
library(optparse)
# library(TCGAutils)

#------------------------------------------------------------------------------#
# handling command line input

arguments <- parse_args(OptionParser(usage = "",
                 description="",
                 option_list=list(
                   make_option(c("-m","--junc_metadata"),
                               default = sprintf("%s",getwd()),
                               help="junc metadata"),
                   make_option(c("-g","--genotypes_files"),
                               default = sprintf("%s",getwd()),
                               help="genotypes files"))))
opt=arguments

genotypes_files <- opt$genotypes_files
junc_metadata_file<- opt$junc_metadata

#------------------------------------------------------------------------------#
# reading in the data

# genotypes_files <- "/media/theron/My_Passport/TCGA_junctions/ext_dat/OptiTypeCallsHLA_20171207.tsv"
# junc_metadata_file <- "/media/theron/My_Passport/TCGA_junctions/TCGA_cancers/BRCA/BRCA_metadata.rds"

genotypes <- read.table(genotypes_files,sep=",",header=T)
junc_metadata <- readRDS(junc_metadata_file)

#------------------------------------------------------------------------------#
# assigning filenames to genotypes

genotypes$A1 <- unname(vapply(genotypes$A1,function(HLA){
  HLA_split<-str_split(HLA,"[*:]")[[1]]
  sprintf("HLA-%s%s-%s",HLA_split[1],HLA_split[2],HLA_split[3])
},character(1)))
genotypes$A2 <- unname(vapply(genotypes$A2,function(HLA){
  HLA_split<-str_split(HLA,"[*:]")[[1]]
  sprintf("HLA-%s%s-%s",HLA_split[1],HLA_split[2],HLA_split[3])
},character(1)))
genotypes$B1 <- unname(vapply(genotypes$B1,function(HLA){
  HLA_split<-str_split(HLA,"[*:]")[[1]]
  sprintf("HLA-%s%s-%s",HLA_split[1],HLA_split[2],HLA_split[3])
},character(1)))
genotypes$B2 <- unname(vapply(genotypes$B2,function(HLA){
  HLA_split<-str_split(HLA,"[*:]")[[1]]
  sprintf("HLA-%s%s-%s",HLA_split[1],HLA_split[2],HLA_split[3])
},character(1)))
genotypes$C1 <- unname(vapply(genotypes$C1,function(HLA){
  HLA_split<-str_split(HLA,"[*:]")[[1]]
  sprintf("HLA-%s%s-%s",HLA_split[1],HLA_split[2],HLA_split[3])
},character(1)))
genotypes$C2 <- unname(vapply(genotypes$C2,function(HLA){
  HLA_split<-str_split(HLA,"[*:]")[[1]]
  sprintf("HLA-%s%s-%s",HLA_split[1],HLA_split[2],HLA_split[3])
},character(1)))

#------------------------------------------------------------------------------#
# matching TCGA samples between mutation counts and genotypes

genotypes_specific <- as.data.frame(matrix(unlist(lapply(junc_metadata$tcga.tcga_barcode,function(barcode){
  a<-which(barcode == genotypes$aliquot_id)
  if (length(a)==0){
    return(rep(NA,6))
  } else {
    return(as.character(genotypes[a[1],seq(6)]))
  }
})),byrow=T,nrow=nrow(junc_metadata)))
colnames(genotypes_specific)<-colnames(genotypes)[seq(6)]

tum_or_norm <- unname(vapply(junc_metadata$tcga.tcga_barcode,function(ID){
  type <- str_split(ID,"[-]")[[1]][4]
  type<-as.numeric(substr(type,1,2))
  if (type >= 1 & type <= 9){
    return("T")
  } else if (type >= 10 & type <= 19) {
    return("N")
  }
  return("C")
},character(1)))
genotypes_specific$aliquot_id <- junc_metadata$tcga.tcga_barcode
genotypes_specific$type <- tum_or_norm

# genotypes_specific$sample_id <- TCGAbarcode(genotypes$aliquot_id, sample=T)
# genotypes_specific$sample_id <- vapply(genotypes$sample_id,function(ID){substr(ID,1,nchar(ID)-1)},character(1))

genotypes_specific$external_id <- junc_metadata$external_id
alleles <- data.frame(alleles=unique(c(genotypes_specific$A1,genotypes_specific$A2,genotypes_specific$B1,genotypes_specific$B2,genotypes_specific$C1,genotypes_specific$C2)))
alleles <- alleles[!is.na(alleles$alleles),]

#------------------------------------------------------------------------------#
# creating genotypes list

geno_names <- genotypes_specific$external_id
genotype_data <- lapply(geno_names,function(ID){
  return(as.character(genotypes_specific[genotypes_specific$external_id==ID,c("A1","A2","B1","B2","C1","C2")]))
})
names(genotype_data) <- geno_names


#------------------------------------------------------------------------------#
# saving the genotypes data
write.table(genotypes_specific,
            file=sprintf("%s/%s_genotypes_specific.txt",dirname(junc_metadata_file),basename(dirname(junc_metadata_file))),
            sep="\t",
            quote=F,
            col.names=T,
            row.names=F)
write.table(alleles,
            file=sprintf("%s/%s_class1_alleles.txt",dirname(junc_metadata_file),basename(dirname(junc_metadata_file))),
            sep="\t",
            quote=F,
            col.names=F,
            row.names=F)
saveRDS(genotype_data,file=sprintf("%s/%s_genotypes.rds",dirname(junc_metadata_file),basename(dirname(junc_metadata_file))))
