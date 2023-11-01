#!/usr/bin/env Rscript

#------------------------------------------------------------------------------#
# loading libraries

library(dplyr)
library(stringr)

#------------------------------------------------------------------------------#
# handling command line input

arguments <- parse_args(OptionParser(usage ="",
  description="create TCGA genotypes using optitype file",
  option_list=list(
   make_option(c("-c","--cancer"), default = sprintf("%s",getwd()), help="The output directory for the kmer data"),
   make_option(c("-g","--genotypes"), default=NULL, help="The txdb object"),
   make_option(c("-o","--output_dir"), default=NULL, help="The txdb object"))))

opt=arguments

output_dir<-opt$output_dir
cancer<-opt$cancer
genotypes_file<-opt$genotypes

#------------------------------------------------------------------------------#
# scripting practice files

output_dir <- "/Users/tpalme15/Desktop/splicemutr_paper/TCGA/intermediates"
cancer <- "CHOL"
genotypes_file <- "/Users/tpalme15/Desktop/splicemutr_paper/SKCM_PRIM_MET/input_data/OptiTypeCallsHLA_20171207.tsv"

#------------------------------------------------------------------------------#
# reading in files

genotypes <- read.csv(genotypes_file,header=T)
rownames(genotypes)<-genotypes$aliquot_id

junc_rse <- recount3::create_rse_manual(
  project = cancer,
  project_home = "data_sources/tcga",
  organism = "human",
  annotation = "gencode_v26",
  type = "jxn"
)

#------------------------------------------------------------------------------#
# extracting cancer metadata

junc_metadata <- as.data.frame(junc_rse@colData@listData)
tcga_barcode <- junc_metadata$tcga.tcga_barcode

#------------------------------------------------------------------------------#
# filtering genotypes

genotypes_specific <- genotypes[tcga_barcode,]
genotypes_specific$aliquot_id<-tcga_barcode

#------------------------------------------------------------------------------#
# creating genotypes_specific

HLA_columns <- seq(6)
aliquot_id_column <- ncol(genotypes)
HLAS <- genotypes_specific[,HLA_columns]
HLAs_reformatted <- vapply(seq(nrow(HLAS)),function(row_val){
  HLA_values <- HLAS[row_val,]
  paste(vapply(HLA_values,function(HLA){str_replace(str_replace(sprintf("HLA-%s",HLA),"[*]",""),":","-")},character(1)),collapse=",")
},character(1))

genotypes_reformatted <- data.frame(HLAS=HLAs_reformatted,sample=genotypes_specific$aliquot_id)
genotypes_reformatted[genotypes_reformatted == "HLA-NA,HLA-NA,HLA-NA,HLA-NA,HLA-NA,HLA-NA"]<-paste(rep(NA,6),collapse=",")

#------------------------------------------------------------------------------#
# saving genotypes_specific

write.table(genotypes_reformatted,
            file=sprintf("%s/%s_genotypes_specfic.txt",output_dir,cancer),
            sep="\t",
            col.names = F,
            row.names = F,
            quote=F)



