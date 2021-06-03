#!/usr/bin/env Rscript

# created: 06/03/2021

# The purpose of this script is to parse and conglomerate the arcashla json files

#------------------------------------------------------------------------------#
# loading libraries

library(rjson)
library(optparse)
library(stringr)
library(stringi)

#------------------------------------------------------------------------------#
# command line input

arguments <- parse_args(OptionParser(usage = "%prog [options] counts_file groups_file",
               description="form transcripts per junction for the given input junction file",
               option_list=list(
                 make_option(c("-o","--output_directory"), default = sprintf("%s",getwd()), help="The output directory for the kmer data"),
                 make_option(c("-g","--genotypes"), default=NULL, help="The file that contains the location of all genotype files"))))

opt=arguments

out_dir<-opt$output_directory
genotype_file<-opt$genotypes

#------------------------------------------------------------------------------#
# internal functions

#' modifies the arcashla genotype list to a vector of alleles of mhcnuggets format
#' @param genotype_list the list of genotypes for the target HLA alleles
#' @return returns a vector of potential class I and class II alleles and pairings
#'
mhcnuggets_format <- function(genotype_list){
  allele_types <- names(genotype_list)
  mhcnuggets_allele_list <- lapply(allele_types,function(allele){
    alleles<-genotype_list[[allele]]
    alleles_split<-str_split(str_replace(alleles,"[*]",""),"[:]")
    lapply(alleles_split,function(splits){
      alleles_mhcnuggets <- sprintf("HLA-%s:%s",splits[1],splits[2])
    })
  })
  DP_pairs <- c()
  DQ_pairs <- c()
  DR_pairs <- c()
  names(mhcnuggets_allele_list)<-allele_types
  if (all(c("DPA1","DPB1") %in% allele_types)){
    pairs <- data.frame(expand.grid(unlist(mhcnuggets_allele_list[["DPA1"]]),unlist(mhcnuggets_allele_list[["DPB1"]])))
    DP_pairs<-sprintf("%s-%s",pairs[,1],pairs[,2])
  }
  if (all(c("DQA1","DQB1") %in% allele_types)){
    pairs <- data.frame(expand.grid(unlist(mhcnuggets_allele_list[["DQA1"]]),unlist(mhcnuggets_allele_list[["DQB1"]])))
    DQ_pairs<-sprintf("%s-%s",pairs[,1],pairs[,2])
  }
  if (all(c("DRA","DRB1") %in% allele_types)){
    pairs <- data.frame(expand.grid(unlist(mhcnuggets_allele_list[["DRA"]]),unlist(mhcnuggets_allele_list[["DRB1"]])))
    DR_pairs<-sprintf("%s-%s",pairs[,1],pairs[,2])
  }
  return(c(unname(unlist(mhcnuggets_allele_list)),DP_pairs,DQ_pairs,DR_pairs))
}


#------------------------------------------------------------------------------#
# reading in the genotypes files

genotypes <- read.table(genotype_file)

genotype_vecs<-lapply(seq(nrow(genotypes)),function(i){
  json_file <- genotypes[i,]
  json_base <- str_split(basename(json_file),"[.]")[[1]][1]
  json_base_file <- sprintf("%s_dir/%s.genotype.json",json_file,json_base)
  genotype <- fromJSON(file=json_base_file)
  genotype_vec <- mhcnuggets_format(genotype)
})

saveRDS(genotype_vecs,file=sprintf("%s/genotype_vecs.rds",out_dir))
genotypes <- data.frame(unique(unlist(genotype_vecs)))
colnames(genotypes)<-"alleles"
class2<- str_detect(genotypes$alleles,"[D]")
class2_genotypes<-genotypes[class2,]
class1_genotypes<-genotypes[!class2,]

write.table(class1_genotypes,file=sprintf("%s/allele_data_class1.txt",out_dir),col.names = F, row.names = F, quote=F)
write.table(class2_genotypes,file=sprintf("%s/allele_data_class2.txt",out_dir),col.names = F, row.names = F, quote=F)

