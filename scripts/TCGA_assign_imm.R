#title: "Daria Class 1 Analysis"
#author: "Theron Palmer"
#date: "08/06/2021"

library(ggplot2)
library(stringr)
library(pheatmap)
library(dplyr)
library(ggpubr)
library(msigdbr)
library(optparse)

#------------------------------------------------------------------------------#
# handling command line input

arguments <- parse_args(OptionParser(usage = "",
                 description="",
                 option_list=list(
                   make_option(c("-d","--dat_file"),
                               default = sprintf("%s",getwd()),
                               help="dat_file"),
                   make_option(c("-s","--splice_dat"),
                               default = sprintf("%s",getwd()),
                               help="splicemutr data"),
                   make_option(c("-m","--mut_count"),
                               default = sprintf("%s",getwd()),
                               help="mutation counts"),
                   make_option(c("-g","--genotypes_files"),
                               default = sprintf("%s",getwd()),
                               help="genotypes files"),
                   make_option(c("-a","--hla_files"),
                               default = sprintf("%s",getwd()),
                               help="hla files"))))
opt=arguments

dat_file <- opt$dat_file
dat_file <- sprintf("%s/data.Rdata",dat_file)
if (!file.exists(dat_file)){quit(save="no")}
load(dat_file)
splicemutr_file <- opt$splice_dat
mutation_count_file <- opt$mut_count
genotypes_files <- opt$genotypes_files
HLA_files <- opt$hla_files

#------------------------------------------------------------------------------#
# Internal functions

create_tcga_splicemutr <- function(introns,splicemutr_dat){
  juncs_to_find <- sprintf("%s:%s:%s",introns$chr,introns$start,introns$end)
  rownames(introns) <- juncs_to_find
  specific_splicemutr_dat <- subset(splicemutr_dat,juncs %in% juncs_to_find)
  specific_splicemutr_dat$cluster <- introns[specific_splicemutr_dat$juncs,"clusterID"]
  specific_splicemutr_dat$verdict <- introns[specific_splicemutr_dat$juncs,"verdict"]
  specific_splicemutr_dat$deltapsi <- introns[specific_splicemutr_dat$juncs,"deltapsi"]
  return(specific_splicemutr_dat)
}

#------------------------------------------------------------------------------#
# reading in the data

splicemutr_dat <- read.table(splicemutr_file, sep = " ",header=T)

mutation_counts <- read.table(mutation_count_file, sep="\t",header=T)

#------------------------------------------------------------------------------#
#Assigning filenames to genotypes

genotypes <- read.table(genotypes_file,sep=",",header=T)
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

mutation_IDs <- vapply(genotypes$aliquot_id,function(ID){
  id_parts <- str_split(ID,"[-]")[[1]]
  sprintf("%s-%s-%s",id_parts[1],id_parts[2],id_parts[3])
},character(1))
genotypes$mutation_IDs <- mutation_IDs

mut_counts <- vapply(genotypes$mutation_IDs,function(ID){
  counts <<- mutation_counts[mutation_counts$Patient.ID==ID,"Mutation.Count"]
  if (length(counts)==0){
    counts <- -1
  } else {
    counts <- max(counts)
  }
  return(counts)
},numeric(1))

cancer_type <- unname(vapply(genotypes$mutation_IDs,function(ID){
  cancer <- str_split(mutation_counts[mutation_counts$Patient.ID==ID,"Study.ID"],"[_]")
  if (length(cancer)==0){
    cancer <- "NONE"
  } else {
    cancer <- cancer[[1]][1]
  }
},character(1)))

tum_or_norm <- unname(vapply(genotypes$aliquot_id,function(ID){
  type <- str_split(ID,"[-]")[[1]][4]
  type<-as.numeric(substr(type,1,2))
  if (type >= 1 & type <= 9){
    return("T")
  } else if (type >= 10 & type <= 19) {
    return("N")
  }
  return("C")
},character(1)))
genotypes$mut_counts <- mut_counts
genotypes$cancer_type <- cancer_type
genotypes$tum_or_norm <- tum_or_norm

genotypes <- genotypes %>% dplyr::filter(cancer_type != "NONE")
cancer_types <- unique(genotypes$cancer_type)

genotypes <- genotypes %>% dplyr::filter(cancer_type == tolower(basename(dirname(dat_file))))

#------------------------------------------------------------------------------#
# calculating average tumor and average normal scores per sample

dat <- data.frame(matrix(0,nrow=nrow(splicemutr_dat),ncol=nrow(genotypes)))
vec <- vector("list",nrow(splicemutr_dat))

for (i in seq(1,nrow(genotypes))){
  row_vec <- rep(F,nrow(splicemutr_dat))
  print(i)
  class1_alleles <- as.character(unname(genotypes[i,seq(6)]))
  iter <- 0
  for (allele in class1_alleles){
    print(allele)
    file_HLA <- sprintf(HLA_files,str_replace(allele,":","-"))
    if (!file.exists(file_HLA)){next}
    info = file.info(file_HLA)
    if (info$size == 0){next}
    iter <- iter + 1
    HLA_dat <- read.table(file_HLA)
    HLA_dat[,1]<-as.numeric(HLA_dat[,1])+1
    counts <- lapply(str_split(HLA_dat[,3],":"),length)
    kmers_split <- str_split(HLA_dat[,2],":")
    rows <- HLA_dat[,1]
    row_vec[rows]<-T
  }
  counts <- counts/iter
  if (length(counts)==0){
    dat[row_vec,i] <- 0
  } else {
    dat[row_vec,i] <- dat[row_vec,i] + counts
  }
}
colnames(dat) <- genotypes$aliquot_id
splicemutr_dat <- cbind(splicemutr_dat,dat)

# saveRDS(dat,file=sprintf("%s/%s",dirname(dat_file),"specific_counts.rds"))

#------------------------------------------------------------------------------#
# creating the specific splicemutr data

specific_splicemutr_dat <- create_tcga_splicemutr(introns,splicemutr_dat)
write.table(specific_splicemutr_dat,
            file=sprintf("%s/%s_splicemutr.txt",dirname(dat_file),basename(dirname(dat_file))),
            sep="\t",
            quote=F,
            col.names=T,
            row.names=F)
