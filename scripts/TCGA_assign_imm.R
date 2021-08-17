#!/usr/bin/env Rscript

#author: "Theron Palmer"
#date: "08/06/2021"

library(stringr)
library(dplyr)
library(optparse)
library(TCGAutils)

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
                               help="hla files"),
                   make_option(c("-r","--groups_file"),
                               default = sprintf("%s",getwd()),
                               help="groups file"),
                   make_option(c("-j","--junc_rse"),
                               default = sprintf("%s",getwd()),
                               help="junc_rse"))))
opt=arguments

dat_file <- opt$dat_file
if (!file.exists(dat_file)){
  print(dat_file)
  quit(save="no")
}
load(dat_file)
splicemutr_file <- opt$splice_dat
mutation_count_file <- opt$mut_count
genotypes_files <- opt$genotypes_files
HLA_files <- opt$hla_files
groups_file <- opt$groups_file
junc_rse_file <- opt$junc_rse

#------------------------------------------------------------------------------#
# Internal functions

create_tcga_splicemutr <- function(introns,splicemutr_dat){
  juncs_to_find <- sprintf("%s:%s:%s",introns$chr,introns$start,introns$end)
  rownames(introns) <- juncs_to_find
  splicemutr_dat$juncs <- sprintf("%s:%s:%s",splicemutr_dat$chr,splicemutr_dat$start,splicemutr_dat$end)
  specific_splicemutr_dat <- subset(splicemutr_dat,juncs %in% juncs_to_find)
  specific_splicemutr_dat$cluster <- introns[specific_splicemutr_dat$juncs,"clusterID"]
  specific_splicemutr_dat$verdict <- introns[specific_splicemutr_dat$juncs,"verdict"]
  specific_splicemutr_dat$deltapsi <- introns[specific_splicemutr_dat$juncs,"deltapsi"]
  return(specific_splicemutr_dat)
}

#------------------------------------------------------------------------------#
# reading in the data

# genotypes_files <- "/media/theron/My_Passport/TCGA_junctions/ext_dat/OptiTypeCallsHLA_20171207.tsv"
# HLA_files <- "/media/theron/My_Passport/TCGA_junctions/summary_files/%s_tx_dict_summary_perc.txt"
# splicemutr_file <- "/media/theron/My_Passport/TCGA_junctions/formed_transcripts/data_splicemutr.txt"
# mutation_count_file <-"/media/theron/My_Passport/TCGA_junctions/cbioportal_data/Mutation_Count.txt"
# dat_file <- "/media/theron/My_Passport/TCGA_junctions/TCGA_cancers/CHOL/data.Rdata"
# load(dat_file)
# groups_file <- "/media/theron/My_Passport/TCGA_junctions/TCGA_cancers/CHOL/groups_file.txt"
# junc_rse_file <- "/media/theron/My_Passport/TCGA_junctions/TCGA_cancers/CHOL/juncrse.rds"

leafcutter_groups <- read.table(groups_file)
colnames(leafcutter_groups) <- c("external_id","type")
splicemutr_dat <- read.table(splicemutr_file, sep = " ",header=T)
mutation_counts <- read.table(mutation_count_file, sep="\t",header=T)
genotypes <- read.table(genotypes_files,sep=",",header=T)
junc_rse <- readRDS(junc_rse_file)
junc_metadata <- as.data.frame(junc_rse@colData@listData)

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

mutation_IDs <- vapply(genotypes$aliquot_id,function(ID){
  id_parts <- str_split(ID,"[-]")[[1]]
  sprintf("%s-%s-%s",id_parts[1],id_parts[2],id_parts[3])
},character(1))
genotypes$mutation_IDs <- mutation_IDs

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
genotypes$cancer_type <- cancer_type
genotypes$tum_or_norm <- tum_or_norm

genotypes <- genotypes %>% dplyr::filter(cancer_type != "NONE")
cancer_types <- unique(genotypes$cancer_type)

# genotypes <- genotypes %>% dplyr::filter(cancer_type == tolower(basename(dirname(dat_file))))

#------------------------------------------------------------------------------#
# determining the cancer type per mutation count sample

mutation_counts$cancer <- vapply(mutation_counts$Study.ID,
                                 function(ID){str_split(ID,"[_]")[[1]][1]},
                                 character(1))

# mutation_counts <- mutation_counts %>% dplyr::filter(cancer == tolower(basename(dirname(dat_file))))

#------------------------------------------------------------------------------#
# creating uniform names

metadata_rows <- unlist(lapply(leafcutter_groups$external_id,
                               function(external_id){
                                 which(junc_metadata$external_id == external_id)
                               }))

leafcutter_groups$sample_id <- vapply(leafcutter_groups$external_id,
                                  function(ID){
                                    barcode <- TCGAbarcode(junc_metadata$tcga.tcga_barcode[which(junc_metadata$external_id == ID)],
                                                sample=T)
                                    substr(barcode,1,nchar(barcode)-1)
                                  },character(1))
genotypes$sample_id <- TCGAbarcode(genotypes$aliquot_id, sample=T)
genotypes$sample_id <- vapply(genotypes$sample_id,function(ID){substr(ID,1,nchar(ID)-1)},character(1))

#------------------------------------------------------------------------------#
# joining leafcutter, genotypes, and mutation_count data to the leafcutter groups information

cols <- c("A1","A2","B1","B2","C1","C2")
HLA_alleles <- lapply(seq(length(leafcutter_groups$sample_id)),function(num){
  ID <- leafcutter_groups$sample_id[num]
  type <- leafcutter_groups$type[num]
  external_id <- leafcutter_groups$external_id[num]
  a<-genotypes[genotypes$sample_id == ID,cols]
  if (nrow(a) == 0){
    a <- data.frame(t(rep(NA,6)))
    colnames(a) <-colnames(genotypes[,seq(6)])
  }
  a$sample_id <- ID
  a$type <- type
  a$external_id <- external_id
  return(a)
})
genotypes_leafcutter <- data.frame(HLA_alleles[[1]])
for (i in seq(2,length(HLA_alleles))){
    genotypes_leafcutter <- rbind(genotypes_leafcutter,
                                  HLA_alleles[[i]])
}
genotypes_leafcutter <- unique(genotypes_leafcutter)

genotypes_leafcutter$mut_counts<- vapply(genotypes_leafcutter$sample_id,function(ID){
  count <- mutation_counts[mutation_counts$Sample.ID == ID,"Mutation.Count"]
  if (length(count) == 0){
    return(0)
  } else {
    return(count)
  }
},numeric(1))


#------------------------------------------------------------------------------#
# calculating average tumor and average normal scores per sample

# dat <- data.frame(matrix(0,nrow=nrow(splicemutr_dat),ncol=nrow(genotypes_leafcutter)))
# vec <- vector("list",nrow(splicemutr_dat))
# gene_row <- nrow(genotypes_leafcutter)
# for (i in seq(1,gene_row)){
#   row_vec <- rep(F,gene_row)
#   print(sprintf("%s:%d:%d",basename(dirname(dat_file)),i,gene_row))
#   class1_alleles <- as.character(unname(genotypes_leafcutter[i,seq(6)]))
#   iter <- 0
#   kmers <- vector(mode = "list", length = nrow(splicemutr_dat))
#   for (allele in class1_alleles){
#     file_HLA <- sprintf(HLA_files,str_replace(allele,":","-"))
#     if (!file.exists(file_HLA)){next}
#     info = file.info(file_HLA)
#     if (info$size == 0){next}
#     iter <- iter + 1
#     HLA_dat <- read.table(file_HLA)
#     HLA_dat[,1]<-as.numeric(HLA_dat[,1])+1
#     rows <- HLA_dat[,1]
#     kmers_split<-str_split(HLA_dat[,2],":")
#     fill <- vapply(seq(length(rows)),function(row){
#       kmers[[rows[row]]]<<-unique(c(kmers[[rows[row]]],kmers_split[[row]]))
#       return(T)
#     },logical(1))
#   }
#   kmers <- unlist(lapply(seq(length(kmers)),function(val){
#     paste(kmers[[val]],collapse=":")
#   }))
#   dat[,i] <- kmers
# }
# cols_to_remove <- which(is.na(genotypes_leafcutter$A1))
# tumor_cols <- genotypes_leafcutter$type == "T"
# tumor_cols[cols_to_remove] <- F
# normal_cols <- genotypes_leafcutter$type == "N"
# normal_cols[cols_to_remove] <- F
#
#
# splicemutr_dat <- cbind(splicemutr_dat, dat)

#------------------------------------------------------------------------------#
# creating the specific splicemutr data

# specific_splicemutr_dat <- create_tcga_splicemutr(introns,splicemutr_dat)
#
# write.table(specific_splicemutr_dat,
#             file=sprintf("%s/%s_splicemutr.txt",dirname(dat_file),basename(dirname(dat_file))),
#             sep="\t",
#             quote=F,
#             col.names=T,
#             row.names=F)

write.table(genotypes_leafcutter,
            file=sprintf("%s/%s_genotypes.txt",dirname(dat_file),basename(dirname(dat_file))),
            sep="\t",
            quote=F,
            col.names=T,
            row.names=F)

write.table(junc_metadata,
            file=sprintf("%s/%s_metadata.txt",dirname(dat_file),basename(dirname(dat_file))),
            sep="\t",
            quote=F,
            col.names=T,
            row.names=F)

