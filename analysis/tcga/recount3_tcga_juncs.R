#!/usr/bin/env Rscript

# created: 07/03/2021

# Taking the recount3 junctions and building .junc files per cancer type and tumor

#------------------------------------------------------------------------------#
# loading libraries

library(recount3)
library(stringr)
library(optparse)

#------------------------------------------------------------------------------#
# handling command line input

arguments <- parse_args(OptionParser(usage = "%prog [options] counts_file groups_file",
                 description="form transcripts per junction for the given input junction file",
                 option_list=list(
                   make_option(c("-o","--output_directory"), 
                               default = sprintf("%s",getwd()), 
                               help="The output directory for the tcga directory"),
                   make_option(c("-i","--iter"), 
                               default = sprintf("%s",getwd()), 
                               help="iter"),
                   make_option(c("-g","--gtex"), 
                               default = sprintf("%s",getwd()), 
                               help="gtex directory"),
                   make_option(c("-m","--tcga_gtex_map"),
                               default = sprintf("%s",getwd()),
                               help="tcga_gtex map"),
                   make_option(c("-r","--include"),
                               default = sprintf("%s",getwd()),
                               help="cancers to include"))))
opt=arguments

tcga_junc_dir <- opt$output_directory
iter <- as.numeric(opt$iter)
gtex_dir <- opt$gtex
tcga_gtex_map_file <- opt$tcga_gtex_map
cancers_to_include_file <- opt$include

#------------------------------------------------------------------------------#
# internal functions

create_junc_file <- function(all_juncs,junc_dat){
  split_data <- as.data.frame(matrix(unlist(str_split(all_juncs,"[:]")),
                                     nrow = nrow(junc_dat),
                                     byrow = T))
  coords <- as.data.frame(matrix(unlist(str_split(split_data[,2],"[-]")),nrow=length(all_juncs),byrow = T))
  colnames(coords)<-c("start","end")
  coords$start <- as.numeric(coords$start)
  coords$end <- as.numeric(coords$end)
  
  split_data <- cbind(split_data, coords)
  junc_file <- data.frame(split_data[,1],
                          split_data$start-1,
                          split_data$end+1,
                          sprintf("J%d",seq(nrow(split_data))),
                          rep(0,nrow(split_data)),split_data[,3])
  colnames(junc_file) <- c("chr","start","end","junc_name","counts","strand")
  return(junc_file)
}

#------------------------------------------------------------------------------#
# reading in the recount3 tcga objects

human_projects <- available_projects()
proj_info <- subset(
  human_projects,
  file_source %in% c("tcga")
)
# saveRDS(proj_info,sprintf("%s/proj_info.rds",tcga_junc_dir))

proj_info<-readRDS(sprintf("%s/proj_info.rds",tcga_junc_dir))

#------------------------------------------------------------------------------#
# reading data in from each tcga project

groups_file_T <- c()
groups_file_N <- c()
junc_files <- c()
cancer <- proj_info[iter,"project"]
cancers_to_include <- read.table(cancers_to_include_file)
cancers_to_include <- cancers_to_include$V1

if (cancer %in% cancers_to_include){
  
  cancer_dir <- sprintf("%s/%s",tcga_junc_dir, cancer)
  if (!dir.exists(cancer_dir)){
    dir.create(cancer_dir)
  }
  junc_rse <- readRDS(sprintf("%s/%s_junc_rse.rds",cancer_dir,cancer))
  junc_dat <- junc_rse@assays@data@listData[["counts"]]
  junc_metadata <- readRDS(sprintf("%s/%s_metadata.rds",cancer_dir,cancer))
  rownames(junc_metadata)<-junc_metadata$external_id
  print(sprintf("%s/%s_metadata.txt",cancer_dir,basename(cancer_dir)))
  
  tcga_barcode <- junc_metadata[colnames(junc_dat),4]
  tcga_barcode_split <- as.data.frame(matrix(unlist(str_split(tcga_barcode,"[-]")),nrow=length(tcga_barcode),byrow=T))
  sample_type <- vapply(tcga_barcode_split[,4],function(sample){
    sample_int <- as.numeric(substr(sample,1,2))
    if (sample_int < 10){
      return("T")
    } else if (sample_int < 20 & sample_int >= 10){
      return("N")
    } else if (sample_int >= 20){
      return("C")
    }
  },character(1))
  
  colnames(junc_dat) <- junc_metadata$external_id
  
  all_juncs <- rownames(junc_dat)
  pre_junc_file <- create_junc_file(all_juncs,junc_dat)
  
  if (!dir.exists(cancer_dir)){
    dir.create(cancer_dir)
  }
  
  for (i in seq(1,length(tcga_barcode))){
    file <- junc_metadata$external_id[i]
    type <- sample_type[i]
    print(sprintf("%d:%s",iter,file))
    junc_file <- pre_junc_file
    junc_file[,"counts"] <- junc_dat[,file]
    tcga_file <- sprintf("%s/%s.junc",cancer_dir,file)
    write.table(junc_file,tcga_file,quote = F,sep="\t",col.names = F,row.names=F)
    junc_files <- c(junc_files,tcga_file)
    if (type == "T"){
      groups_file_T <- rbind(groups_file_T,data.frame(t(c(file,"T"))))
    } else if (type == "N"){
      groups_file_N <- rbind(groups_file_N,data.frame(t(c(file,"N"))))
    }
  }
  
  groups_file <- rbind(groups_file_N, groups_file_T)
  groups_file <- groups_file[!is.na(groups_file[,1]),]
  write.table(junc_files, file=sprintf("%s/junc_file.txt",cancer_dir), quote=F, col.names = F, row.names = F, sep = "\t")
  write.table(groups_file, file=sprintf("%s/groups_file.txt",cancer_dir), quote=F, col.names = F, row.names = F, sep = "\t")
} else {
  print(sprintf("%s is not to be analyzed",cancer))
}
