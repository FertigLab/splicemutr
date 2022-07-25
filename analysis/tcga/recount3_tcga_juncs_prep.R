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
                               help="iter"))))
opt=arguments

tcga_junc_dir <- opt$output_directory
iter <- as.numeric(opt$iter)

#------------------------------------------------------------------------------#
# internal functions

create_junc_file <- function(all_juncs){
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
  file_source %in% c("tcga","gtex")
)
if (iter == 1){
  saveRDS(proj_info,"/home/tpalme15/splicemutr_project/scripts/tcga/proj_info.rds")
}

#------------------------------------------------------------------------------#
# reading data in from each tcga project

groups_file_T <- c()
groups_file_N <- c()
junc_files <- c()
cancer <- proj_info[iter,"project"]
cancer_dir <- sprintf("%s/%s",tcga_junc_dir, cancer)
if (!dir.exists(cancer_dir)){
  dir.create(cancer_dir)
}
junc_rse <- create_rse(proj_info[iter,],type="jxn")
junc_metadata <- as.data.frame(junc_rse@colData@listData)

saveRDS(junc_rse,file=sprintf("%s/%s_junc_rse.rds",cancer_dir,basename(cancer_dir)))
saveRDS(junc_metadata,file=sprintf("%s/%s_metadata.rds",cancer_dir,basename(cancer_dir)))