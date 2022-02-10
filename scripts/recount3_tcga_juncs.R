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
                               help="The output directory for the tcga directory"))))

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
  file_source == "tcga"
)
proj_info <-  readRDS("/home/tpalme15/splicemutr_project/scripts/tcga/proj_info.rds")

#------------------------------------------------------------------------------#
# reading data in from each tcga project


for (iter in seq(2,nrow(proj_info))){
  # iter <- 31
  groups_file_T <- c()
  groups_file_N <- c()
  junc_files <- c()
  junc_rse <- create_rse(proj_info[iter,],type="jxn")
  junc_dat <- junc_rse@assays@data@listData[["counts"]]
  junc_metadata <- as.data.frame(junc_rse@colData@listData)
  rownames(junc_metadata)<-junc_metadata$external_id
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
  colnames(junc_dat) <- tcga_barcode
  cancer <- junc_metadata$study[1]

  all_juncs <- rownames(junc_dat)
  pre_junc_file <- create_junc_file(all_juncs)

  cancer_dir <- sprintf("%s/%s",tcga_junc_dir, cancer)

  if (!dir.exists(cancer_dir)){
    dir.create(cancer_dir)
  }
  for (i in seq(1,length(tcga_barcode))){
    file <- tcga_barcode[i]
    type <- sample_type[i]
    print(sprintf("%d:%s",iter,file))
    junc_file <- pre_junc_file
    junc_file$counts <- junc_dat[,file]
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
}
