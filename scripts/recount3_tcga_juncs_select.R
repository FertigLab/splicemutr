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
                   make_option(c("-t","--study"),
                               default = sprintf("%s",getwd()),
                               help="The TCGA study to analyze"),
                   make_option(c("-o","--output_directory"),
                               default = sprintf("%s",getwd()),
                               help="The output directory for the tcga directory"),
                   make_option(c("-s","--samples"),
                               default = sprintf("%s",getwd()),
                               help="The samples to create files for"))))

opt=arguments

study <- opt$study
tcga_junc_dir <- opt$output_directory
samples <- unlist(read.table(opt$samples,sep="\t"))

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
# reading in the target TCGA cohort data

junc_rse <- recount3::create_rse_manual(
  project = study,
  project_home = "data_sources/tcga",
  organism = "human",
  annotation = "gencode_v26",
  type = "jxn"
)

#------------------------------------------------------------------------------#
# processing the study info

junc_files <- c()
junc_dat <- junc_rse@assays@data@listData[["counts"]]
junc_metadata <- as.data.frame(junc_rse@colData@listData)
rownames(junc_metadata)<-junc_metadata$external_id
tcga_barcode <- junc_metadata[colnames(junc_dat),4]
tcga_barcode_split <- as.data.frame(matrix(unlist(str_split(tcga_barcode,"[-]")),nrow=length(tcga_barcode),byrow=T))
colnames(junc_dat) <- tcga_barcode
cancer <- junc_metadata$study[1]

all_juncs <- rownames(junc_dat)
pre_junc_file <- create_junc_file(all_juncs)

cancer_dir <- sprintf("%s/%s",tcga_junc_dir, cancer)

if (!dir.exists(cancer_dir)){
  dir.create(cancer_dir)
}
for (i in seq(1,length(tcga_barcode))){
  if (tcga_barcode[i] %in% samples){
    file <- tcga_barcode[i]
    print(sprintf("%d:%s",iter,file))
    junc_file <- pre_junc_file
    junc_file$counts <- junc_dat[,file]
    tcga_file <- sprintf("%s/%s.junc",cancer_dir,file)
    write.table(junc_file,tcga_file,quote = F,sep="\t",col.names = F,row.names=F)
  }
}
