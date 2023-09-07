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
# reading in the recount3 tcga objects

human_projects <- available_projects()
proj_info <- subset(
  human_projects,
  file_source == "tcga"
)

#------------------------------------------------------------------------------#
# reading data in from each tcga project

all_cancer_dirs <- c()
for (iter in seq(1,nrow(proj_info))){
  groups_file_T <- c()
  groups_file_N <- c()
  junc_files <- c()
  junc_rse <- create_rse(proj_info[iter,],type="jxn")
  junc_metadata <- as.data.frame(junc_rse@colData@listData)
  rownames(junc_metadata)<-junc_metadata$external_id
  cancer <- junc_metadata$study[1]

  cancer_dir <- sprintf("%s/%s",tcga_junc_dir, cancer)
  all_cancer_dirs<-c(all_cancer_dirs,cancer)

  if (!dir.exists(cancer_dir)){
    dir.create(cancer_dir)
  }

  saveRDS(junc_rse,file=sprintf("%s/%s/junc_rse.rds",tcga_junc_dir,cancer))
  saveRDS(junc_metadata,file=sprintf("%s/%s/junc_metadata.rds",tcga_junc_dir,cancer))
}

all_cancer_dirs <-data.frame(all_cancer_dirs)
colnames(all_cancer_dirs)<-"cancer_dir"
write.table(all_cancer_dirs,file=sprintf("%s/cancer_dirs.txt",tcga_junc_dir),col.names = F,row.names = F,quote = F)


