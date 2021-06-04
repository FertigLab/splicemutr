#!/usr/bin/env Rscript

# created: 06/03/2021

# The purpose of this script is to extract the peptides from the splicemutr output files

#------------------------------------------------------------------------------#
# loading libraries

library(optparse)
library(dplyr)

#------------------------------------------------------------------------------#
# handling command line input

arguments <- parse_args(OptionParser(usage = "",
               description="",
               option_list=list(
                 make_option(c("-o","--output_directory"), default = sprintf("%s",getwd()), help="The output directory for the peptide data"),
                 make_option(c("-d","--data"), default=NULL, help="The splicemutr data"))))

opt=arguments

out_dir<-opt$output_directory
data_file<-opt$data

#------------------------------------------------------------------------------#
# loading in the splicemutr data file

data_files <- read.table(data_file)

for (i in seq(nrow(data_files))){
  if (i == 1){
    data <- read.table(data_files[i,])
  } else {
    data <- rbind(data,read.table(data_files[i,]))
  }
}

data_protein_coding <- data %>% dplyr::filter(protein_coding == "Yes")
proteins <- data_protein_coding$peptide
protein_df <- data.frame(proteins)

#------------------------------------------------------------------------------#
# saving data

write.table(proteins, file=sprintf("%s/proteins.txt",out_dir),
            col.names = F, row.names = F, quote = F)
write.table(data_protein_coding, file=sprintf("%s/data_splicemutr_protein_coding.txt",out_dir),
            col.names = F,row.names = F,quote = F,sep="\t")

