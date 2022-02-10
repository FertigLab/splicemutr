#!/usr/bin/env Rscript

#------------------------------------------------------------------------------#
# loading libraries

library(optparse)

#------------------------------------------------------------------------------#
# handling command line input

arguments <- parse_args(OptionParser(usage = "%prog [options] counts_file groups_file",
               description="form transcripts per junction for the given input junction file",
               option_list=list(
                 make_option(c("-o","--out_dir"), default = sprintf("%s",getwd()), help="The output directory for the kmer data"),
                 make_option(c("-s","--splice_files"), default=NULL, help="The txdb object"))))

opt=arguments
out_dir <- opt$out_dir
splice_files <- opt$splice_files

#------------------------------------------------------------------------------#
# combining and extracting proteins

splicemutr_files<-read.table(splice_files,header=F)
splicemutr_files$nums <- vapply(splicemutr_files$V1,function(file){as.numeric(strsplit(basename(file),"_")[[1]][4])},numeric(1))
splicemutr_files <- splicemutr_files[order(splicemutr_files$nums),]
a<-vapply(seq(nrow(splicemutr_files)),function(row_val){
  print(row_val)
  if (row_val==1){
    splice_data <<- readRDS(splicemutr_files$V1[row_val])
  } else {
    splice_data <<- rbind(splice_data,readRDS(splicemutr_files$V1[row_val]))
  }
  return(T)
},logical(1))
splice_data <- splice_data %>% dplyr::filter(protein_coding == "Yes")
splice_data_ann <- splice_data %>% dplyr::filter(error == "tx" & annotated == "annotated")
splice_data <- splice_data %>% dplyr::filter(annotated != "annotated")
splice_data <- rbind(splice_data,splice_data_ann)

saveRDS(splice_data,file=sprintf("%s/%s",out_dir,"data_splicemutr_all.rds"))
write.table(splice_data,file=sprintf("%s/%s",out_dir,"data_splicemutr_all.txt"),
            col.names = T,quote=F,row.names=F)
splice_data <- splice_data %>% dplyr::filter(!is.na(peptide))

proteins <- data.frame(splice_data$peptide)

write.table(proteins,file=file=sprintf("%s/%s",out_dir,"proteins.txt"),col.names = F,sep="\t",quote=F,row.names = F)
saveRDS(splice_data,file=sprintf("%s/%s",out_dir,"data_splicemutr_all_pep.rds"))
write.table(splice_data,file=sprintf("%s/%s",out_dir,"data_splicemutr_all_pep.txt"),
            col.names=T,quote=F,row.names=F)
