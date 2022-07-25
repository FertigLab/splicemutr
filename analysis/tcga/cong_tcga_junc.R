#!/usr/bin/env Rscript

# created: 05/15/2021

# The purpose of this script is to modify reference transcripts using
# reference transcript information

#------------------------------------------------------------------------------#
# loading libraries

library(stringr)

#------------------------------------------------------------------------------#
# reading in the necessary files

cancer_files <- read.table("/media/theron/My_Passport/TCGA_junctions/files.txt")
cancer_dir <- "/media/theron/My_Passport/TCGA_junctions"

#------------------------------------------------------------------------------#
# conglomerating the data.Rdata objects per cancer


for (i in seq(nrow(cancer_files))){
  print(i)
  data_file <- sprintf("%s/%s/data.Rdata",cancer_dir,cancer_files[i,])
  if (file.exists(data_file)){
    load(data_file)
  }
  filler <- rep(NA,nrow(introns))
  chr <- introns$chr
  start <- introns$start
  end <- introns$end
  clusters <- introns$clusterID
  strand <- unname(vapply(clusters,function(cluster){str_split(cluster,"[_]")[[1]][3]},character(1)))
  clusters_ref <- sprintf("clu_#_%s",strand)
  verdicts <- unname(vapply(introns$verdict,function(verdict){
    if (verdict == "unknown_strand"){
      return("unknown_strand")
    } else {
      return("N")
    }
  },character(1)))
  if (i == 1){
    all_juncs <- data.frame(clusters_ref,filler,filler,chr,start,end,verdicts,filler,filler)
  } else {
    temp_juncs <- data.frame(clusters_ref,filler,filler,chr,start,end,verdicts,filler,filler)
    all_juncs <- rbind(all_juncs,temp_juncs)
    all_juncs <- unique(all_juncs)
  }
}
colnames(all_juncs) <- colnames(introns)
write.table(all_juncs,file=sprintf("%s/all_introns.txt",cancer_dir),
            quote = F,col.names = TRUE, row.names = F)

#------------------------------------------------------------------------------#
# splitting all_juncs into 1000 junction bites

junc_size <- 1000
tot_juncs <- nrow(all_juncs)
iter <- 1
for (i in seq(1,tot_juncs,by=junc_size)){
  if ((i + junc_size) <= tot_juncs){
    small_juncs <- all_juncs[seq(i,(i+junc_size)),]
  } else {
    small_juncs <- all_juncs[seq(i,tot_juncs),]
  }
  split_introns_dir <- sprintf("%s/split_introns",cancer_dir)
  if (!dir.exists(split_introns_dir)){
    dir.create(split_introns_dir)
  }
  saveRDS(small_juncs,file=sprintf("%s/intron%d.rds",split_introns_dir,iter))
  iter <- iter + 1
}

