# creator: Theron Palmer
# created: 07/27/2021

# extracting and transforming the tcga alleles to be used in splicemutr

#------------------------------------------------------------------------------#
# loading necessary libraries

library(stringr)

#------------------------------------------------------------------------------#
# reading in allele file

# file <- "/media/theron/My_Passport/TCGA_junctions/OptiTypeCallsHLA_bams_20170323.tsv"
file <- "/media/theron/My_Passport/TCGA_junctions/OptiTypeCallsHLA_20171207.tsv"
allele_files <- read.table(file,sep=",",header=T)

#------------------------------------------------------------------------------#
# extracting alleles

HLA_A <- c(allele_files$A1,allele_files$A2)
HLA_B <- c(allele_files$B1,allele_files$B2)
HLA_C <- c(allele_files$C1,allele_files$C2)
HLA_all <- c(HLA_A,HLA_B,HLA_C)

#------------------------------------------------------------------------------#
# modifing allles for use with mhcnuggets

mhcnuggets_alleles <- data.frame(unique(unname(vapply(HLA_all,function(HLA){
  HLA_split<-str_split(HLA,"[*:]")[[1]]
  sprintf("HLA-%s%s:%s",HLA_split[1],HLA_split[2],HLA_split[3])
},character(1)))))
colnames(mhcnuggets_alleles)<-"allele"

out_dir <- "/media/theron/My_Passport/TCGA_junctions/alleles"
for (i in seq(nrow(mhcnuggets_alleles))){
  print(i)
  out_file <- sprintf("%s/tcga_alleles%d.txt",out_dir,i)
  write.table(mhcnuggets_alleles[i,],
              file=out_file,
              quote=F,
              col.names = F,
              row.names = F)
}

