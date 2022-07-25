# Theron Palmer
# 07/15/2021

# Looking at the mutation counts per sample per tcga cancer type

#------------------------------------------------------------------------------#
# loading in the appropriate libraries

library(ggplot2)
library(stringr)

#------------------------------------------------------------------------------#
# reading in mutation counts file

mutation_count_file <- "/media/theron/My_Passport/TCGA_junctions/cbioportal_data/Mutation_Count.txt"
mutation_counts <- read.table(mutation_count_file, sep="\t",skip=1)

#------------------------------------------------------------------------------#
# processing each sample for it's cancer type

mutation_counts$cancer <- vapply(mutation_counts$V1,function(sample){
  strsplit(sample,"[_]")[[1]][1]
},character(1))
colnames(mutation_counts) <- c("Study_ID","Patient_ID","Sample_ID","Mutation_Count","Cancer_Type")
  
#------------------------------------------------------------------------------#
# plotting mutation counts per cancer type

pdf(file="/media/theron/My_Passport/TCGA_junctions/cbioportal_data/mutations_per_cancer.pdf",width=10, height=10)

print(ggplot(mutation_counts,aes(x=Cancer_Type,y=log2(Mutation_Count)))+
        geom_boxplot()+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))

dev.off()