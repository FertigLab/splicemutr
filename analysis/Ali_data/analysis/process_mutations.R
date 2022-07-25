#!/usr/bin/env Rscript

#------------------------------------------------------------------------------#
# loading libraries

library(Biostrings)
library(stringr)

#------------------------------------------------------------------------------#
# files to read in

maf_file <- "F:/Ali_data/WES Data for 4T1MIS/MB01JHU503/MB01JHU503_000_analysis/MAF/DNA_4TIMSI_haplotypecaller.maf"
MAF <- read.csv(maf_file,sep="\t",header=T)
protein_fasta_file <- "F:/reference_genomes/SEQUENCES/ENSEMBL/Mus_musculus.GRCm38.pep.all.fa"
protein_fasta <- readAAStringSet(protein_fasta_file)

#------------------------------------------------------------------------------#
# filtering for MAF protein-coding

MAF_filtered <- MAF[MAF$BIOTYPE=="protein_coding" & MAF$Variant_Classification =="Missense_Mutation",]
refs <- unlist(lapply(MAF_filtered$Amino_acids,function(val){
  a<-str_split(val,"/")[[1]]
  if (length(a)==2){
    return(a[1])
  } else {
    return(c("NA"))
  }
}))
mods <- unlist(lapply(MAF_filtered$Amino_acids,function(val){
  a<-str_split(val,"/")[[1]]
  if (length(a)==2){
    return(a[2])
  } else {
    return(c("NA"))
  }
}))

MAF_filtered$ref<-refs
MAF_filtered$mod<-mods

change_pos <- unlist(lapply(MAF_filtered$Protein_position,function(val){
  a<-str_split(val,"/")[[1]]
  if (length(a)==2){
    return(a[1])
  } else {
    return(c("NA"))
  }
}))
length <- unlist(lapply(MAF_filtered$Protein_position,function(val){
  a<-str_split(val,"/")[[1]]
  if (length(a)==2){
    return(a[2])
  } else {
    return(c("NA"))
  }
}))

MAF_filtered$pep_change_pos<-as.numeric(change_pos)
MAF_filtered$pep_length<-as.numeric(length)

#------------------------------------------------------------------------------#
# processing names for gene and transcript name

fasta_names <- names(protein_fasta)
genes_transcripts <- data.frame(t(vapply(fasta_names,function(val){
  a<-str_split(val," ")[[1]][c(4,5)]
  return(c(str_replace(a[1],"gene:",""),str_replace(a[2],"transcript:","")))
},character(2))))
rownames(genes_transcripts) <- seq(nrow(genes_transcripts))
colnames(genes_transcripts) <- c("gene","transcript")
genes_transcripts$gene <- vapply(genes_transcripts$gene,function(val){return(substr(val,1,nchar(val)-2))},character(1))
genes_transcripts$transcript <- vapply(genes_transcripts$transcript,function(val){return(substr(val,1,nchar(val)-2))},character(1))

names(protein_fasta)<- genes_transcripts$transcript
proteins_transcripts <- as.character(protein_fasta)

#------------------------------------------------------------------------------#
# checking transcript length# 2766/15300 mutations have proteins that can't be found in the reference, none have unequal lengths

transcripts <- MAF_filtered$Transcript_ID
lengths<-MAF_filtered$pep_length

testing_lengths<-vapply(seq(length(transcripts)),function(val){
  transcript <- transcripts[val]
  len<-lengths[val]
  peptide <- proteins_transcripts[transcript]
  if (is.na(peptide)){
    return(0)
  } else if(len==nchar(peptide)){
    return(1)
  } else {
    return(2)
  }
},numeric(1))

#------------------------------------------------------------------------------#
# modifying proteins

MAF_further_filtered <- MAF_filtered[testing_lengths==1,]
peptides<-proteins_transcripts[MAF_further_filtered$Transcript_ID]
K<-9
modded_sequences <- data.frame(t(vapply(seq(length(peptides)),function(val){
  peptide<-peptides[val]
  ref<-MAF_further_filtered$ref[val]
  mod<-MAF_further_filtered$mod[val]
  pos<-MAF_further_filtered$pep_change_pos[val]
  len<-MAF_further_filtered$pep_length[val]
  upstream_pos<-1
  downstream_pos<-len
  if (is.na(pos)){
    return(c("NA","NA"))
  }
  if ((pos-(K-1))>upstream_pos){
    upstream_pos<-pos-(K-1)
  }
  if ((pos+(K-1))<downstream_pos){
    downstream_pos<-pos+(K-1)
  }
  ref_check<-substr(peptide,pos,pos)
  if(ref != ref_check){
    print(sprintf("%s:%s",ref,ref_check))
  }
  return(c(paste(c(substr(peptide,upstream_pos,pos-1),mod,substr(peptide,pos+1,downstream_pos)),collapse="-"),
           paste(c(substr(peptide,upstream_pos,pos-1),mod,substr(peptide,pos+1,downstream_pos)),collapse="")))
},character(2))))
colnames(modded_sequences)<-c("highlighted","embedded")
MAF_with_sequences <- cbind(MAF_further_filtered,modded_sequences)

#------------------------------------------------------------------------------#
# saving peptides as fasta

sequences <- MAF_with_sequences$embedded
names(sequences)<-sprintf("mafrow:%d;%s;%s",
                             seq(nrow(MAF_with_sequences)),
                                 MAF_with_sequences$Hugo_Symbol,
                                 MAF_with_sequences$Transcript_ID)
for (i in seq(1,length(sequences),5000)){
  if (i+4999>length(sequences)){
    sequences_small<-AAStringSet(sequences[seq(i,length(sequences))])
  } else {
    sequences_small<-AAStringSet(sequences[seq(i,i+4999)])
  }
  writeXStringSet(sequences_small,
                  sprintf("F:/Ali_data/WES Data for 4T1MIS/MB01JHU503/MB01JHU503_000_analysis/MAF/mutation_peptides_%d.fa",i))
  
}
sequences_all<-AAStringSet(sequences)
writeXStringSet(sequences_all,"F:/Ali_data/WES Data for 4T1MIS/MB01JHU503/MB01JHU503_000_analysis/MAF/mutation_peptides.fa")
write.table(MAF_with_sequences,
            file="F:/Ali_data/WES Data for 4T1MIS/MB01JHU503/MB01JHU503_000_analysis/MAF/MAF_wsequences.maf",
            col.names=T,row.names=F,quote=F,sep="\t")

