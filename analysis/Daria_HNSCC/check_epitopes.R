#!/usr/bin/env Rscript

# purpose of this script is to check which amino acid epitopes used during validation can be found in the amino acid fasta files

#------------------------------------------------------------------------------#
# loading libraries

library(polyester)
library(Biostrings)
library(stringr)
library(seqinr)

#------------------------------------------------------------------------------#
# loading in the amino acid fasta files

target_aa<-"/media/theron/My_Passport/head_and_neck_DARIA/data/IEDB_benchmarking/HLA-A02-01/aa_target.fa"
filler_aa<-"/media/theron/My_Passport/reference_genomes/SEQUENCES/GENCODE/gencode.v36.pc_translations.fa"

target_aa<-readAAStringSet(target_aa)
filler_aa<-readAAStringSet(filler_aa)

#------------------------------------------------------------------------------#
# loading in the epitopes file

epitope_dir<-"/media/theron/My_Passport/head_and_neck_DARIA/data/IEDB_benchmarking/HLA-A02-01/epitopes.txt"
epitopes<-read.table(epitope_dir,header = T)
epitopes<-unique(epitopes$epitope) # turning dataframe into a character vector

#------------------------------------------------------------------------------#
# creating new fasta file

fasta_dir<-"/media/theron/My_Passport/head_and_neck_DARIA/data/IEDB_benchmarking/HLA-A02-01/validation_aa.fa"
writeXStringSet(target_aa,fasta_dir,append=FALSE)
writeXStringSet(filler_aa[1:200],fasta_dir,append=TRUE)

validation_fasta<-readAAStringSet(fasta_dir)

#------------------------------------------------------------------------------#
# creating new fasta file

names<-validation_fasta@ranges@NAMES
sequences<-unname(as.character(validation_fasta)) # sequences: vector

names_updated<-vapply(seq(length(names)),function(x){
  epis<-vapply(epitopes,function(y){
    if(str_detect(sequences[x],y)){
      y
    } else{"-"}
  },character(1))
  name<-sprintf("%s;",names[x])
  paste(c(name,unname(epis)),sep="",collapse="")
},character(1))

sequences<-lapply(sequences,function(x){x}) # sequences: list

#------------------------------------------------------------------------------#
# writing the amino acid validation fasta file

out_dir<-"/media/theron/My_Passport/head_and_neck_DARIA/data/IEDB_benchmarking/HLA-A02-01"
out<-sprintf("%s/validation_aa.fa",out_dir)
write.fasta(sequences,names_updated,file.out=out,open="w",nbchar=80,as.string=TRUE)
