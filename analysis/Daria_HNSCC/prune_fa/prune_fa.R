#!/usr/bin/env Rscript

# Purpose is to extract transcript and amino acid sequences from reference fasta files for a given gene

#------------------------------------------------------------------------------#
# loading libraries

library(Biostrings)
library(stringr)
library(seqinr)
library(optparse)

#------------------------------------------------------------------------------#
# loading in the fasta files

arguments <- parse_args(OptionParser(usage = "%prog [options] counts_file groups_file", 
                   description="form transcripts per junction for the given input junction file",
                   option_list=list(
                     make_option(c("-o","--output_directory"), default=NULL, help="The output directory for the data"),
                     make_option(c("-t","--transcript"), default=NULL, help="The .fa file for the protein-coding transcripts"),
                     make_option(c("-p","--protein"), default=NULL, help="The .fa file for the protein-coding translated transcripts"),
                     make_option(c("-g","--genes"), default=NULL, help=" The gene names in .txt file, one per line"))))

opt=arguments

out_dir<-opt$output_directory
transcript_dir<-opt$transcript
protein_dir<-opt$protein
gene_dir<-opt$genes

#------------------------------------------------------------------------------#
# loading in the transcript and protein fa, extracting names

transcripts<-readDNAStringSet(transcript_dir)
amino_acids<-readAAStringSet(protein_dir)

transcripts_names<-transcripts@ranges@NAMES
amino_acids_names<-amino_acids@ranges@NAMES

#------------------------------------------------------------------------------#
# loading in the reference files

genes<-read.table(gene_dir,header=TRUE)
genes<-genes$gene

#------------------------------------------------------------------------------#
# collecting genes per .fa line

transcripts_genes<-unname(vapply(transcripts_names,function(x){
  str_split(c(x),"[|]")[[1]][6]
},character(1)))

amino_acids_genes<-unname(vapply(amino_acids_names,function(x){
  str_split(c(x),"[|]")[[1]][7]
},character(1)))

#------------------------------------------------------------------------------#
# extracting sequences and names based on genes

transcript_loc<-which(transcripts_genes %in% genes)
amino_acid_loc<-which(amino_acids_genes %in% genes)

sequences_tx<-lapply(transcript_loc,function(x){
  unname(as.character(transcripts[x]))
})
sequences_aa<-lapply(amino_acid_loc,function(x){
  unname(as.character(amino_acids[x]))
})

names_tx<-vapply(amino_acid_loc,function(x){
  unname(as.character(transcripts_names[x]))
},character(1))
names_aa<-vapply(transcript_loc,function(x){
  unname(as.character(amino_acids_names[x]))
},character(1))

#------------------------------------------------------------------------------#
# writing .fa files

out_tx<-sprintf("%s/tx_target.fa",out_dir)
out_aa<-sprintf("%s/aa_target.fa",out_dir)
write.fasta(sequences_tx,names_tx,file.out=out_tx,open="w",nbchar=80,as.string=TRUE)
write.fasta(sequences_aa,names_aa,file.out=out_aa,open="w",nbchar=80,as.string=TRUE)
