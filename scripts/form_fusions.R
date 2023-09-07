#!/usr/bin/env Rscript

# created: 05/15/2021
# updated: 09/29/2021

# The purpose of this script is to modify reference transcripts using
# reference transcript information

# build in exons next to one another function checking if the target exons are directly next to one another in the reference transcript
# this implies that the intron modification is only taking place in a single reference transcript
# look to see if there are different

#------------------------------------------------------------------------------#
# loading libraries

library(GenomicRanges)
library(rtracklayer)
library(annotate)
library(org.Hs.eg.db)
library(ensembldb)
library(biomaRt)
library(DataCombine)
library(stringi)
library(stringr)
library(optparse)
library(dplyr)
library(rlist)
library(AnnotationDbi)
library(Biostrings)
# library(splicemute)

#------------------------------------------------------------------------------#
# handling command line input

arguments <- parse_args(OptionParser(usage = "%prog [options] counts_file groups_file",
                                     description="form transcripts per junction for the given input junction file",
                                     option_list=list(
                                       make_option(c("-o","--out_dir"), default = sprintf("%s",getwd()), help="The output directory for the kmer data"),
                                       make_option(c("-t","--txdb"), default=NULL, help="The txdb object"),
                                       make_option(c("-j","--juncs"), default=NULL, help="The junction file (path and file)"),
                                       make_option(c("-b","--bsgenome_name"),default = "bsgenome_name",help="the bsgenome object name"),
                                       make_option(c("-m","--chr_map"),default=F,help="the chromosome map"),
                                       make_option(c("-f","--funcs"),default=F,help="functions"))))

opt=arguments

out_dir<-opt$out_dir
txdb_file<-opt$txdb
chr_map<-opt$chr_map

bsgenome_name <- opt$bsgenome_name # bsgenome_name so that can create bsgenome object
library(bsgenome_name,character.only=T) # assigning bsgenome object to "bsgenome" variable
assign("bsgenome",get(bsgenome_name))

print("intron_file stage")
intron_file <- opt$juncs
if (str_detect(intron_file,".rds")){
  introns <-readRDS(opt$juncs) # loading in the introns data
} else {
  introns <- read.table(intron_file,header=T)
}

print("chr_map stage")
introns$chr <- str_replace(introns$chr,"chr","") # replacing the chr with ""
if (typeof(chr_map)!="logical"){ # this is done in the case that there is a specific chromosome mapping that needs to be done
  chr_map_file <- opt$chr_map
  chr_map <- read.table(chr_map_file)
  colnames(chr_map)<-c("chr","mapped")
  introns$chr <- vapply(introns$chr,function(chrom){
    return(chr_map$mapped[chr_map$chr==chrom])
  },character(1))
}
# introns<-format_introns(introns)
introns_HPVfusions <-  introns[(introns$chr_donorA == "K02718.1" & introns$chr_acceptorB != "K02718.1") | (introns$chr_donorA != "K02718.1" & introns$chr_acceptorB == "K02718.1"),]

source(opt$funcs)

#------------------------------------------------------------------------------#
# local input

library(BSgenome)
library(BSgenome.Hsapienshpv16.GENCODE.GRCh38.107)
bsgenome <- BSgenome.Hsapienshpv16.GENCODE.GRCh38.107
# intron_file <-"F:/head_and_neck_DARIA/CalifanoHPVOP/star_fusion_out/DGay14-39174_trimedChimeric.out.junction"
intron_file <-"F:/head_and_neck_DARIA/CalifanoHPVOP/star_fusion_out/DGay21-090-1_1_valChimeric.out.junction"
# intron_file <- "F:/head_and_neck_DARIA/CalifanoHPVOP/star_fusion_out/DGay18-047-1_1_valChimeric.out.junction"
introns <- read.table(intron_file,header=T)
# introns_HPVfusions <-  introns[(introns$chr_donorA == "K02718.1" & introns$chr_acceptorB != "K02718.1") | (introns$chr_donorA != "K02718.1" & introns$chr_acceptorB == "K02718.1"),]

# txdb_file <- "F:/reference_genomes/human_hpv16/Homo_sapiens.GRCh38.107.hpv16_mod.txdb"
# source("F:/splicemute/R/functions.R")
# chr_map <- F
# out_dir <- "F:/head_and_neck_DARIA/CalifanoHPVOP/star_fusion_out"

#------------------------------------------------------------------------------#
# preparing the references for transcript formation and kmerization

print("reading in txdb")
txdb<-loadDb(txdb_file) # making the txdb from gtf
if (typeof(chr_map)!="logical"){ # read in the gtf information but map the chromosomes
  all_genes<-map_chroms(genes(txdb),chr_map)
  exons_by_gene<-map_chroms(exonsBy(txdb,by="gene"),chr_map)
  exons_by_tx<-map_chroms(exonsBy(txdb,by=c("tx"),use.names=T),chr_map)
  tx_by_gene<-map_chroms(transcriptsBy(txdb,by="gene"),chr_map)
  five_by_tx<-map_chroms(fiveUTRsByTranscript(txdb,use.names=T),chr_map)
  three_by_tx<-map_chroms(threeUTRsByTranscript(txdb,use.names=T),chr_map)
  cds_by_tx <- map_chroms(cdsBy(txdb,by="tx",use.names=T),chr_map)
}
all_genes<-genes(txdb)
exons_by_gene<-exonsBy(txdb,by="gene")
exons_by_tx<-exonsBy(txdb,by=c("tx"),use.names=T)
tx_by_gene<-transcriptsBy(txdb,by="gene")
five_by_tx<-fiveUTRsByTranscript(txdb,use.names=T)
three_by_tx<-threeUTRsByTranscript(txdb,use.names=T)
cds_by_tx <- cdsBy(txdb,by="tx",use.names=T)

#------------------------------------------------------------------------------#
# preparing the references for transcript formation and kmerization

data_canon<-data.frame()
sequences<-c()
cds_storage <- list()
intron_length<-nrow(introns_HPVfusions)
seq_vals<-data.frame(junction=NA,virus_tx=NA,sequence=NA)
for (i in seq(1,intron_length)){
  print(sprintf("%d introns out of %d total introns",i,intron_length))
  curr_introns<-introns_HPVfusions[i,] # get the current intron
  donor_info <- c(curr_introns$chr_donorA,curr_introns$brkpt_donorA,curr_introns$brkpt_donorA,curr_introns$strand_donorA)
  acceptor_info <- c(curr_introns$chr_acceptorB,curr_introns$brkpt_acceptorB,curr_introns$brkpt_acceptorB,curr_introns$strand_acceptorB)
  if (donor_info[1]==acceptor_info[1]){
    next
  }
  if (donor_info[1]=="K02718.1"){
    exon_donor_identity <- choose_fusion_exons(donor_info,exons_by_tx,T,T)[1] # this is an error
  } else {
    exon_donor_identity <- choose_fusion_exons(donor_info,exons_by_tx,F,T)[1] # this is an error
  }
  if (acceptor_info[1]=="K02718.1"){
    exon_acceptor_identity <- choose_fusion_exons(acceptor_info,exons_by_tx,T,F)[1]
  } else {
    exon_acceptor_identity <- choose_fusion_exons(acceptor_info,exons_by_tx,F,F)[1]
  }
  if (is.na(exon_acceptor_identity) | is.na(exon_donor_identity)){
    next
  }
  donor_flanking_exons <- get_flanking_exons(donor_info,exon_donor_identity,exons_by_tx,T) # get flanking exons
  acceptor_flanking_exons <- get_flanking_exons(acceptor_info,exon_acceptor_identity,exons_by_tx,F) # get flanking exons
  donor_flanker_merged <- merge_granges(donor_flanking_exons,exon_donor_identity,donor_info,T)
  acceptor_flanker_merged <- merge_granges(acceptor_flanking_exons,exon_acceptor_identity,acceptor_info,F)
  seq_vals_temp <- merge_donor_acceptor(bsgenome,donor_flanker_merged,acceptor_flanker_merged,cds_by_tx)
  # print(seq_vals_temp)
  seq_val_name <- sprintf("%s;%s",
                          paste(c(donor_info,names(exon_donor_identity)),collapse=":"),
                          paste(c(acceptor_info,names(exon_acceptor_identity)),collapse=":"))
  if (length(which(!is.na(seq_vals_temp)))==0){
    a<-data.frame(junction=seq_val_name,virus_tx=NA,sequence=NA)
    seq_vals <- rbind(seq_vals,a)
  } else {
    seq_vals_temp_df <- data.frame(junction=rep(seq_val_name,length(seq_vals_temp)),
                                   virus_tx=names(seq_vals_temp),sequence=unname(seq_vals_temp))
    seq_vals <- rbind(seq_vals,seq_vals_temp_df)
  }
}

seq_vals <- seq_vals[complete.cases(seq_vals),]

#------------------------------------------------------------------------------#
# saving seq_vals data

write.table(seq_vals,
            file=sprintf("%s/%s.fusion.peptides",out_dir,str_remove(basename(intron_file),"Chimeric.out.junction")),
            quote=F,sep="\t",col.names = T,row.names = F)
sequences <- seq_vals$sequence
names(sequences) <- sprintf("seq_%d",seq(length(sequences)))
sequences<-AAStringSet(sequences)
writeXStringSet(sequences,sprintf("%s/%s.fa",out_dir,str_remove(basename(intron_file),"Chimeric.out.junction")))
sequences <- data.frame(peps=seq_vals$sequence)
write.table(sequences,
            file=sprintf("%s/%s.peps.fa",out_dir,str_remove(basename(intron_file),"Chimeric.out.junction")),
            quote=F,sep="\t",col.names = T,row.names = F)
