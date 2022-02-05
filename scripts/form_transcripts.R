#!/usr/bin/env Rscript

# created: 05/15/2021
# updated: 09/29/2021

# The purpose of this script is to modify reference transcripts using
# reference transcript information

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
# library(splicemute)
library(AnnotationDbi)

#------------------------------------------------------------------------------#
# handling command line input

arguments <- parse_args(OptionParser(usage = "%prog [options] counts_file groups_file",
  description="form transcripts per junction for the given input junction file",
  option_list=list(
   make_option(c("-o","--out_prefix"), default = sprintf("%s",getwd()), help="The output directory for the kmer data"),
   make_option(c("-t","--txdb"), default=NULL, help="The txdb object"),
   make_option(c("-j","--juncs"), default=NULL, help="The junction file (path and file)"),
   make_option(c("-f","--funcs"), default=NULL, help="The splicemute functions to source"),
   make_option(c("-b","--bsgenome_name"),default = "bsgenome_name",help="the bsgenome object name"),
   make_option(c("-m","--chr_map"),default=F,help="the chromosome map"))))

opt=arguments

out_prefix<-opt$out_prefix
txdb_file<-opt$txdb
funcs<-opt$funcs
source(funcs)
# chr_map_file <- opt$chr_map
# chr_map <- read.table(chr_map_file)
# colnames(chr_map)<-c("chr","mapped")

bsgenome_name <- opt$bsgenome_name

introns <-readRDS(opt$juncs) # loading in the introns data
introns$chr <- str_replace(introns$chr,"chr","")
# if (typeof(chr_map)!="logical"){
#   introns$chr <- vapply(introns$chr,function(chrom){
#     return(chr_map$mapped[chr_map$chr==chrom])
#   },character(1))
# }
introns<-format_introns(introns)

#------------------------------------------------------------------------------#
# assigning bsgenome object to "bsgenome" variable

library(bsgenome_name,character.only = T)
assign("bsgenome",get(bsgenome_name))

#------------------------------------------------------------------------------#
# preparing the references for transcript formation and kmerization

print("reading in txdb")
txdb<-loadDb(txdb_file) # making the txdb from gtf
# if (typeof(chr_map)!="logical"){
#   all_genes<-map_chroms(genes(txdb),chr_map)
#   exons_by_gene<-map_chroms(exonsBy(txdb,by="gene"),chr_map)
#   exons_by_tx<-map_chroms(exonsBy(txdb,by=c("tx"),use.names=T),chr_map)
#   tx_by_gene<-map_chroms(transcriptsBy(txdb,by="gene"),chr_map)
#   five_by_tx<-map_chroms(fiveUTRsByTranscript(txdb,use.names=T),chr_map)
#   three_by_tx<-map_chroms(threeUTRsByTranscript(txdb,use.names=T),chr_map)
#   cds_by_tx <- map_chroms(cdsBy(txdb,by="tx",use.names=T),chr_map)
# }
all_genes<-genes(txdb)
exons_by_gene<-exonsBy(txdb,by="gene")
exons_by_tx<-exonsBy(txdb,by=c("tx"),use.names=T)
tx_by_gene<-transcriptsBy(txdb,by="gene")
five_by_tx<-fiveUTRsByTranscript(txdb,use.names=T)
three_by_tx<-threeUTRsByTranscript(txdb,use.names=T)
cds_by_tx <- cdsBy(txdb,by="tx",use.names=T)

#------------------------------------------------------------------------------#
# annotating junctions using reference introns

introns_by_tx<-data.frame(unlist(intronsByTranscript(txdb,use.names=T)))
ref_juncs<-sprintf("%s:%s:%s:%s",introns_by_tx$seqnames,introns_by_tx$start-1,introns_by_tx$end+1,introns_by_tx$strand)
target_juncs <- sprintf("%s:%s:%s:%s",introns$chr,introns$start,introns$end,introns$strand)

annotated_juncs <- target_juncs %in% ref_juncs
introns$ann <- "unannotated"
introns$ann[annotated_juncs]<-"annotated"
if ("verdict" %in% colnames(introns)){
  introns <- introns %>% dplyr::filter(verdict != "unknown_strand")
}

rm(introns_by_tx)

#------------------------------------------------------------------------------#
# performing transcript formation

data_canon<-data.frame()
sequences<-c()
cds_storage <- list()
intron_length<-nrow(introns)
for (i in seq(intron_length)){
  print(sprintf("%d introns out of %d total introns",i,intron_length))
  curr_introns<-introns[i,]
  ann<-curr_introns$ann
  target_junc <- unname(as.character(introns[i,seq(4)]))
  if (!(target_junc[4] %in% c("+","-"))){next}
  genes<-find_genes(target_junc,all_genes)
  all_tx<-extract_transcripts(genes, tx_by_gene) # the genes and transcripts associated with the junction

  if (length(genes$cis) != 0){
    # form cis and intragene trans-splcing events if they exist
    for (gene in genes$cis){
      # find the exons that are associated with the junction start and junction end, handles cryptic junctions
      exons<-choose_exons(target_junc,exons_by_gene,gene)
      if (length(exons$exons_start) == 0 | length(exons$exons_end) == 0){next}
      # find the transcripts per exon for start and end splice sites associated with each exon
      tx_starts_ends<-choose_transcripts(exons,all_tx$cis[[gene]],c(),exons_by_tx)
      # iterate through all exons for start and all exons for end
      for (start_exon in names(tx_starts_ends$starts)){
        for (end_exon in names(tx_starts_ends$ends)){
            # pair transcripts together
          start_trans<-tx_starts_ends$starts[[start_exon]]
          end_trans<-tx_starts_ends$ends[[end_exon]]
          trans_combos<-data.frame(expand.grid(start_trans,end_trans))
          names(trans_combos)<-c("start_trans","end_trans")
          trans_combos$start_trans<-as.character(trans_combos$start_trans)
          trans_combos$end_trans<-as.character(trans_combos$end_trans)
          trans_combos <- trans_combos %>% dplyr::filter(start_trans == end_trans)
          if (nrow(trans_combos)==0){next}
          trans_combos<-unique(trans_combos)
          # iterate through every combination and pair the transcript GRanges object together
          for (pair in seq(nrow(trans_combos))){
            protein_coding <- "No"
            data_canon_fill<-data.frame()

            trans_fir<-as.character(trans_combos[pair,1])
            trans_sec<-as.character(trans_combos[pair,2])
            trans_pair<-unique(c(trans_fir,trans_sec))
            start_exons<-sort(exons_by_tx[[trans_fir]])
            end_exons<-sort(exons_by_tx[[trans_sec]])

            start_exon_split<-str_split(start_exon,"[-]")[[1]]
            end_exon_split<-str_split(end_exon,"[-]")[[1]]

            start_loc<-which(start(ranges(start_exons))==start_exon_split[1]
                             & end(ranges(start_exons))==start_exon_split[2])

            end_loc<- which(start(ranges(end_exons))==end_exon_split[1]
                            & end(ranges(end_exons))==end_exon_split[2])

            combo_exons<-c(start_exons[1:start_loc],
                           end_exons[end_loc:length(end_exons)])
            tx_junc<-start_loc
            junc<-c(start_loc,start_loc+1)

            # create transcript CDS
            if (!(trans_fir %in% names(five_by_tx) & trans_fir %in% names(three_by_tx)) |
                !(trans_sec %in% names(three_by_tx) & trans_sec %in% names(five_by_tx))){
              cds_mod_info <- list()
              cds_mod_info[[1]] <- data.frame(combo_exons)
              cds_mod_info[[2]] <- junc
              cds_mod<-cds_mod_info[[1]]
              tx_junc_loc<-cds_mod_info[[2]]
              if (cds_mod$strand[1] == "-") { # not intuitive, but necessary for getSeq()
                cds_mod<-cds_mod[order(cds_mod[,2],decreasing=TRUE),]
              }
              sequence<-getSeq(bsgenome,makeGRangesFromDataFrame(cds_mod,keep.extra.columns=TRUE))
              junction<-paste(as.character(junc), collapse=",")
              sequ<-paste(as.character(sequence[1:length(sequence)]), collapse="")
              orf_dat<-find_orfs(sequ,tx_junc_loc)
              names(sequ)<-paste(c(curr_introns$chr,as.character(curr_introns$start),as.character(curr_introns$end),
                                   gene,paste(unique(trans_pair),collapse="-")),collapse=":")
              sequences<-c(sequences, sequ)
              mod<-mod_made(data.frame(combo_exons), cds_mod)
              if ("deltapsi" %in% colnames(curr_introns)){
                deltapsi<-curr_introns$deltapsi # deltapsi
              } else {
                deltapsi<-NA

              }
              if ("verdict" %in% colnames(curr_introns)){
                verdict<-curr_introns$verdict # verdict
              } else{
                verdict<-NA
              }
              if ("clusterID" %in% colnames(curr_introns)){
                clusterID<-curr_introns$clusterID # verdict
              } else{
                clusterID<-NA
              }
              next_row<-c(clusterID,
                          curr_introns$chr,
                          curr_introns$start,
                          curr_introns$end,
                          curr_introns$strand,
                          gene,
                          paste(unique(trans_pair),collapse="-"),
                          mod,
                          junction,
                          paste(tx_junc_loc,collapse=","),
                          orf_dat[4],
                          verdict,
                          deltapsi,
                          orf_dat[1],
                          orf_dat[3],
                          orf_dat[5],
                          orf_dat[6],
                          protein_coding,
                          start_exon,
                          end_exon,
                          ann)
            } else {
              protein_coding <- "Yes"
              if (any(as.character(strand(combo_exons)) == "-")) {
                UTR5<-sort(five_by_tx[[trans_sec]])
                UTR3<-sort(three_by_tx[[trans_fir]])
              } else {
                UTR5<-sort(five_by_tx[[trans_fir]])
                UTR3<-sort(three_by_tx[[trans_sec]])
              }
              junc<-is_UTR(target_junc,UTR5,UTR3,junc)
              if (str_detect(junc[1],"p") & str_detect(junc[2],"p") & junc[1]!=junc[2]){
                cds_mod_info <- list()
                cds_mod_info[[1]] <- data.frame(combo_exons)
                cds_mod_info[[2]] <- junc
                cds_mod<-cds_mod_info[[1]]
                tx_junc_loc<-cds_mod_info[[2]]
                if (cds_mod$strand[1] == "-") { # not intuitive, but necessary for getSeq()
                  cds_mod<-cds_mod[order(cds_mod[,2],decreasing=TRUE),]
                }
                sequence<-getSeq(bsgenome,makeGRangesFromDataFrame(cds_mod,keep.extra.columns=TRUE))
                junction<-paste(as.character(junc), collapse=",")
                sequ<-paste(as.character(sequence[1:length(sequence)]), collapse="")
                if (!check_tx(sequ)){protein_coding<-"tx_error"}
                orf_dat<-find_orfs(sequ,tx_junc_loc)
                names(sequ)<-paste(c(curr_introns$chr,as.character(curr_introns$start),as.character(curr_introns$end),
                                     gene,paste(unique(trans_pair),collapse="-")),collapse=":")
                sequences<-c(sequences, sequ)
                mod<-mod_made(data.frame(combo_cds), cds_mod)
              } else if (str_detect(junc[1],"p") & str_detect(junc[2],"p") & junc[1]==junc[2]) {
                combo_cds<-create_cds(combo_exons,UTR5,UTR3,tx_junc)
                cds_mod_info <- list()
                cds_mod_info[[1]] <- combo_cds
                cds_mod_info[[2]] <- junc
                cds_mod<-cds_mod_info[[1]]
                tx_junc_loc<-cds_mod_info[[2]]
                if (cds_mod$strand[1] == "-") { # not intuitive, but necessary for getSeq()
                  cds_mod<-cds_mod[order(cds_mod[,2],decreasing=TRUE),]
                }
                sequence<-getSeq(bsgenome,makeGRangesFromDataFrame(cds_mod,keep.extra.columns=TRUE))
                junction<-paste(as.character(junc), collapse=",")
                sequ<-paste(as.character(sequence[1:length(sequence)]), collapse="")
                orf_dat<-find_orfs(sequ,tx_junc_loc)
                names(sequ)<-paste(c(curr_introns$chr,as.character(curr_introns$start),as.character(curr_introns$end),
                                     gene,paste(unique(trans_pair),collapse="-")),collapse=":")
                sequences<-c(sequences,sequ)
                mod<-mod_made(data.frame(combo_cds), cds_mod)
              } else {
                combo_cds<-create_cds(combo_exons,UTR5,UTR3,tx_junc)
                cds_mod_info<-modify_cds(combo_cds, target_junc, junc, start_exon, end_exon)
                cds_mod<-cds_mod_info[[1]]
                tx_junc_loc<-cds_mod_info[[2]]
                if (cds_mod$strand[1] == "-") { # not intuitive, but necessary for getSeq()
                  cds_mod<-cds_mod[order(cds_mod[,2],decreasing=TRUE),]
                  combo_cds<-combo_cds[order(combo_cds[,2],decreasing=TRUE),]
                }
                sequence<-getSeq(bsgenome,makeGRangesFromDataFrame(cds_mod,keep.extra.columns=TRUE))
                junction<-paste(as.character(junc), collapse=",")
                sequ<-paste(as.character(sequence[1:length(sequence)]), collapse="")
                orf_dat<-find_orfs(sequ,tx_junc_loc)
                names(sequ)<-paste(c(curr_introns$chr,as.character(curr_introns$start),as.character(curr_introns$end),
                                     gene,paste(unique(trans_pair),collapse="-")),collapse=":")
                sequences<-c(sequences,sequ)
                mod<-mod_made(data.frame(combo_cds), cds_mod)
              }
              # modify the joined cds using the junction
              if ("deltapsi" %in% colnames(curr_introns)){
                deltapsi<-curr_introns$deltapsi # deltapsi
              } else {
                deltapsi<-NA

              }
              if ("verdict" %in% colnames(curr_introns)){
                verdict<-curr_introns$verdict # verdict
              } else{
                verdict<-NA
              }
              if ("clusterID" %in% colnames(curr_introns)){
                clusterID<-curr_introns$clusterID # verdict
              } else{
                clusterID<-NA
              }
              next_row<-c(clusterID,
                          curr_introns$chr,
                          curr_introns$start,
                          curr_introns$end,
                          curr_introns$strand,
                          gene,
                          paste(unique(trans_pair),collapse="-"),
                          mod,
                          junction,
                          paste(tx_junc_loc,collapse=","),
                          orf_dat[4],
                          verdict,
                          deltapsi,
                          orf_dat[1],
                          orf_dat[3],
                          orf_dat[5],
                          orf_dat[6],
                          protein_coding,
                          start_exon,
                          end_exon,
                          ann)
            }

            # for (r in seq(3)){
            # filling the data_canon dataframe with the junction information
            col_names<-c("cluster","chr","start","end","strand","gene","tx_id","modified","is_UTR",
                         "tx_junc_loc","pep_junc_loc","verdict","deltapsi",
                         "error","peptide","tx_length","start_stop","protein_coding",
                         "start_exon","end_exon","annotated") # introns cols are c("chr","start","end","gene","verdict")
            if (nrow(data_canon_fill)==0){
              data_canon_fill<-rbind(data_canon_fill,next_row)
              colnames(data_canon_fill)<-col_names
            } else {
              data_canon_fill<-rbind(data_canon_fill,next_row)
            }
            # }
            data_canon_fill<-unique(data_canon_fill)
            data_canon<-rbind(data_canon,data_canon_fill)
            cds_mod_dat <- list(cds_mod$start,cds_mod$end)
            names(cds_mod_dat) <- c("start","end")
            cds_storage<-list.append(cds_storage,cds_mod_dat)          }
        }
      }
    }
  }

  gene_combos<-data.frame(expand.grid(genes$trans$start,genes$trans$end))
  colnames(gene_combos)<-c("start_gene","end_gene")
  gene_combos$start_gene<-as.character(gene_combos$start_gene)
  gene_combos$end_gene<-as.character(gene_combos$end_gene)
  gene_combos<-gene_combos %>% dplyr::filter(start_gene != end_gene)
  if (length(genes$trans$start) != 0 & length(genes$trans$end) != 0){
    for (pair in seq(nrow(gene_combos))){
      gene_pair<-as.character(unname(gene_combos[pair,]))
      # find the exons that are associated with the junction start and junction end, handles cryptic junctions
      exons<-choose_exons(target_junc,exons_by_gene,gene_pair)
      if (length(exons$exons_start) == 0 | length(exons$exons_start == 0)){next}
      # find the transcripts per exon for start and end splice sites associated with each exon
      tx_starts_ends<-choose_transcripts(exons,all_tx$trans$start[[gene_pair[1]]],
                                         all_tx$trans$end[[gene_pair[2]]],exons_by_tx)
      # iterate through all exons for start and all exons for end
      for (start_exon in names(tx_starts_ends$starts)){
        for (end_exon in names(tx_starts_ends$ends)){
          # pair transcripts together
          start_trans<-tx_starts_ends$starts[[start_exon]]
          end_trans<-tx_starts_ends$ends[[end_exon]]
          trans_combos<-data.frame(expand.grid(start_trans,end_trans))
          names(trans_combos)<-c("start_trans","end_trans")
          trans_combos$start_trans<-as.character(trans_combos$start_trans)
          trans_combos$end_trans<-as.character(trans_combos$end_trans)
          trans_combos <- trans_combos %>% dplyr::filter(start_trans != end_trans)
          if (nrow(trans_combos)==0){next}
          trans_combos<-unique(trans_combos)
          # iterate through every combination and pair the transcript GRanges object together
          for (pair in seq(nrow(trans_combos))){
            protein_coding <- "No"
            data_canon_fill<-data.frame()

            trans_fir<-as.character(trans_combos[pair,1])
            trans_sec<-as.character(trans_combos[pair,2])
            trans_pair<-unique(c(trans_fir,trans_sec))
            start_exons<-sort(exons_by_tx[[trans_fir]])
            end_exons<-sort(exons_by_tx[[trans_sec]])

            start_exon_split<-str_split(start_exon,"[-]")[[1]]
            end_exon_split<-str_split(end_exon,"[-]")[[1]]

            start_loc<-which(start(ranges(start_exons))==start_exon_split[1]
                             & end(ranges(start_exons))==start_exon_split[2])

            end_loc<- which(start(ranges(end_exons))==end_exon_split[1]
                            & end(ranges(end_exons))==end_exon_split[2])

            combo_exons<-c(start_exons[1:start_loc],
                           end_exons[end_loc:length(end_exons)])
            tx_junc<-start_loc
            junc<-c(start_loc,start_loc+1)

            # create transcript CDS
            if (!(trans_fir %in% names(five_by_tx) & trans_fir %in% names(three_by_tx)) |
                !(trans_sec %in% names(three_by_tx) & trans_sec %in% names(five_by_tx))){
              cds_mod_info <- list()
              cds_mod_info[[1]] <- data.frame(combo_exons)
              cds_mod_info[[2]] <- junc
              cds_mod<-cds_mod_info[[1]]
              tx_junc_loc<-cds_mod_info[[2]]
              if (cds_mod$strand[1] == "-") { # not intuitive, but necessary for getSeq()
                cds_mod<-cds_mod[order(cds_mod[,2],decreasing=TRUE),]
              }
              sequence<-getSeq(bsgenome,makeGRangesFromDataFrame(cds_mod,keep.extra.columns=TRUE))
              junction<-paste(as.character(junc), collapse=",")
              sequ<-paste(as.character(sequence[1:length(sequence)]), collapse="")
              orf_dat<-find_orfs(sequ,tx_junc_loc)
              names(sequ)<-paste(c(curr_introns$chr,as.character(curr_introns$start),as.character(curr_introns$end),
                                   paste(gene_pair,collapse="-"),paste(unique(trans_pair),collapse="-")),collapse=":")
              sequences<-c(sequences, sequ)
              mod<-mod_made(data.frame(combo_exons), cds_mod)
              if ("deltapsi" %in% colnames(curr_introns)){
                deltapsi<-curr_introns$deltapsi # deltapsi
              } else {
                deltapsi<-NA

              }
              if ("verdict" %in% colnames(curr_introns)){
                verdict<-curr_introns$verdict # verdict
              } else{
                verdict<-NA
              }
              if ("clusterID" %in% colnames(curr_introns)){
                clusterID<-curr_introns$clusterID # verdict
              } else{
                clusterID<-NA
              }
              next_row<-c(clusterID,
                          curr_introns$chr,
                          curr_introns$start,
                          curr_introns$end,
                          curr_introns$strand,
                          paste(gene_pair,collapse="-"),
                          paste(unique(trans_pair),collapse="-"),
                          mod,
                          junction,
                          paste(tx_junc_loc,collapse=","),
                          orf_dat[4],
                          verdict,
                          deltapsi,
                          orf_dat[1],
                          orf_dat[3],
                          orf_dat[5],
                          orf_dat[6],
                          protein_coding,
                          start_exon,
                          end_exon,
                          ann)
            } else {
              protein_coding <- "Yes"
              if (any(as.character(strand(combo_exons)) == "-")) {
                UTR5<-sort(five_by_tx[[trans_sec]])
                UTR3<-sort(three_by_tx[[trans_fir]])
              } else {
                UTR5<-sort(five_by_tx[[trans_fir]])
                UTR3<-sort(three_by_tx[[trans_sec]])
              }
              junc<-is_UTR(target_junc,UTR5,UTR3,junc)
              if (str_detect(junc[1],"p") & str_detect(junc[2],"p") & junc[1]!=junc[2]){
                cds_mod_info <- list()
                cds_mod_info[[1]] <- data.frame(combo_exons)
                cds_mod_info[[2]] <- junc
                cds_mod<-cds_mod_info[[1]]
                tx_junc_loc<-cds_mod_info[[2]]
                if (cds_mod$strand[1] == "-") { # not intuitive, but necessary for getSeq()
                  cds_mod<-cds_mod[order(cds_mod[,2],decreasing=TRUE),]
                }
                sequence<-getSeq(bsgenome,makeGRangesFromDataFrame(cds_mod,keep.extra.columns=TRUE))
                junction<-paste(as.character(junc), collapse=",")
                sequ<-paste(as.character(sequence[1:length(sequence)]), collapse="")
                orf_dat<-find_orfs(sequ,tx_junc_loc)
                names(sequ)<-paste(c(curr_introns$chr,as.character(curr_introns$start),as.character(curr_introns$end),
                                     paste(gene_pair,collapse="-"),paste(unique(trans_pair),collapse="-")),collapse=":")
                sequences<-c(sequences, sequ)
                mod<-mod_made(data.frame(combo_cds), cds_mod)
              } else if (str_detect(junc[1],"p") & str_detect(junc[2],"p") & junc[1]==junc[2]) {
                combo_cds<-create_cds(combo_exons,UTR5,UTR3,tx_junc)
                cds_mod_info <- list()
                cds_mod_info[[1]] <- combo_cds
                cds_mod_info[[2]] <- junc
                cds_mod<-cds_mod_info[[1]]
                tx_junc_loc<-cds_mod_info[[2]]
                if (cds_mod$strand[1] == "-") { # not intuitive, but necessary for getSeq()
                  cds_mod<-cds_mod[order(cds_mod[,2],decreasing=TRUE),]
                }
                sequence<-getSeq(bsgenome,makeGRangesFromDataFrame(cds_mod,keep.extra.columns=TRUE))
                junction<-paste(as.character(junc), collapse=",")
                sequ<-paste(as.character(sequence[1:length(sequence)]), collapse="")
                orf_dat<-find_orfs(sequ,tx_junc_loc)
                names(sequ)<-paste(c(curr_introns$chr,as.character(curr_introns$start),as.character(curr_introns$end),
                                     paste(gene_pair,collapse="-"),paste(unique(trans_pair),collapse="-")),collapse=":")
                sequences<-c(sequences,sequ)
                mod<-mod_made(data.frame(combo_cds), cds_mod)
              } else {
                combo_cds<-create_cds(combo_exons,UTR5,UTR3,tx_junc)
                cds_mod_info<-modify_cds(combo_cds, target_junc, junc, start_exon, end_exon)
                cds_mod<-cds_mod_info[[1]]
                tx_junc_loc<-cds_mod_info[[2]]
                if (cds_mod$strand[1] == "-") { # not intuitive, but necessary for getSeq()
                  cds_mod<-cds_mod[order(cds_mod[,2],decreasing=TRUE),]
                  combo_cds<-combo_cds[order(combo_cds[,2],decreasing=TRUE),]
                }
                sequence<-getSeq(bsgenome,makeGRangesFromDataFrame(cds_mod,keep.extra.columns=TRUE))
                junction<-paste(as.character(junc), collapse=",")
                sequ<-paste(as.character(sequence[1:length(sequence)]), collapse="")
                if (!check_tx(sequ)){protein_coding<-"tx_error"}
                orf_dat<-find_orfs(sequ,tx_junc_loc)
                names(sequ)<-paste(c(curr_introns$chr,as.character(curr_introns$start),as.character(curr_introns$end),
                                     paste(gene_pair,collapse="-"),paste(unique(trans_pair),collapse="-")),collapse=":")
                sequences<-c(sequences,sequ)
                mod<-mod_made(data.frame(combo_cds), cds_mod)
              }
              # modify the joined cds using the junction
              if ("deltapsi" %in% colnames(curr_introns)){
                deltapsi<-curr_introns$deltapsi # deltapsi
              } else {
                deltapsi<-NA
              }
              if ("verdict" %in% colnames(curr_introns)){
                verdict<-curr_introns$verdict # verdict
              } else{
                verdict<-NA
              }
              if ("clusterID" %in% colnames(curr_introns)){
                clusterID<-curr_introns$clusterID # verdict
              } else{
                clusterID<-NA
              }
              next_row<-c(clusterID,
                          curr_introns$chr,
                          curr_introns$start,
                          curr_introns$end,
                          curr_introns$strand,
                          paste(gene_pair,collapse="-"),
                          paste(unique(trans_pair),collapse="-"),
                          mod,
                          junction,
                          paste(tx_junc_loc,collapse=","),
                          orf_dat[4],
                          verdict,
                          deltapsi,
                          orf_dat[1],
                          orf_dat[3],
                          orf_dat[5],
                          orf_dat[6],
                          protein_coding,
                          start_exon,
                          end_exon,
                          ann)
            }

            # filling the data_canon dataframe with the junction information
            col_names<-c("cluster","chr","start","end","strand","gene","tx_id","modified","is_UTR",
                         "tx_junc_loc","pep_junc_loc","verdict","deltapsi",
                         "error","peptide","tx_length","start_stop","protein_coding",
                         "start_exon","end_exon","annotated") # introns cols are c("chr","start","end","gene","verdict")
            if (nrow(data_canon_fill)==0){
              data_canon_fill<-rbind(data_canon_fill,next_row)
              colnames(data_canon_fill)<-col_names
            } else {
              data_canon_fill<-rbind(data_canon_fill,next_row)
            }
            data_canon_fill<-unique(data_canon_fill)
            data_canon<-rbind(data_canon,data_canon_fill)
            cds_mod_dat <- list(cds_mod$start,cds_mod$end)
            names(cds_mod_dat) <- c("start","end")
            cds_storage<-list.append(cds_storage,cds_mod_dat)
          }
        }
      }
    }
  }
}

#------------------------------------------------------------------------------------------------------------------------------------------------#
# saving data

out<-sprintf("%s_%s%s",out_prefix,"data_splicemutr",".txt")
out_rds <- sprintf("%s_%s%s",out_prefix,"data_splicemutr",".rds")
out_fasta<-sprintf("%s_%s%s",out_prefix,"sequences",".fa")
out_cds <- sprintf("%s_cds_stored.rds",out_prefix)
saveRDS(data_canon,file=out_rds)
write.table(data_canon,file=out,col.names=T,row.names=F,quote=F)
sequences<-DNAStringSet(sequences)
writeXStringSet(sequences,out_fasta)
saveRDS(cds_storage,file=out_cds)
