#hello
#' @name find_genes
#' @title find_genes
#' @param target_junc the target junction c(chr, strand, start, end)
#' @param all_genes the GenomicRanges::GRanges object of genes
#' @return  list of cis, trans_start, and trans_end genes associated with the junction
#' @export
#'
find_genes<-function(target_junc,all_genes){

  junc_start<-GenomicRanges::GRanges(seqnames = c(target_junc[1]), strand = target_junc[4],
                      ranges = IRanges::IRanges(start = as.numeric(target_junc[2]), width = c(1)))
  junc_end<-GenomicRanges::GRanges(seqnames = target_junc[1], strand = target_junc[4],
                    ranges = IRanges::IRanges(start = as.numeric(target_junc[3]), width = c(1)))
  start<-GenomicRanges::findOverlaps(all_genes,junc_start)
  end<-GenomicRanges::findOverlaps(all_genes,junc_end)
  hits_start<-S4Vectors::queryHits(start)
  hits_end<-S4Vectors::queryHits(end)
  genes_start<-all_genes[hits_start]
  genes_end<-all_genes[hits_end]
  genes<-list()
  genes[["cis"]]<-names(genes_start[which(genes_start %in% genes_end)])
  genes[["trans"]][["start"]]<-names(genes_start[which(!(genes_start %in% genes_end))])
  genes[["trans"]][["end"]]<-names(genes_end[which(!(genes_end %in% genes_start))])
  return(genes)
}

#' @name extract_transcripts
#' @title extract_transcripts
#' @param genes character vector containing genes the target junction can exist in
#' @param tx_by_gene GenomicRanges::GRanges object containing introns per for a specified transcript
#' @return  list containing genes and transcripts
#' @export
#'
extract_transcripts<-function(genes,tx_by_gene){

  genes_cis<-genes$cis
  genes_trans_start<-genes$trans$start
  genes_trans_end<-genes$trans$end

  gene_tx_cis<-lapply(genes_cis, function(gene){
    tx_by_gene[[gene]]@elementMetadata$tx_name
  })
  names(gene_tx_cis)<-genes_cis

  gene_tx_trans_start<-lapply(genes_trans_start, function(gene){
    tx_by_gene[[gene]]@elementMetadata$tx_name
  })
  names(gene_tx_trans_start)<-genes_trans_start

  gene_tx_trans_end<-lapply(genes_trans_end, function(gene){
    tx_by_gene[[gene]]@elementMetadata$tx_name
  })
  names(gene_tx_trans_end)<-genes_trans_end
  gene_tx_list<-list()
  gene_tx_list[["cis"]]<-gene_tx_cis
  gene_tx_list[["trans"]][["start"]]<-gene_tx_trans_start
  gene_tx_list[["trans"]][["end"]]<-gene_tx_trans_end

  return(gene_tx_list)
}

#' @name find_introns
#' @title find_introns
#' @param target_junc target junction character vector as c(chr,strand,start,stop)
#' @param tx_introns GenomicRanges::GRanges object containing introns per for a specified transcript
#' @param exact logical indicating whether the junction should exactly match a reference intron or not
#' @return  chracter vector indicating the intron range the target junction overlaps with
#' @export
#'
find_introns<-function(target_junc,tx_introns,exact=FALSE){
  # finds the introns per transcript that the target junction overlaps with
  # inputs:
  #   target_junction: target junction character vector as c(chr,strand,start,stop)
  #   tx_introns: GenomicRanges::GRanges object containing introns per for a specified transcript
  # outputs:
  #   hits: chracter vector indicating the intron range the target junction overlaps with
  start_adj<-as.numeric(target_junc[3])+1
  end_adj<-as.numeric(target_junc[4])-1
  if (exact){
    junction<-GenomicRanges::GRanges(seqnames = c(target_junc[1]), strand = target_junc[2],
                      ranges = IRanges::IRanges(start = start_adj, width = end_adj-start_adj+1))
    return(any(tx_introns==junction))
  } else {
    junction<-GenomicRanges::GRanges(seqnames = c(target_junc[1]), strand = target_junc[2],
                      ranges = IRanges::IRanges(start = start_adj, width = end_adj-start_adj+1))
    junc_intr<-GenomicRanges::findOverlaps(junction,tx_introns)
    hits<-S4Vectors::subjectHits(junc_intr)
    return(hits)
  }
}

#' @name which_junctions_overlap
#' @title which_junctions_overlap
#' @param curr_introns a list of junctions corresponding to the current gene being analyzed
#' @return  a list of non-overlapping junctions for the gene with row indeces dictating the associated row index of the curr_introns list
#' @export
#'
which_junctions_overlap<-function(curr_introns){

  overlaps<-GenomicRanges::findOverlaps(GenomicRanges::makeGRangesFromDataFrame(curr_introns,keep.extra.columns=TRUE))
  queries<-S4Vectors::queryHits(overlaps)
  subjects<-S4Vectors::subjectHits(overlaps)
  queries_all<-seq(1,nrow(curr_introns))
  shared_sets<-lapply(queries_all,function(x){
    query_test<-!(queries_all %in% subjects[queries==x])
    if (all(query_test==FALSE)){
      0
    }else {
      queries_all[!(queries_all %in% subjects[queries==x])]
    }
  })
  return(shared_sets)
}

#' @name adjust_width
#' @title adjust_width
#' @param cds A cds dataframe with vectors of numeric widths, ends, and starts
#' @return  The cds data with an adjusted width
#' @export
#'
adjust_width<-function(cds){
  # adjusts the width column of a cds dataframe
  # inputs:
  #   cds: a dataframe containing transcript coding sequence coordinates
  # outputs:
  #   cds: a dataframe containing transcript coding sequence coordinates with recalc width column
  for (i in seq_len(nrow(cds))){
    cds[i,"width"]<- (cds[i,"end"]-cds[i,"start"] + 1)
  }
  return(cds)
}

#' @name running_sum
#' @title running_sum
#' @param vector The vector of numbers
#' @return  the vector modified to reflect the incremented sum
#' @export
#'
running_sum<-function(vector){
  # translates genomic relative coordinates to transcript relative coordinates
  # inputs:
  #   vector: a vector of genomic relative coordinates
  # outputs:
  #   vector: a vector of transcript relative coordinates

  for (i in 1:length(vector)){
    if (i != 1){
      vector[i]<-vector[i-1]+vector[i]
    }
  }
  return(vector)
}

#' @name modify_cds
#' @title modify_cds
#' @param combo_cds the combined cds output from create_cds()
#' @param target_junc the intron coordinates used for modifying the cds
#' @param junc the junction location information output from is_UTR
#' @param start_exon start exon (relative to junction) coordinates as character in form start-end
#' @param end_exon end exon (relative to junction) coordinates as character in form start-end
#' @return  a list containing the modified cds and the junction location in the cds cds_mod junc_loc
#' @export
#'
modify_cds <- function(combo_cds, target_junc, junc, start_exon, end_exon){
  # are either junction in a UTR?
  junc_UTR <- which(stringr::str_detect(junc,"p"))

  # is the junction inside either exon region?
  start_pos<-as.numeric(target_junc[2])
  end_pos<-as.numeric(target_junc[3])

  start_exon_split <- as.numeric(stringr::str_split(start_exon,"-")[[1]])
  in_start <- start_pos < start_exon_split[2]
  start_exon_mod <- c(start_exon_split[1],start_pos)

  end_exon_split <- as.numeric(stringr::str_split(end_exon,"-")[[1]])
  in_end <- end_pos > end_exon_split[1]
  end_exon_mod <- c(end_pos,end_exon_split[2])

  cds<-data.frame(combo_cds)

  locations<-c(NA,NA)
  if (length(junc_UTR)==0){
    locations <- c(NA,NA)
    for (i in seq(nrow(cds))){
      if (start_pos >= cds[i,"start"] & start_pos <= cds[i,"end"]){
        locations[1] <- i
      } else if (i < nrow(cds)){
        if (start_pos > cds[i,"end"] & start_pos < cds[i+1,"start"]){
          locations[1] <- i
        }
      }
      if (end_pos >= cds[i,"start"] & end_pos <= cds[i,"end"]){
        locations[2]<-i
      } else if (i > 1){
        if (end_pos > cds[i-1,"end"] & end_pos < cds[i,"start"]){
          locations[2]<-i
        }
      }
    }
    cds_mod <- cds
    cds_mod[locations[1],"end"] <- start_exon_mod[2]
    cds_mod[locations[2],"start"] <- end_exon_mod[1]
    junc_start <- locations[1]
    junc_end <- locations[2]
  } else if (junc_UTR == 1){
    # find location of end exon in cds
    for (i in seq(nrow(cds))){
      if (end_pos >= cds[i,"start"] & end_pos <= cds[i,"end"]){
        locations[2]<-i
      } else if (i > 1){
        if (end_pos > cds[i-1,"end"] & end_pos < cds[i,"start"]){
          locations[2]<-i
        }
      }
    }
    cds_mod <- cds[seq(locations[2],nrow(cds)),]
    row <- cds_mod[1,]
    row$start <- start_exon_mod[1]
    row$end <- start_exon_mod[2]
    cds_mod <- rbind(row,cds_mod)
    cds_mod[,c("exon_id","exon_name","exon_rank")]<-NA
    junc_start <- 1
    junc_end <- 2
  } else { # junc_UTR == 2
    # find location of start exon in cds
    for (i in seq(nrow(cds))){
      if (start_pos >= cds[i,"start"] & start_pos <= cds[i,"end"]){
        locations[1] <- i
      } else if (i < nrow(cds)){
        if (start_pos > cds[i,"end"] & start_pos < cds[i+1,"start"]){
          locations[1] <- i
        }
      }
    }
    cds_mod <- cds[seq(locations[1]),]
    row <- cds_mod[1,]
    row$start <- end_exon_mod[1]
    row$end <- end_exon_mod[2]
    junc_start <- nrow(cds_mod)
    cds_mod <- rbind(cds_mod,row)
    cds_mod[,c("exon_id","exon_name","exon_rank")]<-NA
    junc_end <- junc_start+1
  }
  cds_mod<-adjust_width(cds_mod)
  width_sum<-running_sum(cds_mod$width)
  junc_loc<-c(width_sum[junc_start],width_sum[junc_start]+1)
  return(list(cds_mod,junc_loc))
}

#' @name mod_made
#' @title mod_made
#' @param cds reference cds as sorted dataframe
#' @param cds_mod junction-modified cds as sorted dataframe
#' @return "changed" if the starts and ends are different, "normal" if the starts and ends are the same, "add_reg" if more elements in cds_mod than cds, sub_reg if less elements in cds_mod than cds
#' @export
#'
mod_made<-function(cds,cds_mod){
  # determines whether there is a difference between two cds dataframes
  # inputs:
  #   cds: reference cds
  #   cds_mod: changed cds
  # output:
  #   out: "changed" if the starts and ends are different, "normal" if the starts and ends are the same, "add_reg" if an split exon or orf overwrite (UTR)

  if(nrow(cds)==nrow(cds_mod)){
    cds[,c("exon_id","exon_name","exon_rank")]<-1
    cds_mod[,c("exon_id","exon_name","exon_rank")]<-1
    if(any((cds==cds_mod)==FALSE)==TRUE){
      out<-"changed"
      return(out)
    }else{
      out<-"normal"
      return(out)
    }
  }else if (nrow(cds)>nrow(cds_mod)){
    out<-"sub_reg"
    return(out)
  } else {
    out<-"add_reg"
    return(out)
  }
}

#' @name which_codon
#' @title which_codon
#' @param cod DNAstring containing codon start and end coordinates
#' @param tx_junc_loc The junction location in the transcript
#' @param start_loc start codon location in peptide coordinates relative to whole transcript
#' @param stop_loc stop codon location in peptide coordinates relative to whole transcript
#' @return junction location in translated protein coordinates or tx_junc_loc if in UTR
#' @export
#'
which_codon<-function(cod,tx_junc_loc,start_loc,stop_loc){

  starts<-BiocGenerics::start(cod)
  ends<-BiocGenerics::end(cod)
  pep_junc_pos<-which(starts %in% tx_junc_loc | ends %in% tx_junc_loc)

  if(length(pep_junc_pos)==0){
    pep_junc_pos<-tx_junc_loc-start_loc-1 # -1 to account for codon length
    if (all(pep_junc_pos < 0)){
      missing <- 3 - (abs(pep_junc_pos) %% 3)
      pep_junc_pos<-(pep_junc_pos-missing)/3
    } else {
      pep_junc_pos<-tx_junc_loc-stop_loc-2 # -2 to account for codon length
      missing <- 3 - ((pep_junc_pos)%%3)
      pep_junc_pos<-(pep_junc_pos+missing)/3
      pep_junc_pos<-sprintf("+%d,+%d",pep_junc_pos[1],pep_junc_pos[2])
    }
  }
  return(pep_junc_pos)
}

#' @name junc_loc_pep
#' @title junc_loc_pep
#' @param tx_junc_loc The junction location in the transcript
#' @param tx The transcript sequence
#' @param start_loc start codon location in peptide coordinates relative to whole transcript
#' @param stop_loc stop codon location in peptide coordinates relative to whole transcript
#' @return junction location in translated protein coordinates or tx_junc_loc if in UTR
#'
junc_loc_pep<-function(tx_junc_loc,tx,start_loc,stop_loc){ #write
  # determines the junction location in translated protein coordinates
  # inputs:
  #   tx_junc_loc: junction location in transcript coordinates
  #   tx: transcript sequence
  # outputs:
  # pep_junc_loc: junction location in translated protein coordinates or tx_junc_loc if in UTR
  if(any(is.na(tx_junc_loc))){return(tx_junc_loc)}
  if (length(tx_junc_loc)==1 | any(stringr::str_detect(tx_junc_loc,"p"))){
    return(tx_junc_loc)
  } else {
    if ( ((nchar(tx)) %%3 > 0) | ((tx_junc_loc[1]+1) != tx_junc_loc[2])){
      return(F)
    } else {
      cod<-Biostrings::codons(Biostrings::DNAString(tx))
    }
    codon<-which_codon(cod,tx_junc_loc,start_loc,stop_loc)
    return(codon)
  }
}

#' @name find_orfs
#' @title find_orfs
#' @param sequ The DNA sequence in character form to analyze
#' @param tx_junc_loc The junction location in the transcript
#' @return a vector of error type, protein if translatable or DNA sequence, ensembl transcript name , junction location relative to peptide start, modified transcript length, start and stop codon locations
#' @export
#'
find_orfs<-function(sequ,tx_junc_loc){
  start_codon<-"ATG"
  stop_codons<-c("TAG","TAA","TGA")
  errors <- NA
  sequ_len <- nchar(sequ)
  x<-1
  if (x >= nchar(sequ)){
    return(c("seq_too_small",NA,NA,NA,NA,NA,NA,sequ_len))
  } else {
    sequence_split<-stringi::stri_sub(str = sequ,
                             from = seq(x, nchar(sequ) - 3 + 1,by = 3),
                             length = 3)

    start_cods<-which(sequence_split %in% start_codon)
    # no translation
    if (length(start_cods)==0){
      stop_cods<-which(sequence_split %in% stop_codons)
      if (length(stop_cods) == 0){
        return(c("no_starts:no_stops",sequ,NA,NA,NA,NA))
      } else if (stop_cods[1] != length(sequence_split)){
        orf_loc<-c(NA, 1 + (stop_cods[1]*3)-1 -1+x)
        return(c("no_starts:stop_not_end",sequ,NA,NA,sequ_len,paste(c(orf_loc[1],orf_loc[2]),collapse=":")))
      } else {
        orf_loc<-c(NA, 1 + (stop_cods[1]*3)-1 -1+x)
        return(c("no_starts",sequ,NA,NA,sequ_len,paste(c(orf_loc[1],orf_loc[2]),collapse=":")))
      }
    }
    sequence_split<-sequence_split[start_cods[1]:length(sequence_split)]
    stop_cods<-which(sequence_split %in% stop_codons)
    if (length(stop_cods)==0){
      if (start_cods[1] != 1){
        orf_loc<-c((start_cods[1]*3)-2-1+x, NA)
        orf_seq<-paste(sequence_split[1:length(sequence_split)],collapse="")
        return(c("starts_not_beg:no_stops",sequ,NA,NA,sequ_len,paste(c(orf_loc[1],orf_loc[2]),collapse=":")))
      } else {
        orf_loc<-c((start_cods[1]*3)-2-1+x, NA)
        return(c("no_stops",sequ,NA,NA,sequ_len,paste(c(orf_loc[1],orf_loc[2]),collapse=":")))
      }
    }
    if (start_cods[1] != 1 & stop_cods[1] != length(sequence_split)){
      error <- "start_not_beg:stop_not_end"
    } else if (start_cods[1] == 1 & stop_cods[1] != length(sequence_split)){
      error <- "stop_not_end"
  } else if (start_cods[1] != 1 & stop_cods[1] == length(sequence_split)){
    error <- "start_not_beg"
  } else {
    error <- "tx"
  }
  orf_seq<-paste(sequence_split[1:stop_cods[1]],collapse="")
  seq_str_set<-Biostrings::DNAStringSet(orf_seq)
  trans<-Biostrings::translate(seq_str_set,no.init.codon = F)
  orf_loc<-c((start_cods[1]*3)-2, (start_cods[1]*3)-2 + (stop_cods[1]*3)-1)-1+x
  pep_junc_loc<-junc_loc_pep(tx_junc_loc,orf_seq,orf_loc[1],orf_loc[2])
  out<-c(error,
         orf_seq,
         as.character(trans),
         paste(pep_junc_loc,collapse=","),
         sequ_len,
         paste(c(orf_loc[1],orf_loc[2]),collapse=":"))
  return(out)
  }
}

#' @name check_tx
#' @title check_tx
#' @param seq The DNA sequence in character form to analyze
#' @return TRUE or FALSE dependent upon whether the transcript has start codon at first codon and stop codon at last codon
#' @export
#'
check_tx<-function(seq){
  # This function takes as input a DNA sequence, and tests to see whether it contains a start and stop codon.
  # inputs:
  #   seq: the DNA string sequence
  # outputs:
  #   out: logical saying if both starts and stop codons exist in sequence

  if ( all(unique(strsplit(seq,"")[[1]]) %in% c("A","T","C","G")) ) {
    seq<-Biostrings::DNAStringSet(substr(seq,1,nchar(seq)-((nchar(seq)-1+1) %% 3)))
    trans<-Biostrings::translate(seq,no.init.codon = FALSE)
  } else {
    trans<-seq
  }
  out<-FALSE
  a<-strsplit(as.character(trans),"")
  stop<-which((a[[1]] %in% "*") == TRUE)
  start<-which((a[[1]] %in% "M") == TRUE)
  translatable <- (length(start)>0) & (length(stop)>0) & (start[1]==1) & (stop[1]==length(a[[1]])) & (length(stop)==1)
  if (translatable){out<-TRUE}
  return(out)
}

#' @name choose_exons
#' @title choose_exons
#' @param target_junc the target junction to use for analysis
#' @param exons_by_gene the exons per gene in the genome
#' @param gene_tar a vector pair of gene targets, left is first or only gene, right is second gene. Single element vector if only one gene.
#' @return the unique start and end junction exons
#' @export
#'
choose_exons<-function(target_junc, exons_by_gene, gene_tar){
  junc_start<-GenomicRanges::GRanges(seqnames = c(target_junc[1]), strand = target_junc[4],
                      ranges = IRanges::IRanges(start = as.numeric(target_junc[2]), width = c(1)))
  junc_end<-GenomicRanges::GRanges(seqnames = target_junc[1], strand = target_junc[4],
                    ranges = IRanges::IRanges(start = as.numeric(target_junc[3]), width = c(1)))
  if (length(gene_tar)==1){
    exons_small<-sort(exons_by_gene[[gene_tar]])
    start<-GenomicRanges::findOverlaps(exons_small,junc_start)
    hits_start<-S4Vectors::queryHits(start)
    if (length(start)==0){
      hits_start<-which(IRanges::ranges(junc_start) > IRanges::ranges(exons_small))
      hits_start<-hits_start[length(hits_start)]
    }
    end<-GenomicRanges::findOverlaps(exons_small,junc_end)
    hits_end<-S4Vectors::queryHits(end)
    if (length(end)==0){
      hits_end<-which(IRanges::ranges(junc_end) < IRanges::ranges(exons_small))[1]
    }
    exons_start<-unique(exons_small[hits_start])
    exons_end<-unique(exons_small[hits_end])
  } else {
    exons_small_start<-exons_by_gene[[gene_tar[1]]]
    start<-GenomicRanges::findOverlaps(exons_small_start,junc_start)
    hits_start<-S4Vectors::queryHits(start)
    if (length(start)==0){
      hits_start<-which(IRanges::ranges(junc_start) > IRanges::ranges(exons_small_start))
      hits_start<-hits_start[length(hits_start)]
    }
    exons_small_end<-exons_by_gene[[gene_tar[2]]]
    end<-GenomicRanges::findOverlaps(exons_small_end,junc_end)
    hits_end<-S4Vectors::queryHits(end)
    if (length(end)==0){
      hits_end<-which(IRanges::ranges(junc_end) < IRanges::ranges(exons_small_end))[1]
    }
    exons_start<-unique(exons_small_start[hits_start])
    exons_end<-unique(exons_small_end[hits_end])
  }

  exons_start<-exons_start[which(!(exons_start %in% exons_end))]
  exons_end<-exons_end[which(!(exons_end %in% exons_start))]

  # exons_small<-data.frame(unique(sort(exons_by_gene[[gene]])))
  # exons_small<-unique(exons_small[,1:3])
  exons_start_end<-list(exons_start,exons_end)
  names(exons_start_end)<-c("exons_start","exons_end")
  return(exons_start_end)
}

#' @name choose_transcripts
#' @title choose_transcripts
#' @param exons the exons for the start and end of the junction
#' @param transcripts_1 the character vector of transcripts for the start gene
#' @param transcripts_2 the character vector of transcripts for the end gene
#' @param exons_by_tx the exons per transcript
#' @return a list of start junction and end junction transcripts that can be paired
#' @export
#'
choose_transcripts<-function(exons, transcripts_1, transcripts_2=c(), exons_by_tx){
  transcripts<-transcripts_1
  starts<-lapply(seq(length(exons$exons_start)),function(exon){
    vapply(transcripts,function(trans){
      if (exons$exons_start[exon] %in% exons_by_tx[[trans]]){
        return(T)
      } else {
        return(F)
      }
    },logical(1))})

  tx_starts<-lapply(starts,function(start){
    transcripts[start]
  })
  start<-BiocGenerics::start(IRanges::ranges(exons$exons_start))
  end<-BiocGenerics::end(IRanges::ranges(exons$exons_start))
  names(tx_starts)<-sprintf("%d-%d",start,end)
  if (length(transcripts_2)>0){
    transcripts<-transcripts_2
  }
  ends<-lapply(seq(length(exons$exons_end)),function(exon){
    vapply(transcripts,function(trans){
      if (exons$exons_end[exon] %in% exons_by_tx[[trans]]){
        return(T)
      } else {
        return(F)
      }
    },logical(1))})
  tx_ends<-lapply(ends,function(end){
    transcripts[end]
  })
  start<-BiocGenerics::start(IRanges::ranges(exons$exons_end))
  end<-BiocGenerics::end(IRanges::ranges(exons$exons_end))
  names(tx_ends)<-sprintf("%d-%d",start,end)
  tx_starts_ends<-list(tx_starts,tx_ends)
  names(tx_starts_ends)<-c("starts","ends")
  return(tx_starts_ends)
}

#' @name create_cds
#' @title create_cds
#' @param combo_exons the combined exons for the two transcripts associated with the junction
#' @param UTR5 the 5' UTR for the first transcript
#' @param UTR3 the 3' UTR for the final transcript
#' @param tx_junc the junction block coordinates in the transcript
#' @return created cds
#' @export
#'
create_cds<-function(combo_exons,UTR5,UTR3,tx_junc){
  combo_start<-combo_exons[1:tx_junc]
  comb_start<-data.frame(combo_start)
  combo_end<-combo_exons[(tx_junc+1):length(combo_exons)]
  comb_end<-data.frame(combo_end)
  start_UTR<-UTR5
  end_UTR<-UTR3
  if (any(as.character(BiocGenerics::strand(combo_exons)) == "-")){
    start_UTR<-UTR3
    end_UTR<-UTR5
  } else {
    start_UTR<-UTR5
    end_UTR<-UTR3
  }

  start_pos<-BiocGenerics::end(IRanges::ranges(start_UTR[length(start_UTR)]))+1
  end_pos<-BiocGenerics::start(IRanges::ranges(end_UTR[1]))-1

  # need to modify this section such that UTRs can eb found to precede and be found after
  # the overlaps are not working because of this
  # do what you did in choose_exons
  start_cds_block<-GenomicRanges::findOverlaps(combo_start,start_UTR[length(start_UTR)])
  if (length(start_cds_block)==0){
    start_UTR_in<-T
    hits_start<-tx_junc
  } else {
    start_UTR_in<-F
    hits_start<-S4Vectors::queryHits(start_cds_block[length(start_cds_block)])
  }
  end_cds_block<-GenomicRanges::findOverlaps(combo_end,end_UTR[1])
  if (length(end_cds_block)==0){
    end_UTR_in<-T
    hits_end<-1
  } else {
    end_UTR_in<-F
    hits_end<-S4Vectors::queryHits(end_cds_block[1])
  }

  if (!start_UTR_in){ # subtracting the UTR from the cds
    comb_start<-DataCombine::InsertRow(comb_start,comb_start[hits_start,],hits_start)
    comb_start[hits_start+1,"start"]<-start_pos
    comb_start<-comb_start[(hits_start+1):nrow(comb_start),]
  } else {
    comb_start<-comb_start[1:hits_start,]
  }

  if (!end_UTR_in){
    comb_end<-DataCombine::InsertRow(comb_end,comb_end[hits_end,],hits_end)
    comb_end[(hits_end),"end"]<-end_pos
    comb_end<-comb_end[1:(hits_end),]
  } else {
    comb_end<-comb_end[hits_end:nrow(comb_end),]
  }
  comb<-rbind(comb_start,comb_end)
  comb<-adjust_width(comb)
  return(comb)
}

#' @name is_UTR
#' @title is_UTR
#' @param target_junc the target junction in form chr,strand,start,end
#' @param UTR5 the 5' UTR
#' @param UTR3 the 3' UTR
#' @param junc the original junc call
#' @return two logicals indicating whether the junction is UTR or not
#' @export
#'
is_UTR<-function(target_junc, UTR5, UTR3, junc){
  if (target_junc[4] == "-"){
    flip=T
    start_UTR<-UTR3
    end_UTR<-UTR5
  } else {
    flip=F
    start_UTR<-UTR5
    end_UTR<-UTR3
  }
  start_UTR_start<-BiocGenerics::start(IRanges::ranges(start_UTR[1]))
  start_UTR_end<-BiocGenerics::end(IRanges::ranges(start_UTR[length(start_UTR)]))
  end_UTR_start<-BiocGenerics::start(IRanges::ranges(end_UTR[1]))
  end_UTR_end<-BiocGenerics::end(IRanges::ranges(end_UTR[length(end_UTR)]))
  junc_start<-as.numeric(target_junc[2])
  junc_end<-as.numeric(target_junc[3])
  junc_1<-junc[1]
  junc_2<-junc[2]
  if (junc_start >= start_UTR_start & junc_start <= start_UTR_end){
    if (flip){
      junc_1<-"3p"
    } else {
      junc_1<-"5p"
    }
  } else if (junc_start >= end_UTR_start & junc_start <= end_UTR_end){
    if (flip){
      junc_1<-"5p"
    } else {
      junc_1<-"3p"
    }
  }
  if (junc_end >= start_UTR_start & junc_end <= start_UTR_end){
    if (flip){
      junc_2<-"3p"
    } else {
      junc_2<-"5p"
    }
  } else if (junc_end >= end_UTR_start & junc_end <= end_UTR_end){
    if (flip){
      junc_2<-"5p"
    } else {
      junc_2<-"3p"
    }
  }
  return (c(junc_1,junc_2))
}

#' @name countKmers
#' @title countKmers
#' @param seqs_vec vector of sequences (i.e. data_canon[,2])
#' @param K the kmer length
#' @return a list of integer vectors with kmer character vectors as names
#' @export
#'
countKmers <- function(seqs_vec, K){

  kmers <- lapply(seqs_vec, function(x){
    table(stringi::stri_sub(str = x,
                   from = seq(1, nchar(x) - K + 1,
                              by = 1),
                   length = K))
  })

  return(kmers) # return is not necessary in R, but I like to make it explicit
}

#' @name getUniqKmers
#' @title getUniqKmers
#' @param kmers character vector of non-unique kmers
#' @return vector of character vectors representing unique list of kmers from the union of all sequences
#'
getUniqKmers <- function(kmers){

  return(unique(unlist(lapply(kmers, names))))
}

#' @name split_kmers
#' @title split_kmers
#' @param kmers The kmers character vector
#' @param save_dir The directory to save the kmer vector to.
#' @return A vector detailing the number of kmers saved per file saved.
#'
split_kmers<-function(kmers,save_dir){
  kmers_test<-c()
  len<-length(kmers)
  mil<-2000000
  vapply(seq_len(ceiling(length(kmers)/mil)),function(x){
    a<-seq_len(mil)
    if (x > 1){
      a<-a + mil*(x-1)
      if (max(a) > len){
        a<-c(a[1]:len)
      }
    }
    kmers_test[1:length(a)]<-kmers[a]
    kmers_test<-data.frame(kmers_test)
    utils::write.table(kmers_test,file=sprintf(save_dir,x),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
    return(nrow(kmers_test))
  },numeric(1))
}

#' @name star_to_leaf
#' @title star_to_leaf
#' @param STAR_file The star file to be converted
#' @return converted junc_dat dataframe
#' @export
#'
star_to_leaf <- function(STAR_file){
  star_dat <- utils::read.table(STAR_file,header=F)
  chrom <- star_dat[,1]
  chromStart <- as.numeric(star_dat[,2]-1)
  chromEnd <- as.numeric(star_dat[,3]+1)
  name <- sprintf("%s%d",rep("JUNC",nrow(star_dat)),seq(nrow(star_dat)))
  score <- as.numeric(star_dat[,7])
  strand <- vapply(star_dat[,4],function(val){
    if (val == 0){
      return("*")
    } else if (val == 1){
      return("+")
    } else {
      return("-")
    }
  },character(1))
  junc_dat <- data.frame(chrom,chromStart,chromEnd,name,score,strand)
  return(junc_dat)
}

#' @name format_introns
#' @title format_introns
#' @param introns the introns file
#' @return splicemutr formatted introns file
#' @export
#'
format_introns <- function(introns){
  if ("clusterID" %in% colnames(introns)){
    strand <- as.character(matrix(unlist(str_split(introns$clusterID,"_")),byrow=T,nrow=nrow(introns))[,3])
    introns$strand <- strand
  }
  intron_names<-which(colnames(introns)=="chrom")
  if(length(intron_names)>0){
    colnames(introns)[intron_names]<-"chr"
  }
  chr_mask <- colnames(introns) %in% c("chr")
  start_mask <- colnames(introns)=="start"
  end_mask <- colnames(introns)=="end"
  strand_mask <- colnames(introns)=="strand"
  important_info <- chr_mask | start_mask | end_mask | strand_mask
  introns_fill <- introns[,important_info]
  introns_fill <- cbind(introns_fill,introns[,!important_info])
}

