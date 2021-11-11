#!/usr/bin/env Rscript

library(dplyr)
library(stringr)
library(optparse)
library(GO.db)
library(KEGGREST)
library(KEGG.db)
library(fgsea)
library(org.Hs.eg.db)

#------------------------------------------------------------------------------#
# handling command line input

arguments <- parse_args(OptionParser(usage = "",
                         description="",
                         option_list=list(
                           make_option(c("-d","--diff_dat_file"),
                                       default = "",
                                       help="diff_dat"),
                           make_option(c("-o","--out"),
                                       default = "",
                                       help="the out file"))))
opt=arguments
diff_dat_file <- opt$diff_dat_file
out<-opt$out

#------------------------------------------------------------------------------#
# loading in necessary data

diff_dat_all <- readRDS(diff_dat_file)

#------------------------------------------------------------------------------#
# running fgsea

fgsea_all <- list()
for (comp in names(diff_dat_all)){
  print(comp)
  entrez_to_ensembl <- data.frame(org.Hs.egENSEMBL)
  diff_dat_all$entrez <- unlist(lapply(rownames(diff_dat_all),function(gene){
    a<-which(entrez_to_ensembl$ensembl_id == gene)
    if(length(a)==0){
      return(NA)
    } else {
      return(entrez_to_ensembl[a[1],"gene_id"])
    }
  }))
  
  ENSEMBLToGO <- mapIds(org.Hs.eg.db,keys=rownames(diff_dat_all),column='GO',keytype = 'ENSEMBL',multiVals = list)
  GOToENSEMBL <- sapply(reverseSplit(ENSEMBLToGO),unique)
  GOToENSEMBL <- GOToENSEMBL[sapply(GOToENSEMBL,length)>5]
  
  ENSEMBLToKEGG <- mapIds(org.Hs.eg.db,keys=rownames(diff_dat_all),column='PATH',keytype = 'ENSEMBL',multiVals = list)
  KEGGToENSEMBL <- sapply(reverseSplit(ENSEMBLToKEGG),unique)
  KEGGToENSEMBL <- KEGGToENSEMBL[sapply(KEGGToENSEMBL,length)>5]
  
  diff_dat_filt <- diff_dat_all %>% dplyr::filter(!is.na(log2FoldChange) & !is.na(pvalue))
  ranks <- diff_dat_filt$log2FoldChange*-log10(diff_dat_filt$pvalue)
  names(ranks)<-rownames(diff_dat_filt)
  ranks<-ranks[unname(!is.na(ranks))]
  ranks_pos<-ranks[ranks>=0]
  ranks_neg<-ranks[ranks<0]
  ranks_pos<-rank(ranks_pos,ties.method = "random")
  ranks_neg<--1*rank(-1*ranks_neg,ties.method = "random")
  ranks<-c(ranks_pos,ranks_neg)
  fgseaRes_GO <- fgsea(GOToENSEMBL, ranks, eps=0,scoreType="std")
  fgseaRes_GO<-cbind(sapply(c('ONTOLOGY','TERM','DEFINITION'),
                            function(x){mapIds(GO.db,keys=fgseaRes_GO$pathway,keytype='GOID',column=x)}),fgseaRes_GO)
  fgseaRes_KEGG <- fgsea(KEGGToENSEMBL, ranks, eps=0,scoreType="std")
  fgseaRes_KEGG <- cbind(KEGG.Name=unlist(as.list(KEGGPATHID2NAME)[fgseaRes_KEGG$pathway]),
                         fgseaRes_KEGG)
  
  GO_table <- fgseaRes_GO %>% dplyr::filter(padj <= 0.05)
  GO_table<-GO_table[,c("TERM","DEFINITION","padj","NES","size")]
  KEGG_table <- fgseaRes_KEGG %>% dplyr::filter(padj <= 0.05)
  KEGG_table<-KEGG_table[,c("KEGG.Name","padj","NES","size")]
  
  a<-list()
  a[["GO"]]<-GO_table
  a[["KEGG"]]<-KEGG_table
  fgsea_all[[comp]]<-a
}

#------------------------------------------------------------------------------#
# saving fgsea data

saveRDS(fgsea_all,file=out)
