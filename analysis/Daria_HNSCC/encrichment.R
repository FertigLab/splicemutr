
library(stringr)
library(msigdbr)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(knitr)
library(ggplot2)
library(dplyr)

load("/media/theron/My_Passport/head_and_neck_DARIA/data/leafcutter_10_21_2020/data.Rdata")

genes<-unique(vapply(clusters$gene,function(x){
  str_replace(str_replace(x,"<i>",""),"</i>","")
},character(1)))

# FDR<-clusters$FDR
# annotation<-clusters$annotation
# dat<-data.frame(genes,annotation,FDR)
# colnames(dat)<-c("genes","annotation","FDR")

#------------------------------------------------------------------------------#
# performing enrichment analysis

hallmark<-msigdbr(species="Homo sapiens",category="H")
types<-unique(hallmark$gs_name)
genes_symbol<-hallmark$gene_symbol[which(hallmark$gene_symbol %in% genes)]
samps_types<-hallmark$gs_name[which(hallmark$gene_symbol %in% genes)]
vals_origin<-vapply(types,function(y){
  length(which(samps_types==y))
},numeric(1))

genes_samps<-lapply(unique(samps_types),function(samp){
  genes_samp<-genes_symbol[which(samps_types==samp)]
})
names(genes_samps)<-unique(samps_types)

genes_comp<-lapply(seq(length(genes_samps)),function(num){
  genes_num<-genes_samps[[num]]
  comp<-vapply(seq(length(genes_samps)),function(inner_num){
    length(which(genes_num %in% genes_samps[[inner_num]]))/length(genes_num)
  },numeric(1))
})
names(genes_comp)<-names(genes_samps)

comp_dat <- data.frame(matrix(unlist(genes_comp), nrow=length(genes_comp), byrow=TRUE))
colnames(comp_dat)<-names(genes_samps)
rownames(comp_dat)<-names(genes_samps)

dat_pvals<-read.table("/media/theron/My_Passport/head_and_neck_DARIA/data/leafcutter_10_21_2020/enriched_genes.txt",
                      header=T)
pathways<-dat_pvals[seq(26),1]
comp_dat_sig<-comp_dat[pathways,pathways]

pheatmap::pheatmap(comp_dat_sig,cluster_rows = F,cluster_cols = F,
                   show_rownames = T,fontsize_col = 7,fontsize_row=7,
                   clustering_method="complete",display_numbers=T)

#------------------------------------------------------------------------------#
# tsimulating for significance

txdb<-TxDb.Hsapiens.UCSC.hg19.knownGene
genes_entrez<-as.numeric(names(genes(txdb)))

types<-unique(hallmark$gs_name)
sample_size<-nrow(clusters)
vals<-lapply(seq(10000),function(x){
  samps<-sample(genes_entrez,sample_size)
  samps_types<-hallmark$gs_name[which(hallmark$entrez_gene %in% samps)]
  vals<-vapply(types,function(y){
    length(which(samps_types==y))
  },numeric(1))
  names(vals)<-types
  return(vals)
})

vals_frame <- data.frame(matrix(unlist(vals), nrow=length(vals), byrow=TRUE))
colnames(vals_frame)<-types

#------------------------------------------------------------------------------#
# testing for significance

pvals<-vapply(types,function(x){
  test_val<-unname(vals_origin[x])
  length(which(vals_frame[,x]>=test_val))/length(vals_frame[,x])
},numeric(1))
names(pvals)<-types

dat_pvals<-data.frame(names(pvals))
dat_pvals$pvals<-pvals
dat_pvals$counts<-vals_origin
colnames(dat_pvals)<-c("pathways","pvals","counts")
dat_pvals<-dat_pvals[order(dat_pvals$pvals),]

write.table(dat_pvals,
            sep='\t',row.names = F, col.names = T,quote=F,
            file="/media/theron/My_Passport/head_and_neck_DARIA/data/leafcutter_10_21_2020/enriched_genes.txt")

vals_to_plot<-as.data.frame(table(vals_frame$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION))
vals_to_plot$Var1<-as.numeric(as.vector(vals_to_plot$Var1))
vals_to_plot$Freq<-as.numeric(as.vector(vals_to_plot$Freq))

ggplot(vals_to_plot, aes(x=Var1))+
  geom_col(aes(y=Freq))+xlab("Counts")+
  labs(title="HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")+
  geom_vline(xintercept = vals_origin["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"])


vals_to_plot<-as.data.frame(table(vals_frame$HALLMARK_PI3K_AKT_MTOR_SIGNALING))
vals_to_plot$Var1<-as.numeric(as.vector(vals_to_plot$Var1))
vals_to_plot$Freq<-as.numeric(as.vector(vals_to_plot$Freq))

ggplot(vals_to_plot, aes(x=Var1))+
  geom_col(aes(y=Freq))+xlab("Counts")+
  labs(title="HALLMARK_PI3K_AKT_MTOR_SIGNALING")+
  geom_vline(xintercept = vals_origin["HALLMARK_PI3K_AKT_MTOR_SIGNALING"])


vals_to_plot<-as.data.frame(table(vals_frame$HALLMARK_TNFA_SIGNALING_VIA_NFKB))
vals_to_plot$Var1<-as.numeric(as.vector(vals_to_plot$Var1))
vals_to_plot$Freq<-as.numeric(as.vector(vals_to_plot$Freq))

ggplot(vals_to_plot, aes(x=Var1))+
  geom_col(aes(y=Freq))+xlab("Counts")+
  labs(title="HALLMARK_TNFA_SIGNALING_VIA_NFKB")+
  geom_vline(xintercept = vals_origin["HALLMARK_TNFA_SIGNALING_VIA_NFKB"])


vals_to_plot<-as.data.frame(table(vals_frame$HALLMARK_COMPLEMENT))
vals_to_plot$Var1<-as.numeric(as.vector(vals_to_plot$Var1))
vals_to_plot$Freq<-as.numeric(as.vector(vals_to_plot$Freq))

ggplot(vals_to_plot, aes(x=Var1))+
  geom_col(aes(y=Freq))+xlab("Counts")+
  labs(title="HALLMARK_COMPLEMENT")+
  geom_vline(xintercept = vals_origin["HALLMARK_COMPLEMENT"])

