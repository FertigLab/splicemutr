
#title: "Daria Class 1 Analysis"
#author: "Theron Palmer"
#date: "6/4/2021"

library(ggplot2)
library(stringr)
library(pheatmap)
library(dplyr)
library(ggpubr)
library(msigdbr)

#------------------------------------------------------------------------------#
#The directories where all of the data is

genotypes_file <- "/media/theron/My_Passport/TCGA_junctions/ext_dat/OptiTypeCallsHLA_20171207.tsv"
arcas_files <- "/media/theron/My_Passport/head_and_neck_DARIA/data/splicemutr_05_26_2021/processed_data/arcashla_dirs.txt"
counts_psi_file <- read.table("/media/theron/My_Passport/head_and_neck_DARIA/data/splicemutr_05_26_2021/processed_data/data_perind.counts")
HLA_files <- "/media/theron/My_Passport/TCGA_junctions/summary_files/%s_tx_dict_summary_perc.txt"
splicemutr_dat <- read.table("/media/theron/My_Passport/TCGA_junctions/formed_transcripts/data_splicemutr.txt", sep = " ",header=T)
groups_file <- read.table("/media/theron/My_Passport/head_and_neck_DARIA/data/splicemutr_05_26_2021/processed_data/groups_file.txt")

mutation_count_file <- "/media/theron/My_Passport/TCGA_junctions/cbioportal_data/Mutation_Count.txt"
mutation_counts <- read.table(mutation_count_file, sep="\t",header=T)

imm_dat <- read.table("/media/theron/My_Passport/TCGA_junctions/summary_files/HLA_summ.txt")
rownames(imm_dat)<-imm_dat[,1]

#------------------------------------------------------------------------------#
#Assigning filenames to genotypes

genotypes <- read.table(genotypes_file,sep=",",header=T)
genotypes$A1 <- unname(vapply(genotypes$A1,function(HLA){
  HLA_split<-str_split(HLA,"[*:]")[[1]]
  sprintf("HLA-%s%s-%s",HLA_split[1],HLA_split[2],HLA_split[3])
},character(1)))
genotypes$A2 <- unname(vapply(genotypes$A2,function(HLA){
  HLA_split<-str_split(HLA,"[*:]")[[1]]
  sprintf("HLA-%s%s-%s",HLA_split[1],HLA_split[2],HLA_split[3])
},character(1)))
genotypes$B1 <- unname(vapply(genotypes$B1,function(HLA){
  HLA_split<-str_split(HLA,"[*:]")[[1]]
  sprintf("HLA-%s%s-%s",HLA_split[1],HLA_split[2],HLA_split[3])
},character(1)))
genotypes$B2 <- unname(vapply(genotypes$B2,function(HLA){
  HLA_split<-str_split(HLA,"[*:]")[[1]]
  sprintf("HLA-%s%s-%s",HLA_split[1],HLA_split[2],HLA_split[3])
},character(1)))
genotypes$C1 <- unname(vapply(genotypes$C1,function(HLA){
  HLA_split<-str_split(HLA,"[*:]")[[1]]
  sprintf("HLA-%s%s-%s",HLA_split[1],HLA_split[2],HLA_split[3])
},character(1)))
genotypes$C2 <- unname(vapply(genotypes$C2,function(HLA){
  HLA_split<-str_split(HLA,"[*:]")[[1]]
  sprintf("HLA-%s%s-%s",HLA_split[1],HLA_split[2],HLA_split[3])
},character(1)))


#------------------------------------------------------------------------------#
# matching TCGA samples between mutation counts and genotypes

mutation_IDs <- vapply(genotypes$aliquot_id,function(ID){
  id_parts <- str_split(ID,"[-]")[[1]]
  sprintf("%s-%s-%s",id_parts[1],id_parts[2],id_parts[3])
},character(1))
genotypes$mutation_IDs <- mutation_IDs

mut_counts <- vapply(genotypes$mutation_IDs,function(ID){
  counts <<- mutation_counts[mutation_counts$V2==ID,"V4"]
  if (length(counts)==0){
    counts <- -1
  } else {
    counts <- max(counts)
  }
  return(counts)
},numeric(1))

cancer_type <- unname(vapply(genotypes$mutation_IDs,function(ID){
  cancer <- str_split(mutation_counts[mutation_counts$V2==ID,"V1"],"[_]")
  if (length(cancer)==0){
    cancer <- "NONE"
  } else {
    cancer <- cancer[[1]][1]
  }
},character(1)))

tum_or_norm <- unname(vapply(genotypes$aliquot_id,function(ID){
  type <- str_split(ID,"[-]")[[1]][4]
  type<-as.numeric(substr(type,1,2))
  if (type >= 1 & type <= 9){
    return("T")
  } else if (type >= 10 & type <= 19) {
    return("N")
  }
  return("C")
},character(1)))
genotypes$mut_counts <- mut_counts
genotypes$cancer_type <- cancer_type
genotypes$tum_or_norm <- tum_or_norm

genotypes <- genotypes %>% dplyr::filter(cancer_type != "NONE")
cancer_types <- unique(genotypes$cancer_type)
acc_cancers <- function(cancer_file){
  cancers <- read.table(cancer_file)
  cancers <- as.character(cancers[,1])
  cancers <- vapply(cancers,function(cancer){
    tolower(basename(cancer))
  },character(1))
}
 acceptable_cancers <-acc_cancers("/media/theron/My_Passport/TCGA_junctions/TCGA_cancers/filenames.txt") 
 
 genotypes$acceptable <- vapply(genotypes$cancer_type,function(cancer){
   if (cancer %in% acceptable_cancers){
     return("processed")
   } else {
     return("unprocessed")
   }
 },character(1))
 
#------------------------------------------------------------------------------#
# contructing binding per sample and genotype

get_score_2<-function(HLA){
  return(imm_dat[HLA,2])
}
get_score_3<-function(HLA){
  return(imm_dat[HLA,3])
}
get_score_4<-function(HLA){
  return(imm_dat[HLA,4])
}

genotypes$imm_score_median <- apply(apply(genotypes[,seq(6)],2,get_score_2),1,sum)
genotypes$imm_score_mean <- apply(apply(genotypes[,seq(6)],2,get_score_3),1,sum)
genotypes$imm_score_num <- apply(apply(genotypes[,seq(6)],2,get_score_4),1,sum)

#------------------------------------------------------------------------------#
# developing cancer breakdown


cancer_breakdown <- data.frame(t(vapply(cancer_types,function(c_type){
    genotype_cancer <- genotypes %>% dplyr::filter(cancer_type == c_type)
    return(c(mean(genotype_cancer$mut_counts), mean(genotype_cancer$imm_score_num)))
},numeric(2))))
cancer_breakdown$cancer <- cancer_types
colnames(cancer_breakdown)<-c("mut_counts","imm_score_num","cancer_type")


#------------------------------------------------------------------------------#
# plotting genotype data

pdf(file="/media/theron/My_Passport/TCGA_junctions/binders.pdf",width=10, height=10)

print(ggplot(cancer_breakdown,aes(x=log10(imm_score_num),y=log10(mut_counts),label=cancer_type))+
        geom_text()+stat_cor(method = "spearman"))

genotypes_tum <- genotypes %>% dplyr::filter(tum_or_norm == "T")
genotypes_norm <- genotypes %>% dplyr::filter(tum_or_norm == "N")
print(ggplot(genotypes_tum,aes(x=imm_score_median,y=log10(mut_counts),color=acceptable))+
  geom_point())
print(ggplot(genotypes_tum,aes(x=imm_score_mean,y=log10(mut_counts),color=cancer_type))+
  geom_point())
print(ggplot(genotypes_tum,aes(x=log10(imm_score_num),y=log10(mut_counts),color=cancer_type))+
  geom_point())

print(ggplot(genotypes_norm,aes(x=imm_score_median,y=log10(mut_counts),color=cancer_type))+
  geom_point())
print(ggplot(genotypes_norm,aes(x=imm_score_mean,y=log10(mut_counts),color=cancer_type))+
  geom_point())
print(ggplot(genotypes_norm,aes(x=log10(imm_score_num),y=log10(mut_counts),color=cancer_type))+
  geom_point())

print(ggplot(genotypes,aes(x=cancer_type,y=imm_score_median))+geom_boxplot())
print(ggplot(genotypes,aes(x=cancer_type,y=imm_score_mean))+geom_boxplot())
print(ggplot(genotypes,aes(x=cancer_type,y=log10(imm_score_num)))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
  
dev.off()


#------------------------------------------------------------------------------#
# plotting splicemutr data

# pdf(file="/media/theron/My_Passport/head_and_neck_DARIA/data/splicemutr_05_26_2021/processed_data/class1/binders.pdf",width=10, height=10)

print(ggplot(tumor_metadata,aes(deltapsi,y=tx_length,color=verdict))+geom_point())

ggplot(cluster_metrics,aes(x=tum_ov_norm))+geom_histogram(binwidth=0.5)+scale_y_log10()

# png("/media/theron/My_Passport/SARE_2021/project/spliemutr_images/tumor_vs_normal_2.png",res=200)
ggplot(tumor_metadata,aes(x=deltapsi,y=log2((tumor_mean+1)/(normal_mean+1)),color=verdict))+
  geom_point(size=2)+
  # stat_cor(method = "spearman")+
  ylab("log2(mean tumor SB counts/mean normal SB counts)")+
  theme(text = element_text(size = 15))
ggsave("/media/theron/My_Passport/SARE_2021/project/spliemutr_images/tumor_vs_normal_2.png",
       width=8, height = 8, units = "in",dpi=300)
# dev.off()

print(ggplot(tumor_metadata,aes(x=deltapsi,y=log2((tumor_mean_scale_counts+1)/(normal_mean_scale_counts+1)),color=verdict))+
        geom_point(shape=1)+
        stat_cor(method = "spearman")+ylab("log2(mean log2(tumor SB counts)/mean log(normal SB counts)) scaled"))
print(ggplot(tumor_metadata,aes(x=deltapsi,y=((tumor_mean_scale_counts+1)),color=verdict))+
        geom_point(shape=1)+
        stat_cor(method = "spearman")+ylab("mean log(tumor SB counts) scaled)"))
print(ggplot(tumor_metadata,aes(x=deltapsi,y=((normal_mean_scale_counts+1)),color=verdict))+
        geom_point(shape=1)+
        stat_cor(method = "spearman")+ylab("log10(mean normal SB counts scaled)"))
print(ggplot(tumor_metadata,aes(x=deltapsi,y=log2((tumor_mean_numerator+1)/(normal_mean_numerator+1)),color=verdict))+
        geom_point(shape=1)+
        stat_cor(method = "spearman")+ylab("log2(mean tumor junc counts/mean normal junc counts)"))
print(ggplot(tumor_metadata,aes(x=deltapsi,y=log2((tumor_mean_numerator+1)/(normal_mean_numerator+1))))+
        geom_point(shape=1)+
        stat_cor(method = "spearman")+geom_smooth(method=lm)+
        ylab("log2(mean tumor junc counts/mean normal junc counts)"))
print(ggplot(tumor_metadata,aes(x=deltapsi,y=log2((tumor_mean_scale+1)/(normal_mean_scale+1)),color=verdict))+
        geom_point(shape=1)+
        stat_cor(method = "spearman")+
        ylab("log2(mean tumor junc counts/mean normal junc counts) scaled"))
print(ggplot(tumor_metadata,aes(x=deltapsi,y=log2((tumor_mean_scale+1)),color=verdict))+
        geom_point(shape=1)+
        stat_cor(method = "spearman")+
        ylab("log2(mean tumor junc counts) scaled"))
print(ggplot(tumor_metadata,aes(x=deltapsi,y=log2((normal_mean_scale+1)),color=verdict))+
        geom_point(shape=1)+
        stat_cor(method = "spearman")+
        ylab("log2(mean normal junc counts) scaled"))
print(ggplot(tumor_metadata,aes(x=deltapsi,y=log2((tumor_mean_scale+1)/(normal_mean_scale+1))))+
        geom_point(shape=1)+
        stat_cor(method = "spearman")+geom_smooth(method=lm)+
        ylab("log2(mean tumor junc counts/mean normal junc counts) scaled"))

print(ggplot(tumor_metadata_nonann,aes(x=deltapsi,y=log2((tumor_mean+1)/(normal_mean+1)),color=verdict))+
        geom_point(shape=1)+
      stat_cor(method = "spearman")+
        ylab("log2(mean tumor SB counts/mean normal SB counts)"))
print(ggplot(tumor_metadata_nonann,aes(x=deltapsi,y=log2((tumor_mean_scale_counts+1)/(normal_mean_scale_counts+1)),color=verdict))+
        geom_point(shape=1)+
        stat_cor(method = "spearman")+
        ylab("log2(mean tumor SB counts/mean normal SB counts) scaled"))
print(ggplot(tumor_metadata_nonann,aes(x=deltapsi,y=log2((tumor_mean_numerator+1)/(normal_mean_numerator+1)),color=verdict))+
        geom_point(shape=1)+
        stat_cor(method = "spearman")+
        ylab("log2(mean tumor junc counts/mean normal junc counts)"))
print(ggplot(tumor_metadata_nonann,aes(x=deltapsi,y=log2((tumor_mean_numerator+1)/(normal_mean_numerator+1))))+
        geom_point(shape=1)+
        stat_cor(method = "spearman")+geom_smooth(method=lm)+
        ylab("log2(mean tumor junc counts/mean normal junc counts)"))
print(ggplot(tumor_metadata_nonann,aes(x=deltapsi,y=log2((tumor_mean_scale+1)/(normal_mean_scale+1)),color=verdict))+
        geom_point(shape=1)+
        stat_cor(method = "spearman")+
        ylab("log2(mean tumor junc counts/mean normal junc counts) scaled"))
print(ggplot(tumor_metadata_nonann,aes(x=deltapsi,y=log2((tumor_mean_scale+1)/(normal_mean_scale+1))))+
        geom_point(shape=1)+
        stat_cor(method = "spearman")+geom_smooth(method=lm)+
        ylab("log2(mean tumor junc counts/mean normal junc counts) scaled"))

print(ggplot(sample_df,aes(x=sample,y=SB_count))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylab("log2(tumor SB counts/mean normal SB_counts")+
    xlab("Tumor Sample"))

print(ggplot(sample_df_tum,aes(x=sample,y=SB_count))+
        geom_boxplot()+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
        ylab("log2(tumor SB counts/mean normal SB_counts")+
        xlab("Tumor Sample"))


my_comparisons <- list( c("annotated", "cryptic_fiveprime"),
                        c("annotated", "cryptic_threeprime"),
                        c("annotated", "cryptic_unanchored"),
                        c("annotated", "novel annotated pair"),
                        c("cryptic_fiveprime", "cryptic_threeprime"),
                        c("cryptic_fiveprime", "cryptic_unanchored"),
                        c("cryptic_fiveprime", "novel annotated pair"),
                        c("cryptic_threeprime", "cryptic_unanchored"),
                        c("cryptic_threeprime", "novel annotated pair"),
                        c("cryptic_unanchored","novel annotated pair"))
print(ggplot(tumor_metadata_correct,aes(x=verdict,y=log2((tumor_mean_scale_counts+1)/(normal_mean_scale_counts+1))))+
        geom_boxplot()+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
        ylab("log2(mean tumor SB counts/mean normal SB_counts")+
        xlab("Verdict")+
        stat_compare_means(comparisons = my_comparisons,method="wilcox.test"))    # Add global p-value)

print(ggplot(tumor_metadata_correct,aes(x=verdict,y=log2((tumor_mean_scale+1)/(normal_mean_scale+1))))+
        geom_boxplot()+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
        ylab("log2(mean tumor junc counts/mean normal junc counts) scaled")+
        xlab("Verdict")+
        stat_compare_means(comparisons = my_comparisons,method="wilcox.test"))    # Add global p-value)

my_comparisons <- list( c("start_not_beg", "start_not_beg:stop_not_end"),
                        c("start_not_beg", "stop_not_end"),
                        c("start_not_beg", "tx"),
                        c("start_not_beg:stop_not_end", "stop_not_end"),
                        c("start_not_beg:stop_not_end", "tx"),
                        c("stop_not_end", "tx"))
tumor_metadata_filt <- tumor_metadata_correct %>% dplyr::filter(!(error %in% c("no_starts","no_starts:stop_not_end","no_stops","starts_not_beg:no_stops")))
print(ggplot(tumor_metadata_filt,aes(x=error,y=log2((tumor_mean_scale_counts+1)/(normal_mean_scale_counts+1))))+
        geom_boxplot()+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
        ylab("log2(mean tumor SB counts/mean normal SB_counts")+
        xlab("Errors")+
        stat_compare_means(comparisons = my_comparisons,method="wilcox.test"))
print(ggplot(tumor_metadata_filt,aes(x=error,y=log2((tumor_mean_scale+1)/(normal_mean_scale+1))))+
        geom_boxplot()+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
        ylab("log2(mean tumor junc counts/mean normal junc counts) scaled")+
        xlab("Errors")+
        stat_compare_means(comparisons = my_comparisons,method="wilcox.test"))


tumor_metadata_filt <- tumor_metadata_correct %>% dplyr::filter(deltapsi > 0)

my_comparisons <- list( c("annotated", "cryptic_fiveprime"),
                        c("annotated", "cryptic_threeprime"),
                        c("annotated", "cryptic_unanchored"),
                        c("annotated", "novel annotated pair"),
                        c("cryptic_fiveprime", "cryptic_threeprime"),
                        c("cryptic_fiveprime", "cryptic_unanchored"),
                        c("cryptic_fiveprime", "novel annotated pair"),
                        c("cryptic_threeprime", "cryptic_unanchored"),
                        c("cryptic_threeprime", "novel annotated pair"),
                        c("cryptic_unanchored","novel annotated pair"))

ggplot(tumor_metadata_filt,aes(x=verdict,y=log2((tumor_mean_scale_counts+1)/(normal_mean_scale_counts+1))))+
        geom_boxplot()+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
        ylab("log2(mean tumor SB counts/mean normal SB_counts)")+
        xlab("Verdict")+
        stat_compare_means(comparisons = my_comparisons,method="wilcox.test")+
  theme(text = element_text(size = 15))
# Add global p-value)
ggsave("/media/theron/My_Passport/SARE_2021/project/spliemutr_images/binders_by_junc_type_3.png",
       width=7, height = 7, units = "in",dpi=300)

print(ggplot(tumor_metadata_filt,aes(x=verdict,y=log2((tumor_mean_scale+1)/(normal_mean_scale+1))))+
        geom_boxplot()+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
        ylab("log2(mean tumor junc counts/mean normal junc counts) scaled Tumor Junctions")+
        xlab("Verdict")+
        stat_compare_means(comparisons = my_comparisons,method="wilcox.test"))

my_comparisons <- list( c("start_not_beg", "start_not_beg:stop_not_end"),
                        c("start_not_beg", "stop_not_end"),
                        c("start_not_beg", "tx"),
                        c("start_not_beg:stop_not_end", "stop_not_end"),
                        c("start_not_beg:stop_not_end", "tx"),
                        c("stop_not_end", "tx"))
tumor_metadata_filt <- tumor_metadata_filt %>% dplyr::filter(!(error %in% c("no_starts","no_starts:stop_not_end","no_stops","starts_not_beg:no_stops")))

ggplot(tumor_metadata_filt,aes(x=error,y=log2((tumor_mean_scale_counts+1)/(normal_mean_scale_counts+1))))+
        geom_boxplot()+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
        ylab("log2(mean tumor SB counts/mean normal SB_counts)")+
        xlab("Errors")+
        stat_compare_means(comparisons = my_comparisons,method="wilcox.test")+
  theme(text = element_text(size = 15))
ggsave("/media/theron/My_Passport/SARE_2021/project/spliemutr_images/binders_by_error_3.png",
       width=7, height = 7, units = "in",dpi=300)

print(ggplot(tumor_metadata_filt,aes(x=error,y=log2((tumor_mean_scale+1)/(normal_mean_scale+1))))+
        geom_boxplot()+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
        ylab("log2(mean tumor junc counts/mean normal junc counts) scaled Tumor Junctions")+
        xlab("Errors")+
        stat_compare_means(comparisons = my_comparisons,method="wilcox.test"))

print(ggplot(pathways,aes(x=pathway,y=num_genes,colour=type))+
  geom_point(shape=1)+geom_line()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))

# print(ggplot(pathways_adj,aes(x=pathway,y=log10((num_tumor+1)/(num_normal+1))))+
#   geom_point(shape=1)+labs(y = "log10(tumor genes in pathway/normal genes in pathway)", color = "log10(tumor pathway juncs/normal pathway juncs)")+
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
#
# print(ggplot(pathways_adj,aes(x=pathway,y=log10((tum_juncs+1)/(norm_juncs+1))))+
#   geom_point(shape=1)+labs(y = "log10(tumor pathway juncs/normal pathway juncs)", color = "log10(tumor pathway juncs/normal pathway juncs)")+
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))

print(ggplot(pathways_nonann,aes(x=pathway,y=num_genes,colour=type))+
        geom_point(shape=1)+geom_line()+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))

print(ggplot(pathways_adj_nonann,aes(x=pathway,y=log2((num_tumor+1)/(num_normal+1))))+
        geom_point(shape=1)+labs(y = "log2(tumor genes in pathway/normal genes in pathway) Non ann", color = "log2(tumor pathway juncs/normal pathway juncs)")+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))

print(ggplot(pathways_adj_nonann,aes(x=pathway,y=log2((tum_juncs+1)/(norm_juncs+1))))+
        geom_point(shape=1)+labs(y = "log2(tumor pathway juncs/normal pathway juncs) Non ann", color = "log2(tumor pathway juncs/normal pathway juncs)")+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))

print(pheatmap::pheatmap(top_vals_pathways,cluster_rows = T, cluster_cols = T, show_rownames = F, show_colnames = T,annotation_row = row_vals))
print(ggplot(PIK3_AKT_MTOR,aes(x=deltapsi,y=log2((tumor_mean_scale_counts+1)/(normal_mean_scale_counts+1)),color=verdict))+
  geom_point(shape=1)+
  stat_cor(method = "spearman")+ylab("log2(mean tumor SB counts/mean normal SB counts) scaled"))
print(ggplot(PIK3_AKT_MTOR,aes(x=deltapsi,y=log2((tumor_mean_scale+1)/(normal_mean_scale+1))))+
  geom_point(shape=1)+
  stat_cor(method = "spearman")+geom_smooth(method=lm)+
  ylab("log2(mean tumor junc counts/mean normal junc counts) scaled"))

print(ggplot(hotspots_df,aes(x=seq(nrow(hotspots_df)),y=hotspot))+geom_point()+
        xlab("Protein Position")+ylab("MHC Class 2 Binder Count"))
print(ggplot(hotspots_df_2,aes(x=seq(nrow(hotspots_df_2)),y=hotspot))+geom_point()+
        xlab("Protein Position")+ylab("MHC Class 2 Binder Count"))


# for (i in seq(72)){
#   print(i)
#
#   print(ggplot(SB_psi_dat,aes(y=log2(SB_psi_dat[,i]),x=SB_psi_dat[,72+i],colour=splicemutr_dat$verdict))+
#     geom_point(shape=1)+
#     labs(title=sprintf("%s_%s",colnames(SB_dat_ref_small)[i],keys[colnames(SB_dat_ref_small)[i],"type"]))+
#     xlab("psi")+ylab("SB_count"))
#
#   print(ggplot(tumor_meta_samples,aes(x=deltapsi,y=log2((log2(tumor_meta_samples[,72+i]+1)+1)/(normal_mean_scale_counts+1)),color=verdict))+
#           geom_point(shape=1)+
#           stat_cor(method = "spearman")+ylab("log2(tumor SB counts/mean normal SB counts)")+
#           labs(title=sprintf("%s_%s",colnames(SB_dat_ref_small)[i],keys[colnames(SB_dat_ref_small)[i],"type"])))
#
#   print(ggplot(tumor_meta_samples,aes(x=deltapsi,y=log2((tumor_meta_samples[,i]+1)/(normal_mean_scale+1)),color=verdict))+
#           geom_point(shape=1)+
#           stat_cor(method = "spearman")+
#           ylab("log2(tumor junc counts/mean normal junc counts) scaled")+
#           labs(title=sprintf("%s_%s",colnames(SB_dat_ref_small)[i],keys[colnames(SB_dat_ref_small)[i],"type"])))
#
#
#   print(ggplot(SB_psi_dat_nov,aes(y=log2(SB_psi_dat_nov[,i]),x=SB_psi_dat_nov[,72+i],colour=verdict))+
#           geom_point(shape=1)+
#           labs(title=sprintf("%s_%s",colnames(SB_dat_ref_small)[i],keys[colnames(SB_dat_ref_small)[i],"type"]))+
#           xlab("psi")+ylab("SB_count"))
#
#   print(ggplot(tumor_meta_samples_nov,aes(x=deltapsi,y=log2((tumor_meta_samples_nov[,72+i]+1)/(normal_mean+1)),color=verdict))+
#           geom_point(shape=1)+
#           stat_cor(method = "spearman")+ylab("log2(tumor SB counts/mean normal SB counts) scaled")+
#           labs(title=sprintf("%s_%s",colnames(SB_dat_ref_small)[i],keys[colnames(SB_dat_ref_small)[i],"type"])))
#
#   print(ggplot(tumor_meta_samples_nov,aes(x=deltapsi,y=log2((tumor_meta_samples_nov[,i]+1)/(normal_mean_scale+1)),color=verdict))+
#           geom_point(shape=1)+
#           stat_cor(method = "spearman")+
#           ylab("log2(tumor junc counts/mean normal junc counts) scaled")+
#           labs(title=sprintf("%s_%s",colnames(SB_dat_ref_small)[i],keys[colnames(SB_dat_ref_small)[i],"type"])))
#
# }

dev.off()

# for (i in seq(22,71)){
#
#   tumor_metadata_pathways_filt <- tumor_metadata_pathways[tumor_metadata_pathways[,colnames(tumor_metadata_pathways)[i]] != "diff",]
#   print(ggplot(tumor_metadata_pathways_filt,aes(x=deltapsi,y=log10((tumor_mean_scale_counts+1)/(normal_mean_scale_counts+1)),colour=verdict))+
#           geom_point()+
#           stat_cor(method = "spearman")+ylab("log10(mean tumor SB counts/mean normal SB counts) scaled")+
#       ggtitle(sprintf("%s",colnames(tumor_metadata_pathways)[i])))
#   print(ggplot(tumor_metadata_pathways_filt,aes(x=deltapsi,y=log10((tumor_mean_scale+1)/(normal_mean_scale+1)),colour=verdict))+
#           geom_point()+
#           stat_cor(method = "spearman")+
#           ylab("log10(mean tumor junc counts/mean normal junc counts) scaled")+
#     ggtitle(sprintf("%s",colnames(tumor_metadata_pathways)[i])))
#
#   tumor_metadata_nonann <- tumor_metadata_pathways_filt %>% dplyr::filter(verdict != "annotated")
#   print(ggplot(tumor_metadata_nonann,aes(x=deltapsi,y=log10((tumor_mean_scale_counts+1)/(normal_mean_scale_counts+1)),colour=verdict))+
#           geom_point()+
#           stat_cor(method = "spearman")+ylab("log10(mean tumor SB counts/mean normal SB counts) scaled")+
#           ggtitle(sprintf("%s",colnames(tumor_metadata_pathways)[i])))
#   print(ggplot(tumor_metadata_nonann,aes(x=deltapsi,y=log10((tumor_mean_scale+1)/(normal_mean_scale+1)),colour=verdict))+
#           geom_point()+
#           stat_cor(method = "spearman")+
#           ylab("log10(mean tumor junc counts/mean normal junc counts) scaled")+
#           ggtitle(sprintf("%s",colnames(tumor_metadata_pathways)[i])))
# }
#
# dev.off()



