#!/usr/bin/env Rscript

library(dplyr)
library(stringr)
library(optparse)
library(TCGAutils)


#------------------------------------------------------------------------------#
# handling command line input

arguments <- parse_args(OptionParser(usage = "",
                                     description="",
                                     option_list=list(
                                       make_option(c("-t","--tcga_cibersort"),
                                                   default = "",
                                                   help="cibersort data"),
                                       make_option(c("-m","--mutation_load"),
                                                   default = "",
                                                   help="the muatation load"),
                                       make_option(c("-l","--leukocyte_fraction"),
                                                   default = "",
                                                   help="the leokocyte fraction"),
                                       make_option(c("-b","--binding_snv"),
                                                   default = "",
                                                   help="the binding snv"),
                                       make_option(c("-c","--tcr_clonality"),
                                                   default = "",
                                                   help="the tcr clonality"),
                                       make_option(c("-p","--tumor_purity"),
                                                   default = "",
                                                   help="the tumor purity "),
                                       make_option(c("-k","--cancer"),
                                                   default = "",
                                                   help="the tcga cancer subtype"),
                                       make_option(c("-g","--genotypes"),
                                                   default = "",
                                                   help="the genotypes file"),
                                       make_option(c("-s","--splicemutr"),
                                                   default = "",
                                                   help="the splicemutr file containing all splicemutr files output from calc_gene_metric"),
                                       make_option(c("-f","--coding_potential"),
                                                   default = "",
                                                   help="the coding potential file"),
                                       make_option(c("-o","--out_dir"),
                                                   default = "",
                                                   help="the output directory"))))
opt=arguments
tcga_cibersort_file <- opt$tcga_cibersort
mutation_load_file <- opt$mutation_load
leukocyte_fraction_file <- opt$leukocyte_fraction
binding_snv_file <- opt$binding_snv
tcr_clonality_file <- opt$tcr_clonality
tumor_purity_file <- opt$tumor_purity
cancer <- opt$cancer
genotypes_file <- opt$genotypes
splicemutr_file <- opt$splicemutr
coding_potential_file <- opt$coding_potential

out_dir <- opt$out_dir

#------------------------------------------------------------------------------#
# Preparing TCGA correlation dataframes


TMB_all_cor <- data.frame(cancer)
rownames(TMB_all_cor) <- cancer
TMB_all_pvals <- data.frame(cancer)
rownames(TMB_all_pvals) <- cancer
TMB_all_cor_norm <- data.frame(cancer)
rownames(TMB_all_cor_norm) <- cancer
TMB_all_pvals_norm <- data.frame(cancer)
rownames(TMB_all_pvals_norm) <- cancer
DA_all_cor <- data.frame(cancer)
rownames(DA_all_cor) <- cancer
DA_all_pvals <- data.frame(cancer)
rownames(DA_all_pvals) <- cancer
cibersort_all_pvals <- data.frame(cancer)
rownames(cibersort_all_pvals)<-cancer
cibersort_all_cor <- data.frame(cancer)
rownames(cibersort_all_cor)<-cancer
cibersort_all_pvals_norm <- data.frame(cancer)
rownames(cibersort_all_pvals_norm)<-cancer
cibersort_all_cor_norm <- data.frame(cancer)
rownames(cibersort_all_cor_norm)<-cancer
dups <- c()

#------------------------------------------------------------------------------#
# Preparing TCGA cibersort data

TCGA_cibersort_all <- readRDS(tcga_cibersort_file)
TCGA_cibersort_all$barcode <- str_replace_all(TCGA_cibersort_all$SampleID,"[.]","-")
TCGA_cibersort_all$sample_ID <- vapply(TCGAbarcode(TCGA_cibersort_all$barcode,sample=T),
                                       function(val){substr(val,1,nchar(val)-1)},character(1))
TCGA_cibersort_all$patient_ID <- TCGAbarcode(TCGA_cibersort_all$barcode)
TCGA_cibersort_all <- TCGA_cibersort_all[as.numeric(TCGA_cibersort_all$P.value)<=0.05,]
cibersort_cells <- c("SampleID","B.cells.naive","B.cells.memory","Plasma.cells",
                     "T.cells.CD8","T.cells.CD4.naive","T.cells.CD4.memory.resting",
                     "T.cells.CD4.memory.activated","T.cells.follicular.helper",
                     "T.cells.regulatory..Tregs.","T.cells.gamma.delta","NK.cells.resting",
                     "NK.cells.activated","Monocytes","Macrophages.M0","Macrophages.M1",
                     "Macrophages.M2","Dendritic.cells.resting","Dendritic.cells.activated",
                     "Mast.cells.resting","Mast.cells.activated","Eosinophils","Neutrophils","P.value","Correlation","RMSE",
                     "barcode","sample_ID","patient_ID")
actual_cibersort_cells <- c("B.cells.naive","B.cells.memory","Plasma.cells",
                            "T.cells.CD8","T.cells.CD4.naive","T.cells.CD4.memory.resting",
                            "T.cells.CD4.memory.activated","T.cells.follicular.helper",
                            "T.cells.regulatory..Tregs.","T.cells.gamma.delta","NK.cells.resting",
                            "NK.cells.activated","Monocytes","Macrophages.M0","Macrophages.M1",
                            "Macrophages.M2","Dendritic.cells.resting","Dendritic.cells.activated",
                            "Mast.cells.resting","Mast.cells.activated","Eosinophils","Neutrophils")

#------------------------------------------------------------------------------#
# Preparing non-silent mutation load


# contains sample id minus the last letter, or rather the vial code

mutation_load <- read.table(mutation_load_file,sep="\t",header=T)

colnames(mutation_load)<-c("cohort","patient_ID","sample_ID","silent_per_mb","non_silent_per_mb")

TMB_vals <- data.frame(median_TMB = vapply(unique(mutation_load$cohort),
                                           function(cohort){median(mutation_load$non_silent_per_mb[mutation_load$cohort==cohort])},numeric(1)))
TMB_vals <- TMB_vals[order(TMB_vals$median_TMB),,drop=F]
mutation_load$cohort <- factor(mutation_load$cohort,levels=rownames(TMB_vals))

# ggplot(mutation_load,aes(x=cohort,y=log10(non_silent_per_mb+1)))+
#   geom_boxplot()+#geom_violin(alpha=0.2,fill="grey")+
#   theme(axis.text.x = element_text(angle = 90))+
#   ylab("log10(Non silent mutations per Mb)")+
#   geom_hline(yintercept=log10(10))+
#   xlab("TCGA Cohort")

#------------------------------------------------------------------------------#
# Preparing Leuokocyte Fraction Data

leukocyte_fraction <- readRDS(leukocyte_fraction_file)
colnames(leukocyte_fraction)<-c("cohort","barcode","fraction")
leukocyte_fraction$patient_ID <- TCGAbarcode(leukocyte_fraction$barcode)
leukocyte_fraction$sample_ID <- vapply(TCGAbarcode(leukocyte_fraction$barcode,sample=T),
                                       function(val){substr(val,1,nchar(val)-1)},character(1))

leuk_vals <- data.frame(median_LF = vapply(unique(leukocyte_fraction$cohort),
                                           function(cohort){median(leukocyte_fraction$fraction[leukocyte_fraction$cohort==cohort])},numeric(1)))
leuk_vals <- leuk_vals[order(leuk_vals$median_LF),,drop=F]
leukocyte_fraction$cohort <- factor(leukocyte_fraction$cohort,levels=rownames(leuk_vals))

# ggplot(leukocyte_fraction,aes(x=cohort,y=fraction))+
#   geom_boxplot()+#geom_violin(alpha=0.2,fill="grey")+
#   theme(axis.text.x = element_text(angle = 90))+
#   xlab("TCGA Cohort")+
#   ylab("Leukocyte Fraction")

#------------------------------------------------------------------------------#
# Preparing Immunogenic SNV Data

binding_snv <- readRDS(binding_snv)
binding_snv$patient_ID <- TCGAbarcode(binding_snv$barcode)
binding_snv$sample_ID <- vapply(TCGAbarcode(binding_snv$barcode,sample=T),
                                function(val){substr(val,1,nchar(val)-1)},character(1))

#------------------------------------------------------------------------------#
# Preparing TCR clonality data

tcr_clonality <-  readRDS(tcr_clonality_file)
tcr_clonality$sample_ID <- vapply(TCGAbarcode(tcr_clonality$AliquotBarcode,sample=T),
                                  function(val){substr(val,1,nchar(val)-1)},character(1))

tcr_clonality_vals <- data.frame(median_tcr = vapply(unique(tcr_clonality$Study),
                                                     function(Study){median(tcr_clonality$numClones[tcr_clonality$Study==Study])},numeric(1)))
tcr_clonality_vals <- tcr_clonality_vals[order(tcr_clonality_vals$median_tcr),,drop=F]
tcr_clonality$Study <- factor(tcr_clonality$Study,levels=rownames(tcr_clonality_vals))

# ggplot(tcr_clonality,aes(x=Study,y=log10(numClones+1)))+
#   geom_boxplot()+#geom_violin(alpha=0.2,fill="grey")+
#   theme(axis.text.x = element_text(angle = 90))+
#   xlab("TCGA Cohort")+
#   ylab("log10(Number of TCR Clones)")

#------------------------------------------------------------------------------#
# Preparing the tumor purity information

tumor_purity_file <-
tumor_purity <- readRDS("./input_data/TCGA_mastercalls.abs_tables_JSedit.fixed.rds")
tumor_purity_filt <- tumor_purity %>% dplyr::filter(solution=="new")
tumor_purity_filt$sample_ID <- vapply(TCGAbarcode(tumor_purity_filt$sample,sample=T),
                                      function(val){substr(val,1,nchar(val)-1)},character(1))

#------------------------------------------------------------------------------#
# processing the tumor genotype file

tumor_geno <- read.table(genotypes_file,header=T) # reading in the tumor genotypes file which contains sample mappings to TCGA barcodes as well as MHC genotype
tumor_geno$sample_ID <- vapply(TCGAbarcode(tumor_geno$aliquot_id,sample=T),
                               function(val){substr(val,1,nchar(val)-1)},character(1)) # generating the sample ID based on the entire TCGA aliquot ID
tumor_geno$patient_ID <- TCGAbarcode(tumor_geno$aliquot_id) # generating the patient ID using the TCGA sample barcode
normal_samples <- tumor_geno$external_id[tumor_geno$type=="N"] # extracting out the normal samples
tumor_samples <- tumor_geno$external_id[tumor_geno$type=="T"] # extracting out the tumor samples

normal_geno <- tumor_geno %>% dplyr::filter(type=="N")
tumor_geno <- tumor_geno %>% dplyr::filter(type=="T") # subsetting the tumor genotypes file by those samples that are tumor samples
tumor_geno <- tumor_geno[complete.cases(tumor_geno),] # only keeping those samples that have complete genotype information

#------------------------------------------------------------------------------#
# combining the tumor splicemutr files batched based off of splice junction expression

splice_mut_files <- read.table(splice_mut_file)
for (j in seq(length(splice_mut_files$V1))){ # splicing antigenicity
  if (j == 1){
    combined_gene_metric_tumor <- readRDS(new_path(splice_mut_files$V1[j],
                                                   sprintf("%s/%s/GENE_METRIC_CP",TCGA_ROOT_DIR,cancer)))
  } else {
    a<-readRDS(new_path(splice_mut_files$V1[j],
                        sprintf("%s/%s/GENE_METRIC_CP",TCGA_ROOT_DIR,cancer)))
    combined_gene_metric_tumor <- cbind(combined_gene_metric_tumor,a)
  }
}

#------------------------------------------------------------------------------#
# loading in the coding potential

coding_potential_LGC_tumor <- readRDS(coding_potential_file) # reading in the LGC coding potential
HIGH_CP_GENES <- rownames(coding_potential_LGC_tumor)[coding_potential_LGC_tumor$coding_potential_LGC>0] # generating the set of genes with high differential coding potential wrt to tumor splicing events

#------------------------------------------------------------------------------#
# combining the normal splicemutr files batched based off of splice junction expression

splice_mut_files <- read.table(splice_mut_file)
for (j in seq(length(splice_mut_files$V1))){
  if (j == 1){
    combined_gene_metric_normal <- readRDS(new_path(splice_mut_files$V1[j],
                                                    sprintf("%s/%s/GENE_METRIC_CP",TCGA_ROOT_DIR,cancer)))
  } else {
    a<-readRDS(new_path(splice_mut_files$V1[j],
                        sprintf("%s/%s/GENE_METRIC_CP",TCGA_ROOT_DIR,cancer)))
    combined_gene_metric_normal <- cbind(combined_gene_metric_normal,a)
  }
}

#------------------------------------------------------------------------------#
# Filtering the tumor and normal splicing antigenicity by represented tumor and normal samples and generating the combined set of tumor and normal genes used for analysis

combined_gene_metric_tumor <- combined_gene_metric_tumor[,tumor_samples[tumor_samples %in% colnames(combined_gene_metric_tumor)]] # filtering out all normal samples from the tumor splicing antignicity
combined_gene_metric_normal <- combined_gene_metric_normal[,normal_samples[normal_samples %in% colnames(combined_gene_metric_normal)]] # filtering out all tumor samples from the normal splicing antigenicity
conglomerate_genes <- unique(rownames(combined_gene_metric_tumor),rownames(combined_gene_metric_normal)) # generating the set of conglomerate genes for the tumor and normal samples to be used for analysis

#------------------------------------------------------------------------------#
# normalizing for the tumor purity

all_external_ids <- c(colnames(combined_gene_metric_tumor),colnames(combined_gene_metric_normal))
all_geno <- rbind(tumor_geno,normal_geno)
external_ids <- all_geno$external_id
tcga_barcodes <- all_geno$aliquot_id
external_dict <- as.list(tcga_barcodes)
names(external_dict)<-external_ids

tcga_barcodes_colnames <- unname(unlist(external_dict[all_external_ids]))
sample_IDs <- vapply(TCGAbarcode(tcga_barcodes_colnames,sample=T),
                     function(val){substr(val,1,nchar(val)-1)},character(1))

tumor_purity_estimate <- vapply(unname(sample_IDs),function(samp){
  a<-which(tumor_purity_filt$sample_ID==samp) # searching for the sample location
  if(length(a)==0){ # if there is no sample ID hit
    return(NA) # return NA
  } else { # if there is a tumor purity hit
    return(tumor_purity_filt[a[1],"purity"]) # return the first sample's tumor purity
  }
},numeric(1))

normal_purity <- tail(tumor_purity_estimate,ncol(combined_gene_metric_normal))
normal_purity[is.na(normal_purity)]<-1
tumor_purity <- head(tumor_purity_estimate,ncol(combined_gene_metric_tumor))

num_tumor_genes <- nrow(combined_gene_metric_tumor)
tumor_purity_matrix <- data.frame(t(matrix(rep(tumor_purity,num_tumor_genes),nrow=ncol(combined_gene_metric_tumor),ncol=num_tumor_genes)))

splicing_antigenicity_tumor_norm <- combined_gene_metric_tumor/tumor_purity_matrix
splicing_antigenicity_normal_norm <- combined_gene_metric_normal

#------------------------------------------------------------------------------#
# Linearizing the normal sample splicing antigenicity per gene for the set of normal samples

normal_all <- splicing_antigenicity_normal_norm[conglomerate_genes,
                                                normal_samples[normal_samples %in% colnames(splicing_antigenicity_normal_norm)]] # filtering by available normal samples and the unique set of conglomerate genes
normal_all[is.na(normal_all)]<-0 # setting all genes that have NA values to zero
normal_all_vec <- c() # creating the linearization vector
normal_filler<-vapply(seq(ncol(normal_all)),function(col_val){ # generating the linearized normal splicing antigenicity
  normal_all_vec <<- c(normal_all_vec,normal_all[,col_val])
  return(T)
},logical(1))
normal_all <- data.frame(SA=normal_all_vec,gene=rep(conglomerate_genes,ncol(normal_all))) # generating dataframe for the linearized normal splicing antigneicity
normal_all$type<-"normal" # marking all as normal gene splicing antigenicity

#------------------------------------------------------------------------------#
# Linearizing the tumor sample splicing antigenicity per gene for the set of tumor samples

tumor_all <- splicing_antigenicity_tumor_norm[conglomerate_genes,tumor_samples[tumor_samples %in% colnames(splicing_antigenicity_tumor_norm)]] # filtering by available tumor samples and the unique set of conglomerate genes
tumor_all[is.na(tumor_all)]<-0 # setting all genes that have NA values to zero
tumor_all_vec <- c() # creating the linearization vector
tumor_filler<-vapply(seq(ncol(tumor_all)),function(col_val){ # generating the linearized tumor splicing antigenicity
  tumor_all_vec <<- c(tumor_all_vec,tumor_all[,col_val])
  return(T)
},logical(1))
tumor_all <- data.frame(SA=tumor_all_vec,gene=rep(conglomerate_genes,ncol(tumor_all))) # generating the dataframe for the linearized tumor splicing antigenicity
tumor_all$type<-"tumor" # marking all as tumor gene splicing antigenicity

#------------------------------------------------------------------------------#
# Combining the tumor and normal per gene and per sample splicing antigenicity and performing a wilcoxon differential test per gene

tumor_normal_all <- rbind(tumor_all,normal_all) # combining the tumor and the normal per gene and sample splicing antigenicity
wilcox_test_val_all <- compare_means(SA ~ type, p.adjust.method = "BH", data=tumor_normal_all,group.by = c("gene")) # performing the wilcoxon differential test per gene based on tumor or normal specification, 1 value per gene per sample
diff_sig_genes <- wilcox_test_val_all$gene[wilcox_test_val_all$p.adj<0.05] # extracting the differential genes using the adjusted p-value
conglomerate_genes <- intersect(conglomerate_genes,diff_sig_genes) # finding the intersection of the conglomerate genes and the genes that are differential

#------------------------------------------------------------------------------#
# Filtering tumor and normal splicing antigenicity by those genes that have differential splicing antigenicity and calculating the mean splicing antigenicity across tumor and normal samples

splicing_antigenicity_tumor_norm <- splicing_antigenicity_tumor_norm[conglomerate_genes,] # filtering the non-linearized tumor splicing antigenicity by conglomerate and differential genes
splicing_antigenicity_normal_norm <- splicing_antigenicity_normal_norm[conglomerate_genes,] # filtering the non-linearized normal splicing antigenicity by conglomerate and differential genes
splicing_antigenicity_normal_norm[is.na(splicing_antigenicity_normal_norm)]<-0

splicing_antigenicity_tumor_norm_mean <- apply(splicing_antigenicity_tumor_norm,1,mean,na.rm=T)
splicing_antigenicity_normal_norm_mean <- apply(splicing_antigenicity_normal_norm,1,mean,na.rm=T)

#------------------------------------------------------------------------------#
# calculating the differential agretopicity between tumor and normal splicing antigenicity

splicing_antigenicity_DA <- data.frame(DA=(splicing_antigenicity_tumor_norm_mean)/(splicing_antigenicity_normal_norm_mean),
                                       DA_inv=(splicing_antigenicity_normal_norm_mean)/(splicing_antigenicity_tumor_norm_mean))

splicing_antigenicity_DA$cancer <- cancer
test_1 <- splicing_antigenicity_DA$DA==Inf & # which have zero normal splicing antigenicity
  !is.na(splicing_antigenicity_DA$DA) # which are not NAs
test_2 <- splicing_antigenicity_DA$DA!=Inf & # which have non-zero normal splicing antigenicity
  !is.na(splicing_antigenicity_DA$DA) & # which are not NAs
  splicing_antigenicity_DA$DA > 1 # which have DA greater than one (because no log transform)
test_3 <- splicing_antigenicity_DA$DA==Inf & # which have zero normal splicing antigenicity
  !is.na(splicing_antigenicity_DA$DA) # which are not NAs

splicing_antigenicity_DA$DA[test_1] <- sample(splicing_antigenicity_DA$DA[test_2],length(which(test_3)),replace=T) # replace those genes that have zero normal splicing antigenicity with sampled genes that have nonzero normal splicing antigenicity and therefore have a differential agretopicity score

test_1 <- splicing_antigenicity_DA$DA_inv==Inf & # which have zero tumor splicing antigenicity
  !is.na(splicing_antigenicity_DA$DA_inv) # which are not NAs
test_2 <- splicing_antigenicity_DA$DA_inv!=Inf & # which have non-zero tumor splicing antigenicity
  !is.na(splicing_antigenicity_DA$DA_inv) & # which are not NAs
  splicing_antigenicity_DA$DA_inv > 1 # which have DA greater than one (becauase no log transform)
test_3 <- splicing_antigenicity_DA$DA_inv==Inf & # which have zero tumor splicing antigenicity
  !is.na(splicing_antigenicity_DA$DA_inv) # which are not NAs

splicing_antigenicity_DA$DA[test_1]<-sample(splicing_antigenicity_DA$DA[test_2],length(which(test_3)),replace=T)# replace those genes that have zero tumor splicing antigenicity with sampled genes that have nonzero tumor splicing antigenicity and therefore have a differential agretopicity score
splicing_antigenicity_DA$genes <- rownames(splicing_antigenicity_DA) # provide gene names for all differential agretopicity scores

#----------------------------------------------------#
# placing differential agretopicity data into container to be saved for later

splicing_antigenicity_DA_cohort <- splicing_antigenicity_DA

#------------------------------------------------------------------------------#
# placing tumor splicing antigenicity data into container to be saved for later

splicing_antigenicity_tumor_genes_cohort <- data.frame(SA=unname(splicing_antigenicity_tumor_norm_mean),
                                                       cancer=rep(cancer,length(splicing_antigenicity_tumor_norm_mean)),
                                                       gene=conglomerate_genes) # taking the average tumor splicing antigenicity per gene and marking with gene names and cancer subtype
splicing_antigenicity_tumor_samples_cohort <- data.frame(SA = apply(splicing_antigenicity_tumor_norm,2,mean,na.rm=T),
                                                         cancer=rep(cancer,ncol(splicing_antigenicity_tumor_norm)),
                                                         type=rep("tumor",ncol(splicing_antigenicity_tumor_norm)),
                                                         sample=colnames(splicing_antigenicity_tumor_norm)) # calculating the average splicing antigenicity per sample across genes and marking with cancer subtype, tumor type, and sample name
splicing_antigenicity_tumor_samples_cohort$sample_ID <- vapply(splicing_antigenicity_tumor_samples_cohort$sample,function(sample){tumor_geno$sample_ID[which(tumor_geno$external_id == sample)[1]]},character(1)) # filling the sample_ID per sample using the tumor genotypes file

#------------------------------------------------------------------------------#
# placing normal splicing antigenicity data into container to be saved for later

splicing_antigenicity_normal_genes_cohort <- data.frame(SA=unname(splicing_antigenicity_normal_norm_mean),
                                                        cancer=rep(cancer,length(splicing_antigenicity_normal_norm_mean)),
                                                        gene=conglomerate_genes)# taking the average normal splicing antigenicity per gene and marking with gene names and cancer subtype
splicing_antigenicity_normal_samples_cohort <- data.frame(SA = apply(splicing_antigenicity_normal_norm,2,mean,na.rm=T),
                                                          cancer=rep(cancer,ncol(splicing_antigenicity_normal_norm)),
                                                          type=rep("normal",ncol(splicing_antigenicity_normal_norm)),
                                                          sample=colnames(splicing_antigenicity_normal_norm)) # calculating the average splicing antigenicity per sample across genes and marking with cancer subtype, tumor type, and sample name
splicing_antigenicity_normal_samples_cohort$sample_ID <- vapply(splicing_antigenicity_normal_samples_cohort$sample,function(sample){normal_geno$sample_ID[which(normal_geno$external_id == sample)[1]]},character(1)) # filling the sample_ID per sample using the normal genotypes file

#------------------------------------------------------------------------------#
# filtering out low DA genes and finding the intersection in differential and high coding potential genes

high_DA_genes <-  splicing_antigenicity_DA[splicing_antigenicity_DA$DA>1,"genes"] # filtering out low DA genes, this was done incorrectly previously, I did not filter out low DA genes it as > 0, this entire section needs to be rerun
high_DA_genes <- high_DA_genes[!is.na(high_DA_genes)] # filtering out all genes with NA DA
HIGH_DA_HIGH_CP_genes <- intersect(intersect(high_DA_genes,HIGH_CP_GENES),diff_sig_genes) # finding the intersection between differential genes, high coding potential genes, and differential genes
splicing_antigenicity_tumor_HIGH_DA_HIGH_CP_diff_sig<-splicing_antigenicity_tumor_norm[HIGH_DA_HIGH_CP_genes,] # saving the high differential agretopicity (DA), high coding potential (CP), and significantly differerential splicing antigenicity genes (between tumor and normal samples)

#------------------------------------------------------------------------------#
# renaming the splicing antigenicity using the tumor sample ID

colnames_splicing_antigenicity <-colnames(splicing_antigenicity_tumor_HIGH_DA_HIGH_CP_diff_sig) # extracting out the external IDs from the splicing antigenicity
colnames(splicing_antigenicity_tumor_HIGH_DA_HIGH_CP_diff_sig) <- vapply(colnames_splicing_antigenicity,function(col_name){
  tumor_geno$sample_ID[which(tumor_geno$external_id == col_name)[1]]
},character(1)) # renaming the external IDs using the sample IDs
splicing_antigenicity_barcode <- vapply(colnames_splicing_antigenicity,function(col_name){
  tumor_geno$aliquot_id[which(tumor_geno$external_id == col_name)[1]]
},character(1)) # extracting the splicing antigenicity barcodes per external ID

#------------------------------------------------------------------------------#
# filtering the metadata using the sample IDs from the colnames of the splicing antigenicity

mutation_load_small <- mutation_load %>% dplyr::filter(sample_ID %in% colnames(splicing_antigenicity_tumor_HIGH_DA_HIGH_CP_diff_sig))
leukocyte_fraction_small <- leukocyte_fraction %>% dplyr::filter(sample_ID %in% colnames(splicing_antigenicity_tumor_HIGH_DA_HIGH_CP_diff_sig))
binding_snv_small <- binding_snv %>% dplyr::filter(sample_ID %in% colnames(splicing_antigenicity_tumor_HIGH_DA_HIGH_CP_diff_sig))
tcr_clonality_small <- tcr_clonality %>% dplyr::filter(sample_ID %in% colnames(splicing_antigenicity_tumor_HIGH_DA_HIGH_CP_diff_sig))
TCGA_cibersort_small <- TCGA_cibersort_all %>% dplyr::filter(sample_ID %in% colnames(splicing_antigenicity_tumor_HIGH_DA_HIGH_CP_diff_sig))

#------------------------------------------------------------------------------#
# preparing the splicing antigenicity for metadata annotations by taking the average per sample

splicing_antigenicity_per_sample<-data.frame(sample_ID=colnames(splicing_antigenicity_tumor_HIGH_DA_HIGH_CP_diff_sig),
                                             barcode=splicing_antigenicity_barcode)
splicing_antigenicity_per_sample$splicing_antigenicity <- vapply(seq(ncol(splicing_antigenicity_tumor_HIGH_DA_HIGH_CP_diff_sig)),
                                                                 function(col_val){mean(as.numeric(splicing_antigenicity_tumor_HIGH_DA_HIGH_CP_diff_sig[,col_val]),
                                                                                        na.rm=T)},numeric(1))

#------------------------------------------------------------------------------#
# loading the splicing antigenicity per sample with sample metadata

# mutation load

splicing_antigenicity_per_sample$non_silent_mutations_per_mb <- unlist(lapply(splicing_antigenicity_per_sample$sample_ID,function(ID){
  a<-which(mutation_load_small$sample_ID==ID) # searching metadata for sample ID
  if (length(a)==1){ # if there are no duplicate samples
    return(mutation_load_small$non_silent_per_mb[a]) # return the metadata for the specific sample
  } else if (length(a)>1){ # if there are duplicates
    dups<<-c(dups,ID) # placing all duplicated samples into a dups object for later
    return(mutation_load_small$non_silent_per_mb[a[1]]) # return the first sample
  } else { # if there are no hits
    return(NA)
  }
}))

# leukocyte_fraction

splicing_antigenicity_per_sample$leukocyte_fraction <- unlist(lapply(seq(length(splicing_antigenicity_per_sample$sample_ID)),function(ID_num){
  ID<-splicing_antigenicity_per_sample$sample_ID[ID_num] # extracting the sample ID
  barcode<-splicing_antigenicity_per_sample$barcode[ID_num] # extracting the barcode
  a<-which(leukocyte_fraction_small$sample_ID==ID) # searching for the sample ID
  if (length(a)==1){ # if there are no duplicates
    return(leukocyte_fraction_small$fraction[a]) # extract metadata for sample
  } else if (length(a)>1){ # if there are duplicates
    aa<-which(leukocyte_fraction$barcode==barcode) # search for the specific barcode
    if (length(aa)==1){ # if there are no barcode duplicates
      return(leukocyte_fraction_small$fraction[aa]) # extract metadata
    } else if (length(aa)>1){ # if there are barcode duplicates
      dups<<-c(dups,ID) # place all duplicated samples into a dups object for filtering later
      return(leukocyte_fraction_small$fraction[aa[1]]) # return the metadata for the first sample
    } else { # if there are no barcode hits
      return(NA)
    }
  } else { # if there are no sample ID hits
    return(NA)
  }
}))

# numberOfImmunogenicMutation

splicing_antigenicity_per_sample$numberOfImmunogenicMutation <- unlist(lapply(seq(length(splicing_antigenicity_per_sample$sample_ID)),function(ID_num){
  ID<-splicing_antigenicity_per_sample$sample_ID[ID_num] # extract out the sample ID
  barcode<-splicing_antigenicity_per_sample$barcode[ID_num] # extract out the barcode
  a<-which(binding_snv_small$sample_ID==ID) # search for the sample ID
  if (length(a)==1){ # if there are no duplicates
    return(binding_snv_small$numberOfImmunogenicMutation[a]) # return the metadata
  } else if (length(a)>1){ # if there are duplicates
    aa<-which(binding_snv_small$barcode==barcode) # search the barcode
    if (length(aa)==1){ # if there are no duplicate barcodes
      return(binding_snv_small$numberOfImmunogenicMutation[aa]) # return the metadata
    } else if (length(aa)>1){ # if there are barcode duplicates
      dups<<-c(dups,ID) # add the duplciated sample ID to the duplicated sample list
      return(binding_snv_small$numberOfImmunogenicMutation[aa[1]]) # return the first sample metadata
    } else { # if there are no barcode hits
      return(NA)
    }
  } else  { # if there are no sample ID hits
    return(NA)
  }
}))

# numClones

splicing_antigenicity_per_sample$numClones <- unlist(lapply(seq(length(splicing_antigenicity_per_sample$sample_ID)),function(ID_num){
  ID<-splicing_antigenicity_per_sample$sample_ID[ID_num] # extracting the sample ID
  barcode<-splicing_antigenicity_per_sample$barcode[ID_num] # extracting the barcode ID
  a<-which(tcr_clonality$sample_ID==ID) # searching for sample ID location
  if (length(a)==1){ # if there are no dupliates
    return(tcr_clonality$numClones[a]) # return the metadata
  } else if (length(a)>1){ # if there are duplicates
    aa<-which(tcr_clonality$AliquotBarcode==barcode) # search the barcode location
    if (length(aa)==1){ # if there are no duplicate barcodes
      return(tcr_clonality$numClones[aa])  # return the metadata
    } else if (length(aa)>1){ # if there are barcode duplicates
      dups<<-c(dups,ID) # mark the sample ID as being duplicated
      return(tcr_clonality$numClones[aa[1]]) # return the first sample metadata
    } else { # if there are no barcode hits
      return(NA)
    }
  } else { # if there are no sample ID hits
    return(NA)
  }
}))


#------------------------------------------------------------------------------#
# marking the average splicing antigenicity per sample with the cancer subtype and labelling the columns

splicing_antigenicity_per_sample$cancer <- cancer
colnames(splicing_antigenicity_per_sample) <- c("sample_ID","barcode","splicing_antigenicity","non_silent_mutations_per_mb",
                                                "leukocyte_fraction","numberOfImmunogenicMutation",
                                                "num_TCR_Clones","cancer")

#------------------------------------------------------------------------------#
# marking the average splicing antigenicity per sample as TMB high vs TMB low

splicing_antigenicity_per_sample <- splicing_antigenicity_per_sample[!(splicing_antigenicity_per_sample$sample_ID %in% dups),]
splicing_antigenicity_per_sample$TMB_TYPE <- "LOW"
splicing_antigenicity_per_sample$TMB_TYPE <- unlist(lapply(splicing_antigenicity_per_sample$non_silent_mutations_per_mb,function(val){
  if (is.na(val)){ # if the tumor purity is NA
    return(NA) # return NA for the tumor purity
  } else if(val>=10){ # if the tumor purity is greater than or equal to 10
    return("TMB HIGH") # return TMB high
  } else { # the the tumor purity if less than 10
    return("TMB LOW") # return TMB low
  }}))

#------------------------------------------------------------------------------#
# marking the average splicing antigenicity per sample with the tumor purity

tumor_purity_estimate <- vapply(splicing_antigenicity_per_sample$sample_ID,function(samp){
  a<-which(tumor_purity_filt$sample_ID==samp) # searching for the sample location
  if(length(a)==0){ # if there is no sample ID hit
    return(NA) # return NA
  } else { # if there is a tumor purity hit
    return(tumor_purity_filt[a[1],"purity"]) # return the first sample's tumor purity
  }
},numeric(1))
splicing_antigenicity_per_sample$purity=tumor_purity_estimate # fill the purity column with the tumor purity metadata

#------------------------------------------------------------------------------#
# generating the CIBERSORTx data per sample

# print("cibersort_per_samp")
cibersort_per_samp <- lapply(seq(length(splicing_antigenicity_per_sample$sample_ID)),
                             function(samp_num){
                               samp<-TCGA_cibersort_small$sample_ID[samp_num] # extracting the sample ID
                               barcode <- TCGA_cibersort_small$barcode[samp_num] # extracting the barcode
                               a<-which(TCGA_cibersort_small$sample_ID==samp) # selecting the sample location(s)
                               if (length(a)==1){ # if there are no duplicates
                                 return(TCGA_cibersort_small[a,cibersort_cells,drop=F]) # return the cibersort data for the sample
                               } else if (length(a)>1){ # if there are ducplicates
                                 aa<-which(TCGA_cibersort_small$barcode==barcode) # search the cibersortx data for the barcode
                                 if (length(aa)==1){ # if there is not barcode duplicate
                                   return(TCGA_cibersort_small[aa,cibersort_cells,drop=F]) # return the cibersort data
                                 } else if (length(aa)>1){ # if there are duplicates
                                   dups<<-c(dups,samp) # mark the samples as duplicates
                                   return(TCGA_cibersort_small[aa[1],cibersort_cells,drop=F]) # return the first sample's cibersort data
                                 } else { # if there are no barcode hits
                                   a<-data.frame(t(rep(NA,length(cibersort_cells)))) # return NAs for the cibersort columns
                                   colnames(a)<-cibersort_cells
                                   return(a)
                                 }
                               } else { # if there are no sample ID hits
                                 a<-data.frame(t(rep(NA,length(cibersort_cells)))) # return NAs for the cibersort columns
                                 colnames(a)<-cibersort_cells
                                 return(a)
                               }
                             })

#------------------------------------------------------------------------------#
# filling the CIBERSORTx dataframe with the CIBERSORTx list information

cibersort_per_samp_df<-cibersort_per_samp[[1]] # fill the cibersortx dataframe with the first sample dataframe
for (k in seq(2,length(cibersort_per_samp))){ # fill the cibersortx dataframe with all consecutive sample dataframes
  cibersort_per_samp_df <- rbind(cibersort_per_samp_df,cibersort_per_samp[[k]]) # rbind combination
}
cibersort_per_samp_df <- cibersort_per_samp_df[complete.cases(cibersort_per_samp_df),] # remove all samples that have NA for cibersort data

#------------------------------------------------------------------------------#
# adding the splicing antigenicity data to the CIBERSORTx dataframe

cibersort_cols_to_add <- c("splicing_antigenicity","purity","non_silent_mutations_per_mb","leukocyte_fraction","numberOfImmunogenicMutation","num_TCR_Clones","cancer")
cibersort_per_samp_df_cols <- colnames(cibersort_per_samp_df)

cibersort_per_samp_df <- cbind(cibersort_per_samp_df,as.data.frame(matrix(unlist(lapply(seq(length(cibersort_per_samp_df$sample_ID)),function(ID_num){
  ID <- cibersort_per_samp_df$sample_ID[ID_num] # extracting out the sample ID
  barcode <- cibersort_per_samp_df$barcode[ID_num] # extracting out the barcode
  a<-which(splicing_antigenicity_per_sample$sample_ID==ID) # searching the splicing antigenicity for the sample ID
  if (length(a)==1){ # if there are no duplicates
    return(as.character(splicing_antigenicity_per_sample[a,cibersort_cols_to_add])) # return the splicing antigenicity column data
  } else if (length(a)>1){ # if there are duplicates
    aa<-which(splicing_antigenicity_per_sample$barcode==barcode) # search the data using the barcode
    if (length(aa)==1){ # if there are no barcode duplicates
      return(splicing_antigenicity_per_sample[aa,cibersort_cols_to_add,drop=F]) # return the splicing antigenicity data
    } else if (length(aa)>1){ # if there are barcode duplicates
      dups<<-c(dups,ID) # store the sample ID as a duplicated sample ID
      return(splicing_antigenicity_per_sample[aa[1],cibersort_cols_to_add,drop=F]) # return the first sample splicing antigenicity data
    } else { # if there are no barcode hits
      return(rep(NA,length(cibersort_cols_to_add))) # return NAs for the splicing antigenicity data
    }
  } else { # if there are no sample ID hits
    return(rep(NA,length(cibersort_cols_to_add))) # return NAs for the splicing antigenicity data
  }
})),byrow=T,nrow=nrow(cibersort_per_samp_df))))
colnames(cibersort_per_samp_df) <- c(cibersort_per_samp_df_cols,cibersort_cols_to_add) # naming the cibersortx data columns so that the splicing antigenicity columns are included
cibersort_per_samp_df <- cibersort_per_samp_df[!(cibersort_per_samp_df$SampleID %in% dups),] # removing all duplicated cibersortx data

#------------------------------------------------------------------------------#
# performing kendall tau statistical tests between the splicing antigenicity and various markers of resoponse to immune checkpoint inhibition

a<-cor.test(splicing_antigenicity_per_sample$splicing_antigenicity,log10(splicing_antigenicity_per_sample$non_silent_mutations_per_mb+1),method="kendall")
TMB_all_cor[cancer,"SA_vs_TMB_cor"]<-a$estimate
TMB_all_pvals[cancer,"SA_vs_TMB_pval"]<-a$p.value

a<-cor.test(splicing_antigenicity_per_sample$splicing_antigenicity,splicing_antigenicity_per_sample$leukocyte_fraction,method="kendall")
TMB_all_cor[cancer,"SA_vs_leuk_frac_cor"]<-a$estimate
TMB_all_pvals[cancer,"SA_vs_leuk_frac_pval"]<-a$p.value


a<-cor.test(splicing_antigenicity_per_sample$splicing_antigenicity,log10(splicing_antigenicity_per_sample$numberOfImmunogenicMutation+1),method="kendall")
TMB_all_cor[cancer,"SA_vs_numImmMut_cor"]<-a$estimate
TMB_all_pvals[cancer,"SA_vs_numImmMut_pval"]<-a$p.value


a<-cor.test(splicing_antigenicity_per_sample$splicing_antigenicity,log2(splicing_antigenicity_per_sample$num_TCR_Clones+1),method="kendall")
TMB_all_cor[cancer,"SA_vs_num_TCR_Clones_cor"]<-a$estimate
TMB_all_pvals[cancer,"SA_vs_num_TCR_Clones_pval"]<-a$p.value


a<-cor.test(log10(splicing_antigenicity_per_sample$numberOfImmunogenicMutation+1),log10(splicing_antigenicity_per_sample$non_silent_mutations_per_mb+1),method="kendall")
TMB_all_cor[cancer,"numImmMut_vs_TMB_cor"]<-a$estimate
TMB_all_pvals[cancer,"numImmMut_vs_TMB_pval"]<-a$p.value


for (cell in actual_cibersort_cells){
  a<-cor.test(as.numeric(cibersort_per_samp_df$splicing_antigenicity),as.numeric(cibersort_per_samp_df[,cell]),method="kendall")
  cibersort_all_cor[cancer,sprintf("%s",cell)]<-as.numeric(a$estimate)
  cibersort_all_pvals[cancer,sprintf("%s",cell)]<-as.numeric(a$p.value)

}

#------------------------------------------------------------------------------#
# saving Rdata for all relevant files

save(splicing_antigenicity_tumor_HIGH_DA_HIGH_CP_diff_sig,
     splicing_antigenicity_per_sample,
     cibersort_per_samp_df,
     TMB_all_cor,
     TMB_all_pvals,
     cibersort_all_cor,
     cibersort_all_pvals,
     file=sprintf("%s/%s_splicing_antigenicity.Rdata",out_dir,cancer))
