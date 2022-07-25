#------------------------------------------------------------------------------#
# loading libraries
print("loading libraries")

library(stringr)
library(optparse)
library(dplyr)
library(ensembldb)
library(IRanges)
library(GenomicFeatures)
library(GenomicRanges)
library(ggplot2)
library(umap)
library(pheatmap)

#------------------------------------------------------------------------------#
# plotting things before 03/30/2021

count_file<-"/media/theron/My_Passport/head_and_neck_DARIA/data/leafcutter_10_21_2020/data_perind_numers.counts"
gtf_file<-"/media/theron/My_Passport/reference_genomes/GTF_GFF/GENCODE/gencode.v35.annotation.gtf"
tot_files<-207
pep_dir<-"/media/theron/My_Passport/head_and_neck_DARIA/data/leafcutter_10_21_2020/leafcutter_introns_split"
out_dir<-"/media/theron/My_Passport/head_and_neck_DARIA/data/leafcutter_10_21_2020/leafcutter_introns_split/scores"
ic50<-"/media/theron/My_Passport/head_and_neck_DARIA/data/leafcutter_10_21_2020/leafcutter_introns_split/predictions_1/HLA-B07-02_kmers_9.txt.scores.rds"
groups_file<-"/media/theron/My_Passport/head_and_neck_DARIA/data/leafcutter_10_21_2020/groups_file_noid.txt"
mhc_alleles_file<-"/media/theron/My_Passport/head_and_neck_DARIA/data/leafcutter_10_21_2020/leafcutter_introns_split/mhc_I_to_process.txt"
mhc_alleles_dir<-"/media/theron/My_Passport/head_and_neck_DARIA/data/leafcutter_10_21_2020/leafcutter_introns_split/predictions_1/rows_and_targets"


txdb<-makeTxDbFromGFF(gtf_file) # making the txdb from gtf
cds_by_tx <- cdsBy(txdb,by="tx",use.names=T)
cds_by_tx_df <- data.frame(unlist(cds_by_tx))
cds_by_tx_df$tx_id <- names(unlist(cds_by_tx))

#------------------------------------------------------------------------------#
# loading in the necessary data

print("compiling peps")

for (i in seq(as.numeric(tot_files))){
  file<-sprintf("%s/data_splicemutr%s.txt",pep_dir,i)
  if (file.exists(file)){
    if (i == 1){
      data_splicemutr<-read.table(file)
      data_splicemutr<-unique(data_splicemutr)
    } else {
      peps_fill<-read.table(file)
      data_splicemutr<-rbind(data_splicemutr,peps_fill)
    }
  }
}
rm(peps_fill)

#------------------------------------------------------------------------------#
# calculating the binding count and binders per line

data_splicemutr_filt <- data_splicemutr

# data_splicemutr_filt <- data_splicemutr %>% dplyr::filter(!is.na(peptide))

print("binders")

SB_thresh=50
WB_thresh=500

binders<-lapply(seq(nrow(data_splicemutr_filt)),function(x){
  SB<-scores[[x]][which(scores[[x]] <= SB_thresh)]
  WB<-scores[[x]][which(scores[[x]] <= WB_thresh & scores[[x]] > SB_thresh)]
  binders<-list()
  binders["SB"]<-list(SB)
  binders["WB"]<-list(WB)
  return(binders)
})

binders <- readRDS(sprintf("%s/%s",out_dir,"binders.rds"))

#------------------------------------------------------------------------------#
# determining which rows have binders and which don't

rows_with_binds<-vapply(seq(length(binders)),function(bind_row){
  binds<-binders[[bind_row]]
  if (length(binds$SB)>0 | length(binds$WB)>0){
    return(bind_row)
  } else {
    return(0)
  }
},numeric(1))

data_splicemutr_filt$binders<-NA
data_splicemutr_filt[!is.na(data_splicemutr_filt$peptide),"binders"]<-rows_with_binds

#------------------------------------------------------------------------------#
# placing binders in the splicemutr filt file

# data_splicemutr_filt <- data_splicemutr_filt %>% dplyr::filter(binders>0)
binder_vals<-lapply(data_splicemutr_filt$binders,function(location){
  if(location == 0 | is.na(location)){
    return(NA)
  }
  SB<-paste(unname(binders[[location]]$SB),sep=",",collapse=",")
  WB<-paste(unname(binders[[location]]$WB),sep=",",collapse=",")
  paste(SB,WB,sep="-")
})
data_splicemutr_filt$binder_vals<-NA
data_splicemutr_filt[,"binder_vals"]<-unlist(binder_vals)

#------------------------------------------------------------------------------#
# summarizing the abundance estimates across transcripts tumor

# abund_est_counts_tumor<-abund_est_counts[,tumor_samples]
# abund_est_counts_normal<-abund_est_counts[,normal_samples]
# 
# abund_mean_tumor<-apply(abund_est_counts_tumor,1,mean)
# abund_mean_normal<-apply(abund_est_counts_normal,1,mean)

abund_tpm_tumor<-abund_est_counts[,tumor_samples]
abund_tpm_normal<-abund_est_counts[,normal_samples]

abund_mean_tumor<-apply(abund_tpm_tumor,1,mean)
abund_mean_normal<-apply(abund_tpm_normal,1,mean)

abund_dat<-data.frame(abund_mean_tumor)
abund_dat<-cbind(abund_dat,abund_mean_normal)

colnames(abund_dat)<-c("mean_tumor","mean_normal")
rownames(abund_dat)<-rownames(abund_tpm_tumor)

# filtering for those transcripts from the data
trans<-lapply(data_splicemutr_filt$tx_id,function(tx){
  str_split(tx,"[-]")[[1]]
})
trans<-unique(unlist(trans))

abund_dat_trans<-abund_dat[trans,]

png("/media/theron/My_Passport/head_and_neck_DARIA/data/leafcutter_10_21_2020/leafcutter_introns_split/kallisto_runs/abund_mean_tumor.png",width=700,height=700)
ggplot(abund_dat_trans,aes(x=seq(nrow(abund_dat_trans)),y=mean_tumor))+geom_point()
dev.off()

png("/media/theron/My_Passport/head_and_neck_DARIA/data/leafcutter_10_21_2020/leafcutter_introns_split/kallisto_runs/abund_mean_normal.png",width=700,height=700)
ggplot(abund_dat_trans,aes(x=seq(nrow(abund_dat_trans)),y=mean_normal))+geom_point()
dev.off()

#------------------------------------------------------------------------------#
# filling splicemutr with abundance data

min_abundance_tumor<-vapply(data_splicemutr_filt$tx_id,function(tx){
  min_abund<-min(abund_dat[str_split(tx,"[-]")[[1]],"mean_tumor"])
  if(any(is.na(min_abund))){
    return(-1)
  } else {
    return(min_abund)
  }
},numeric(1))

min_abundance_normal<-vapply(data_splicemutr_filt$tx_id,function(tx){
  min_abund<-min(abund_dat[str_split(tx,"[-]")[[1]],"mean_normal"])
  if(any(is.na(min_abund))){
    return(-1)
  } else {
    return(min_abund)
  }
},numeric(1))

saveRDS(min_abundance_tumor,file="/media/theron/My_Passport/head_and_neck_DARIA/data/leafcutter_10_21_2020/leafcutter_introns_split/kallisto_runs/min_abundance_tumor.rds")
saveRDS(min_abundance_normal,file="/media/theron/My_Passport/head_and_neck_DARIA/data/leafcutter_10_21_2020/leafcutter_introns_split/kallisto_runs/min_abundance_normal.rds")

min_abundance_tumor<-readRDS("/media/theron/My_Passport/head_and_neck_DARIA/data/leafcutter_10_21_2020/leafcutter_introns_split/kallisto_runs/min_abundance_tumor.rds")
min_abundance_normal<-readRDS("/media/theron/My_Passport/head_and_neck_DARIA/data/leafcutter_10_21_2020/leafcutter_introns_split/kallisto_runs/min_abundance_normal.rds")
# add to splicemutr file

rows <- !is.na(data_splicemutr_filt$peptide) & data_splicemutr_filt$binders > 0
data_splicemutr_filt$abund_tumor<-NA
data_splicemutr_filt$abund_normal<-NA
data_splicemutr_filt[rows,"abund_tumor"]<-min_abundance_tumor
data_splicemutr_filt[rows,"abund_normal"]<-min_abundance_normal
# data_splicemutr_filt_abund<-data_splicemutr_filt %>% 
#   dplyr::filter(abund_tumor>=0 & abund_normal >=0)

png("/media/theron/My_Passport/head_and_neck_DARIA/data/leafcutter_10_21_2020/leafcutter_introns_split/kallisto_runs/abund_mean_diff.png",width=700,height=700)
ggplot(data_splicemutr_filt_abund,
       aes(x=seq(nrow(data_splicemutr_filt_abund)),
           y=abund_normal-abund_tumor)) + geom_point()
dev.off()

#------------------------------------------------------------------------------#
# filtering splicemutr data for intragenic trans splicing

intragenic_trans<-vapply(seq(nrow(data_splicemutr_filt)),function(splice_row){
  ann<-data_splicemutr_filt[splice_row,"verdict"]=="annotated"
  double_trans<-str_detect(data_splicemutr_filt[splice_row,"tx_id"],"-")
  if (ann & double_trans){
    return(T)
  } else {
    return(F)
  }
},logical(1))

data_splicemutr_filt$intra_trans<-intragenic_trans

#------------------------------------------------------------------------------#
# processing data_splicemutr counts

counts<-read.table(count_file,header = T,check.names=F)
tumor_counts_mean<-apply(counts[,tumor_samples],1,mean)
tumor_counts_std<-apply(counts[,tumor_samples],1,sd)
normal_counts_mean<-apply(counts[,normal_samples],1,mean)
normal_counts_std<-apply(counts[,normal_samples],1,sd)

counts_dat<-data.frame(tumor_counts_mean)
counts_dat$tumor_sd<-tumor_counts_std
counts_dat$normal_mean<-normal_counts_mean
counts_dat$normal_sd<-normal_counts_std
colnames(counts_dat)<-c("tumor_mean",
                        "tumor_sd",
                        "normal_mean",
                        "normal_sd")
junc_rownames<-rownames(counts_dat)

data_splicemutr_names<-data_splicemutr_filt[,c(1,2,3)]
query_names<-sprintf("%s:%s:%s",
                     data_splicemutr_names[,1],
                     data_splicemutr_names[,2],
                     data_splicemutr_names[,3])
subject_names<-rownames(counts_dat)
subject_names_cut<-unname(vapply(subject_names,function(sname){
  fill<-str_split(sname,"[:]")[[1]]
  return(sprintf("%s:%d:%d",fill[1],as.numeric(fill[2]),as.numeric(fill[3])))
},character(1)))

subject_names_cut_test<-(unname(subject_names_cut) %in% query_names)
counts_dat_cut<-counts_dat[subject_names_cut_test,]
subject_names_full_cut<-subject_names_cut[subject_names_cut_test]
# rownames(counts_dat_cut)<-subject_names_cut[subject_names_cut_test]
counts_tum_nor<-t(vapply(query_names,function(qname){
  loc<-which(subject_names_full_cut == qname)
  if (length(loc)==1){
    return(c(counts_dat_cut[loc,1],
             counts_dat_cut[loc,2],
             counts_dat_cut[loc,3],
             counts_dat_cut[loc,4]))
  } else {
    return(rep(NA,4))
  }
},numeric(4)))
rownames(counts_tum_nor)<-seq(nrow(counts_tum_nor))
colnames(counts_tum_nor)<-c("tumor_mean",
                            "tumor_sd",
                            "normal_mean",
                            "normal_sd")
saveRDS(counts_tum_nor,file="/media/theron/My_Passport/head_and_neck_DARIA/data/leafcutter_10_21_2020/leafcutter_introns_split/kallisto_runs/counts_tum_nor_junc.rds")

counts_tum_nor<-readRDS("/media/theron/My_Passport/head_and_neck_DARIA/data/leafcutter_10_21_2020/leafcutter_introns_split/kallisto_runs/counts_tum_nor_junc.rds")
data_splicemutr_filt_counts<-cbind(data_splicemutr_filt,
                                   data.frame(counts_tum_nor))

#------------------------------------------------------------------------------#
# extracting av binder vals, binder counts SB and WB

binder_vals<-vapply(data_splicemutr_filt$binder_vals,function(vals){
  if (is.na(vals)){
    return(rep(0,4))
  }
  split_vals<-str_split(vals,"[-]")[[1]]
  SB_vals<-split_vals[1]
  SB_vals_split<-str_split(SB_vals,"[,]")
  if (SB_vals_split[[1]][1] == ""){
    SB_out<-c(0,0)
  } else {
    SB_out<-c(mean(as.numeric(SB_vals_split[[1]])),length(SB_vals_split[[1]]))
  }
  WB_vals<-split_vals[2]
  WB_vals_split<-str_split(WB_vals,"[,]")
  if (WB_vals_split[[1]][1] == ""){
    WB_out<-c(0,0)
  } else {
    WB_out<-c(mean(as.numeric(WB_vals_split[[1]])),length(WB_vals_split[[1]]))
  }
  return(c(SB_out,WB_out))
},numeric(4))

binder_vals<-as.data.frame(t(binder_vals))
colnames(binder_vals)<-c("SB_AV","SB_NUM","WB_AV","WB_NUM")

data_splicemutr_filt<-cbind(data_splicemutr_filt,binder_vals)
data_splicemutr_fin<-data_splicemutr_filt %>%
  dplyr::filter(!intragenic_trans)

#------------------------------------------------------------------------------#
# obtaining junction names

splicemutr_names <- sprintf("%s:%s:%s",
                            data_splicemutr_fin$chr,
                            data_splicemutr_fin$start,
                            data_splicemutr_fin$end)
data_splicemutr_fin$junc_names <- splicemutr_names

#------------------------------------------------------------------------------#
# obtaining junction abundance w.r.t. deltapsi

psi_vals <- data_splicemutr_fin$deltapsi
abundance_tumor <- data_splicemutr_fin$abund_tumor
abundance_normal <- data_splicemutr_fin$abund_normal
abundances <- vapply(seq(nrow(data_splicemutr_fin)), function(row_val){
  if (is.na(data_splicemutr_fin[row_val,"binder_vals"])) {
    return(0)
  } else if (psi_vals[row_val] < 0) {
    return(abundance_normal[row_val])
  } else if (psi_vals[row_val] > 0) {
    return(abundance_tumor[row_val])
  } else {
    return(mean(c(abundance_normal[row_val],abundance_tumor[row_val])))
  }
},numeric(1))
data_splicemutr_fin$abundances <- abundances

#------------------------------------------------------------------------------#
# calculating peptide length

pep_len<-vapply(data_splicemutr_fin$peptide,function(pep){
  if (is.na(pep)){
    return(NA)
  } else {
    return(nchar(pep))
  }
},numeric(1))
data_splicemutr_fin$pep_len <- pep_len

#------------------------------------------------------------------------------#
# calculating reference transcript peptide length

tx_ids <- unique(data_splicemutr_fin$tx_id)
# start <- Sys.time()
pep_len_ref <-rep(0,length(tx_ids))
for (i in seq(length(tx_ids))){
  if (i %% 1000 == 0){
    print(i)
  }
  ID_split <- str_split(tx_ids[i],"[-]")[[1]]
  small_cds <- cds_by_tx_df %>% dplyr::filter(tx_id %in% ID_split)
  if (nrow(small_cds) == 0){
    pep_len_ref[i]<-0
  } else {
    pep_len_ref[i]<-sum(small_cds$width)/3
  }
}
pep_len_ref_df <- data.frame(pep_len_ref)
rownames(pep_len_ref_df) <- tx_ids
colnames(pep_len_ref_df) <- pep_len_ref
data_splicemutr_fin$pep_len_ref <- pep_len_ref_df[data_splicemutr_fin$tx_id,"pep_len_ref"]

saveRDS(pep_len_ref_df,file="/media/theron/My_Passport/head_and_neck_DARIA/data/leafcutter_10_21_2020/leafcutter_introns_split/pep_len_ref.rds")

#------------------------------------------------------------------------------#
# calculating junction location in transcript

junc_locs <- vapply(data_splicemutr_fin$tx_junc_loc,function(loc){
  if (is.na(loc)){
    return(NA)
  }
  else {
    locs <- as.numeric(str_split(loc,"[,]")[[1]][1])
  }
},numeric(1))

starts_stops <- vapply(data_splicemutr_fin$start_stop,function(pair){
  if (is.na(pair)){
    return(c(NA,NA))
  }
  else {
    locs <- as.numeric(str_split(pair,"[:]")[[1]])
  }
},numeric(2))
starts_stops <- data.frame(t(starts_stops))
colnames(starts_stops) <- c("beg","stop")


data_splicemutr_fin <- cbind(data_splicemutr_fin,starts_stops)
data_splicemutr_fin$junc_loc_tx <- junc_locs
colnames(data_splicemutr_fin) <- c("chr","start","end","gene","tx_id","modified",
                                   "is_UTR","tx_junc_loc","pep_junc_loc",
                                   "verdict","deltapsi","error","peptide",
                                   "start_stop","orf","binders","binder_vals",
                                   "abund_tumor","abund_normal","intra_trans",
                                   "SB_AV","SB_NUM","WB_AV","WB_NUM",
                                   "junc_names","abundances","pep_len",
                                   "pep_len_ref","beg","stop","junc_loc_tx")

#------------------------------------------------------------------------------#
# calculating the junction-specific scores

unique_names <- unique(data_splicemutr_fin$junc_names)

# # using SB and WB averages for my score
# junction_scores <- vapply(unique_names,function(name){
#   data_splice_small <- data_splicemutr_filt %>% 
#     dplyr::filter(junc_names == name)
# 
#   if (nrow(data_splice_small) == 0){
#     SB_score <- 0
#   } else {
#     data_splice_small[data_splice_small$abundances==-1,"abundances"] = 0
#     abund_SB <- data_splice_small$abundances
#     SB_AV <- -1*(data_splice_small$SB_AV - 50)
#     SB_NUM <- data_splice_small$SB_NUM
#     SB_score <- max(log2(SB_AV)*log2(abund_SB+1))
#   }
#   if (nrow(data_splice_small) == 0) {
#     WB_score <- 0
#   } else {
#     data_splice_small[data_splice_small$abundances==-1,"abundances"] = 0
#     abund_WB <- data_splice_small$abundances
#     WB_AV <- -1*(data_splice_small$WB_AV - 500)
#     WB_NUM <- data_splice_small$WB_NUM
#     WB_score <- max(log2(WB_AV)*log2(abund_WB+1))
#   }
#   return(c(SB_score,WB_score))
# },numeric(2))
# 
# junction_scores <- data.frame(t(junction_scores))
# colnames(junction_scores) <- c("SB_scores","WB_scores")
# 
# row_vals <- rownames(junction_scores)
# deltapsi <- unname(vapply(row_vals,function(row_name){
#   data_splicemutr_small <- data_splicemutr_filt %>% 
#     dplyr::filter(junc_names == row_name)
#   return(data_splicemutr_small$deltapsi[1])
# },numeric(1)))
# verdict <- unname(vapply(row_vals,function(row_name){
#   data_splicemutr_small <- data_splicemutr_filt %>% 
#     dplyr::filter(junc_names == row_name)
#   return(data_splicemutr_small$verdict[1])
# },character(1)))
# pep_len <- unname(vapply(row_vals,function(row_name){
#   data_splicemutr_small <- data_splicemutr_filt %>% 
#     dplyr::filter(junc_names == row_name)
#   return(data_splicemutr_small$pep_len[1])
# },numeric(1)))
# modified <- unname(vapply(row_vals,function(row_name){
#   data_splicemutr_small <- data_splicemutr_filt %>% 
#     dplyr::filter(junc_names == row_name)
#   return(data_splicemutr_small$modified[1])
# },character(1)))
# 
# junction_scores$deltapsi <- deltapsi
# junction_scores$verdict <- verdict
# junction_scores$pep_len <- pep_len
# junction_scores$modified <- modified
# rownames(junction_scores) <- unique_names
# junction_scores_prev <- junction_scores

# using SB and WB counts for my score
junction_counts_SB <- vapply(unique_names,function(name){
  data_splice_small <- data_splicemutr_fin %>% 
    dplyr::filter(junc_names == name)
  
  if (nrow(data_splice_small) == 0){
    SB_score <- 0
    loc <- 0
    orf <- 0
    pep_len <- 0
    pep_len_ref <- 0
    beg <- 0
    stop <- 0
    junc_loc <- 0
    delta_psi <- 0
    verdict <- 0
    modified <- 0
  } else {
    data_splice_small[data_splice_small$abundances==-1,"abundances"] = 0
    abund_SB <- data_splice_small$abundances
    # SB_AV <- -1*(data_splice_small$SB_AV - 50)
    SB_NUM <- data_splice_small$SB_NUM
    SB_scores <- log2(SB_NUM+1)*log2(abund_SB+1)
    SB_score <- max(SB_scores)
    loc <- which(SB_scores == SB_score)
    orf <- data_splice_small[loc[1],"orf"]
    pep_len <- data_splice_small[loc[1],"pep_len"]
    pep_len_ref <- data_splice_small[loc[1],"pep_len_ref"]
    beg <- data_splice_small[loc[1],"beg"]
    stop <- data_splice_small[loc[1],"stop"]
    junc_loc <- data_splice_small[loc[1],"junc_loc_tx"]
    delta_psi <- data_splice_small[loc[1],"deltapsi"]
    verdict <- data_splice_small[loc[1],"verdict"]
    modified <- data_splice_small[loc[1],"modified"]
  }
  return(as.character(c(SB_score,orf,pep_len,pep_len_ref,beg,stop,junc_loc,delta_psi,verdict,modified)))
},character(10))

junction_counts_WB <- vapply(unique_names,function(name){
  data_splice_small <- data_splicemutr_fin %>% 
    dplyr::filter(junc_names == name)
  if (nrow(data_splice_small) == 0) {
    WB_score <- 0
    loc <- 0
    orf <- 0
    pep_len <- 0
    pep_len_ref <- 0
    beg <- 0
    stop <- 0
    junc_loc <- 0
    delta_psi <- 0
    verdict <- 0
    modified <- 0
  } else {
    data_splice_small[data_splice_small$abundances==-1,"abundances"] = 0
    abund_WB <- data_splice_small$abundances
    # WB_AV <- -1*(data_splice_small$WB_AV - 500)
    WB_NUM <- data_splice_small$WB_NUM
    WB_scores <- log2(WB_NUM+1)*log2(abund_WB+1)
    WB_score <- max(WB_scores)
    loc <- which(WB_scores == WB_score)
    orf <- data_splice_small[loc[1],"orf"]
    pep_len <- data_splice_small[loc[1],"pep_len"]
    pep_len_ref <- data_splice_small[loc[1],"pep_len_ref"]
    beg <- data_splice_small[loc[1],"beg"]
    stop <- data_splice_small[loc[1],"stop"]
    junc_loc <- data_splice_small[loc[1],"junc_loc_tx"]
    delta_psi <- data_splice_small[loc[1],"deltapsi"]
    verdict <- data_splice_small[loc[1],"verdict"]
    modified <- data_splice_small[loc[1],"modified"]
  }
  return(as.character(c(WB_score,orf,pep_len,pep_len_ref,beg,stop,junc_loc,delta_psi,verdict,modified)))
},character(10))

junction_counts_SB <- data.frame(t(junction_counts_SB))
colnames(junction_counts_SB) <- c("SB_scores","orf","pep_len","pep_len_ref","beg",
                                  "stop","junc_loc","delta_psi","verdict","modified")

junction_counts_WB <- data.frame(t(junction_counts_WB))
colnames(junction_counts_WB) <- c("WB_scores","orf","pep_len","pep_len_ref","beg",
                               "stop","junc_loc","delta_psi","verdict","modified")

saveRDS(junction_counts_SB,file="/media/theron/My_Passport/head_and_neck_DARIA/data/leafcutter_10_21_2020/leafcutter_introns_split/junction_counts_SB.rds")
saveRDS(junction_counts_WB,file="/media/theron/My_Passport/head_and_neck_DARIA/data/leafcutter_10_21_2020/leafcutter_introns_split/junction_counts_WB.rds")

# # using SB and WB averages for mean score
# junction_scores_mean <- vapply(unique_names,function(name){
#   data_splice_small <- data_splicemutr_filt %>% 
#     dplyr::filter(junc_names == name)
#   
#   if (nrow(data_splice_small) == 0){
#     SB_score <- 0
#   } else {
#     data_splice_small[data_splice_small$abundances==-1,"abundances"] = 0
#     abund_SB <- data_splice_small$abundances
#     SB_AV <- -1*(data_splice_small$SB_AV - 50)
#     SB_NUM <- data_splice_small$SB_NUM
#     SB_score <- mean(log2(SB_AV)*log2(abund_SB+1))
#   }
#   
#   if (nrow(data_splice_small) == 0) {
#     WB_score <- 0
#   } else {
#     data_splice_small[data_splice_small$abundances==-1,"abundances"] = 0
#     abund_WB <- data_splice_small$abundances
#     WB_AV <- -1*(data_splice_small$WB_AV - 500)
#     WB_NUM <- data_splice_small$WB_NUM
#     WB_score <- mean(log2(WB_AV)*log2(abund_WB+1))
#   }
#   return(c(SB_score,WB_score))
# },numeric(2))
# 
# junction_scores_mean <- data.frame(t(junction_scores_mean))
# colnames(junction_scores_mean) <- c("SB_scores","WB_scores")
# 
# row_vals <- rownames(junction_scores_mean)
# deltapsi <- unname(vapply(row_vals,function(row_name){
#   data_splicemutr_small <- data_splicemutr_filt %>% 
#     dplyr::filter(junc_names == row_name)
#   return(data_splicemutr_small$deltapsi[1])
# },numeric(1)))
# verdict <- unname(vapply(row_vals,function(row_name){
#   data_splicemutr_small <- data_splicemutr_filt %>% 
#     dplyr::filter(junc_names == row_name)
#   return(data_splicemutr_small$verdict[1])
# },character(1)))
# pep_len <- unname(vapply(row_vals,function(row_name){
#   data_splicemutr_small <- data_splicemutr_filt %>% 
#     dplyr::filter(junc_names == row_name)
#   return(data_splicemutr_small$pep_len[1])
# },numeric(1)))
# modified <- unname(vapply(row_vals,function(row_name){
#   data_splicemutr_small <- data_splicemutr_filt %>% 
#     dplyr::filter(junc_names == row_name)
#   return(data_splicemutr_small$modified[1])
# },character(1)))
# 
# junction_scores_mean$deltapsi <- deltapsi
# junction_scores_mean$verdict <- verdict
# junction_scores_mean$pep_len <- pep_len
# junction_scores_mean$modified <- modified
# rownames(junction_scores_mean) <- unique_names
# junction_scores_mean_prev <- junction_scores_mean

# # using SB and WB count for mean score 
# junction_counts_mean <- vapply(unique_names,function(name){
#   data_splice_small <- data_splicemutr_filt %>% 
#     dplyr::filter(junc_names == name)
#   
#   if (nrow(data_splice_small) == 0){
#     SB_score <- 0
#   } else {
#     data_splice_small[data_splice_small$abundances==-1,"abundances"] = 0
#     abund_SB <- data_splice_small$abundances
#     SB_AV <- -1*(data_splice_small$SB_AV - 50)
#     SB_NUM <- data_splice_small$SB_NUM
#     SB_score <- mean(log2(SB_NUM+1)*log2(abund_SB+1))
#   }
#   
#   if (nrow(data_splice_small) == 0) {
#     WB_score <- 0
#   } else {
#     data_splice_small[data_splice_small$abundances==-1,"abundances"] = 0
#     abund_WB <- data_splice_small$abundances
#     WB_AV <- -1*(data_splice_small$WB_AV - 500)
#     WB_NUM <- data_splice_small$WB_NUM
#     WB_score <- mean(log2(WB_NUM+1)*log2(abund_WB+1))
#   }
#   return(c(SB_score,WB_score))
# },numeric(2))
# 
# junction_counts_mean <- data.frame(t(junction_counts_mean))
# colnames(junction_counts_mean) <- c("SB_scores","WB_scores")
# 
# row_vals <- rownames(junction_counts_mean)
# deltapsi <- unname(vapply(row_vals,function(row_name){
#   data_splicemutr_small <- data_splicemutr_filt %>% 
#     dplyr::filter(junc_names == row_name)
#   return(data_splicemutr_small$deltapsi[1])
# },numeric(1)))
# verdict <- unname(vapply(row_vals,function(row_name){
#   data_splicemutr_small <- data_splicemutr_filt %>% 
#     dplyr::filter(junc_names == row_name)
#   return(data_splicemutr_small$verdict[1])
# },character(1)))
# pep_len <- unname(vapply(row_vals,function(row_name){
#   data_splicemutr_small <- data_splicemutr_filt %>% 
#     dplyr::filter(junc_names == row_name)
#   return(data_splicemutr_small$pep_len[1])
# },numeric(1)))
# modified <- unname(vapply(row_vals,function(row_name){
#   data_splicemutr_small <- data_splicemutr_filt %>% 
#     dplyr::filter(junc_names == row_name)
#   return(data_splicemutr_small$modified[1])
# },character(1)))
# 
# junction_counts_mean$deltapsi <- deltapsi
# junction_counts_mean$verdict <- verdict
# junction_counts_mean$pep_len <- pep_len
# junction_counts_mean$modified <- modified
# rownames(junction_counts_mean) <- unique_names
# junction_counts_mean_prev <- junction_counts_mean

#------------------------------------------------------------------------------#
# plotting the data ########: 05/02/2021

junction_counts_SB<-readRDS("/media/theron/My_Passport/head_and_neck_DARIA/data/leafcutter_10_21_2020/leafcutter_introns_split/junction_counts_SB.rds")
junction_counts_WB<-readRDS("/media/theron/My_Passport/head_and_neck_DARIA/data/leafcutter_10_21_2020/leafcutter_introns_split/junction_counts_WB.rds")


colnames(junction_counts_SB) <- c("SB_scores","orf","pep_len","pep_len_ref","beg",
                                  "stop","junc_loc","delta_psi","verdict","modified")
row_names <- rownames(junction_counts_SB)
junction_counts_SB$row_names <- row_names
rownames(junction_counts_SB) <- seq(nrow(junction_counts_SB))
junction_counts_SB <- junction_counts_SB %>% dplyr::filter(verdict != "unknown_strand")

colnames(junction_counts_WB) <- c("WB_scores","orf","pep_len","pep_len_ref","beg",
                                  "stop","junc_loc","delta_psi","verdict","modified")
row_names <- rownames(junction_counts_WB)
junction_counts_WB$row_names <- row_names
rownames(junction_counts_WB) <- seq(nrow(junction_counts_WB))
junction_counts_WB <- junction_counts_WB %>% dplyr::filter(verdict != "unknown_strand")

for (i in seq(8)){
  junction_counts_SB[,i] <- as.numeric(junction_counts_SB[,i])
  junction_counts_WB[,i] <- as.numeric(junction_counts_WB[,i])
}

save_dir <- "/media/theron/My_Passport/head_and_neck_DARIA/data/leafcutter_10_21_2020/leafcutter_introns_split/figures/junction_scores"
pdf(sprintf("%s/%s",save_dir,"exploring_scores_05_02_2021.pdf"),width=10,height=10)

print(ggplot(junction_counts_SB,aes(x=delta_psi,y=SB_scores,color=verdict))+
        geom_point()+labs(title="SB dpsi vs score"))
print(ggplot(junction_counts_SB,aes(x=delta_psi,y=log2(pep_len),color=verdict))+
        geom_point()+labs(title="SB dpsi vs pep_len"))
print(ggplot(junction_counts_SB,aes(x=delta_psi,y=log2(pep_len_ref),color=verdict))+
        geom_point()+labs(title="SB dpsi vs pep_len_ref"))
print(ggplot(junction_counts_SB,aes(x=delta_psi,y=log2(pep_len/pep_len_ref),color=verdict))+
        geom_point()+labs(title="SB dpsi vs pep_len/pep_len_ref"))

print(ggplot(junction_counts_SB,aes(x=delta_psi,y=(junc_loc-beg+1)/(stop-beg+1),color=verdict))+
        geom_point()+labs(title="SB dpsi vs Junction Location Relative to start and end of transcript"))
print(ggplot(junction_counts_SB,aes(x=delta_psi,y=(junc_loc-beg+1)/(stop-beg+1),color=verdict))+
        geom_point()+
        labs(title="SB dpsi vs Junction Location Relative to start and end of transcript"))+
        ylim(c(-100,1000))
print(ggplot(junction_counts_SB,aes(y=SB_scores,x=(junc_loc-beg+1)/(stop-beg+1),color=verdict))+
        geom_point()+labs(title="SB Junc location vs score"))
junction_counts_SB_wscore <- junction_counts_SB %>% dplyr::filter(SB_scores > 0)
print(ggplot(junction_counts_SB_wscore,aes(y=SB_scores,x=(junc_loc-beg+1)/(stop-beg+1),color=verdict))+
        geom_point()+labs(title="SB Junc location vs score > 0"))

junction_counts_SB$orf <- as.character(junction_counts_SB$orf)
print(ggplot(junction_counts_SB, aes(x=orf, y=delta_psi, color=orf, group))+ geom_violin())
print(ggplot(junction_counts_SB, aes(x=orf, y=SB_scores, color=orf, group))+ geom_violin())
print(ggplot(junction_counts_SB, aes(x=orf, y=log2(pep_len), color=orf, group))+ geom_violin())
print(ggplot(junction_counts_SB, aes(x=orf, y=log2(pep_len_ref), color=orf, group))+ geom_violin())
print(ggplot(junction_counts_SB, aes(x=orf, y=pep_len/pep_len_ref, color=orf, group))+ geom_violin())
print(ggplot(junction_counts_SB, aes(x=orf, y=pep_len/pep_len_ref, color=orf, group))+ geom_violin())+ylim(c(0,5))
print(ggplot(junction_counts_SB, aes(x=orf, y=(junc_loc-beg+1)/(stop-beg+1), color=orf, group))+ geom_violin())
print(ggplot(junction_counts_SB, aes(x=orf, y=(junc_loc-beg+1)/(stop-beg+1), color=orf, group))+ geom_violin())+ylim(c(-100,1000))

junction_counts_SB_orf3 <- junction_counts_SB %>% dplyr::filter(orf == "3")
print(ggplot(junction_counts_SB_orf3,aes(x=delta_psi,y=SB_scores,color=verdict))+
        geom_point()+labs(title="orf 3 deltapsi vs SB score"))
print(ggplot(junction_counts_SB_orf3,aes(x=delta_psi,y=pep_len/pep_len_ref,color=verdict))+
        geom_point()+labs(title="orf 3 deltapsi vs pep_len ratio"))
junction_counts_SB_orf2 <- junction_counts_SB %>% dplyr::filter(orf == "2")
print(ggplot(junction_counts_SB_orf2,aes(x=delta_psi,y=SB_scores,color=verdict))+
        geom_point()+labs(title="orf 2 deltapsi vs SB score"))
print(ggplot(junction_counts_SB_orf2,aes(x=delta_psi,y=pep_len/pep_len_ref,color=verdict))+
        geom_point()+labs(title="orf 2 deltapsi vs pep_len ratio"))
junction_counts_SB_orf1 <- junction_counts_SB %>% dplyr::filter(orf == "1")
print(ggplot(junction_counts_SB_orf1,aes(x=delta_psi,y=SB_scores,color=verdict))+
        geom_point()+labs(title="orf 1 deltapsi vs SB score"))
print(ggplot(junction_counts_SB_orf1,aes(x=delta_psi,y=pep_len/pep_len_ref,color=verdict))+
        geom_point()+labs(title="orf 1 deltapsi vs pep_len ratio"))

dev.off()

#------------------------------------------------------------------------------#
# plotting the data ########: previous week from ^^^^^

junction_scores <- junction_scores %>% dplyr::filter(verdict != "unknown_strand")
junction_scores$SB_scores_norm <- vapply(seq(nrow(junction_scores)),function(row_val){
  if (junction_scores[row_val,"pep_len"] == 0){return(0)}
  else {return(junction_scores[row_val,"SB_scores"]/junction_scores[row_val,"pep_len"])}
},numeric(1))
junction_scores$WB_scores_norm <- vapply(seq(nrow(junction_scores)),function(row_val){
  if (junction_scores[row_val,"pep_len"] == 0){return(0)}
  else {return(junction_scores[row_val,"WB_scores"]/junction_scores[row_val,"pep_len"])}
},numeric(1))
junction_scores$SB_scores_comp <- vapply(seq(nrow(junction_scores)),function(row_val){
  if (junction_scores[row_val,"SB_scores"] == 0){return(0)}
  else {return(junction_scores[row_val,"SB_scores_norm"]/junction_scores[row_val,"SB_scores"])}
},numeric(1))
junction_scores$WB_scores_comp <- vapply(seq(nrow(junction_scores)),function(row_val){
  if (junction_scores[row_val,"WB_scores"] == 0){return(0)}
  else {return(junction_scores[row_val,"WB_scores_norm"]/junction_scores[row_val,"WB_scores"])}
},numeric(1))

junction_counts <- junction_scores %>% dplyr::filter(verdict != "unknown_strand")
junction_counts$SB_scores_norm <- vapply(seq(nrow(junction_counts)),function(row_val){
  if (junction_counts[row_val,"pep_len"] == 0){return(0)}
  else {return(junction_counts[row_val,"SB_scores"]/junction_counts[row_val,"pep_len"])}
},numeric(1))
junction_counts$WB_scores_norm <- vapply(seq(nrow(junction_counts)),function(row_val){
  if (junction_counts[row_val,"pep_len"] == 0){return(0)}
  else {return(junction_counts[row_val,"WB_scores"]/junction_counts[row_val,"pep_len"])}
},numeric(1))
junction_counts$SB_scores_comp <- vapply(seq(nrow(junction_counts)),function(row_val){
  if (junction_counts[row_val,"SB_scores"] == 0){return(0)}
  else {return(junction_counts[row_val,"SB_scores_norm"]/junction_counts[row_val,"SB_scores"])}
},numeric(1))
junction_counts$WB_scores_comp <- vapply(seq(nrow(junction_counts)),function(row_val){
  if (junction_counts[row_val,"WB_scores"] == 0){return(0)}
  else {return(junction_counts[row_val,"WB_scores_norm"]/junction_counts[row_val,"WB_scores"])}
},numeric(1))

junction_scores_mean <- junction_scores_mean %>% dplyr::filter(verdict != "unknown_strand")
junction_scores_mean$SB_scores_norm <- vapply(seq(nrow(junction_scores_mean)),function(row_val){
  if (junction_scores_mean[row_val,"pep_len"] == 0){return(0)}
  else {return(junction_scores_mean[row_val,"SB_scores"]/junction_scores_mean[row_val,"pep_len"])}
},numeric(1))
junction_scores_mean$WB_scores_norm <- vapply(seq(nrow(junction_scores_mean)),function(row_val){
  if (junction_scores_mean[row_val,"pep_len"] == 0){return(0)}
  else {return(junction_scores_mean[row_val,"WB_scores"]/junction_scores_mean[row_val,"pep_len"])}
},numeric(1))

iter <- 0
junction_scores_mean$SB_scores_comp <- vapply(seq(nrow(junction_scores_mean)),function(row_val){
  iter <<- iter +1
  if (junction_scores_mean[row_val,"SB_scores"] == 0){return(0)}
  else {return(junction_scores_mean[row_val,"SB_scores_norm"]/junction_scores_mean[row_val,"SB_scores"])}
},numeric(1))

junction_scores_mean$WB_scores_comp <- vapply(seq(nrow(junction_scores_mean)),function(row_val){
  if (junction_scores_mean[row_val,"WB_scores"] == 0){return(0)}
  else {return(junction_scores_mean[row_val,"WB_scores_norm"]/junction_scores_mean[row_val,"WB_scores"])}
},numeric(1))

junction_counts_mean <- junction_counts_mean %>% dplyr::filter(verdict != "unknown_strand")
junction_counts_mean$SB_scores_norm <- vapply(seq(nrow(junction_counts_mean)),function(row_val){
  if (junction_counts_mean[row_val,"pep_len"] == 0){return(0)}
  else {return(junction_counts_mean[row_val,"SB_scores"]/junction_counts_mean[row_val,"pep_len"])}
},numeric(1))
junction_counts_mean$WB_scores_norm <- vapply(seq(nrow(junction_counts_mean)),function(row_val){
  if (junction_counts_mean[row_val,"pep_len"] == 0){return(0)}
  else {return(junction_counts_mean[row_val,"WB_scores"]/junction_counts_mean[row_val,"pep_len"])}
},numeric(1))
junction_counts_mean$SB_scores_comp <- vapply(seq(nrow(junction_counts_mean)),function(row_val){
  if (junction_counts_mean[row_val,"SB_scores"] == 0){return(0)}
  else {return(junction_counts_mean[row_val,"SB_scores_norm"]/junction_counts_mean[row_val,"SB_scores"])}
},numeric(1))
junction_counts_mean$WB_scores_comp <- vapply(seq(nrow(junction_counts_mean)),function(row_val){
  if (junction_counts_mean[row_val,"WB_scores"] == 0){return(0)}
  else {return(junction_counts_mean[row_val,"WB_scores_norm"]/junction_counts_mean[row_val,"WB_scores"])}
},numeric(1))

save_dir <- "/media/theron/My_Passport/head_and_neck_DARIA/data/leafcutter_10_21_2020/leafcutter_introns_split/figures/junction_scores"
pdf(sprintf("%s/%s",save_dir,"exploring_scores.pdf"),width=10,height=10)

#

print(ggplot(junction_scores,aes(x=deltapsi,y=SB_scores,color=verdict))+
        geom_point()+labs(title="Junction Scores, Max"))

print(ggplot(junction_scores,aes(x=deltapsi,y=WB_scores,color=verdict))+
        geom_point()+labs(title="Junction Scores, Max"))

print(ggplot(junction_scores,aes(x=deltapsi,y=SB_scores_norm,color=verdict))+
        geom_point()+labs(title="Junction Scores, Max"))

print(ggplot(junction_scores,aes(x=deltapsi,y=WB_scores_norm,color=verdict))+
        geom_point()+labs(title="Junction Scores, Max"))

print(ggplot(junction_scores,aes(x=deltapsi,y=SB_scores_comp,color=verdict))+
        geom_point()+labs(title="Junction Scores, Max"))

print(ggplot(junction_scores,aes(x=deltapsi,y=WB_scores_comp,color=verdict))+
        geom_point()+labs(title="Junction Scores, Max"))

print(ggplot(junction_scores,aes(x=deltapsi,y=log2(pep_len),color=verdict))+
        geom_point()+labs(title="Junction Scores, Max"))

print(ggplot(junction_scores,aes(x=deltapsi,y=SB_scores,color=modified))+
        geom_point()+labs(title="Junction Scores, Max"))

print(ggplot(junction_scores,aes(x=deltapsi,y=WB_scores,color=modified))+
        geom_point()+labs(title="Junction Scores, Max"))

print(ggplot(junction_scores,aes(x=deltapsi,y=SB_scores_norm,color=modified))+
        geom_point()+labs(title="Junction Scores, Max"))

print(ggplot(junction_scores,aes(x=deltapsi,y=WB_scores_norm,color=modified))+
        geom_point()+labs(title="Junction Scores, Max"))

print(ggplot(junction_scores,aes(x=deltapsi,y=SB_scores_comp,color=modified))+
        geom_point()+labs(title="Junction Scores, Max"))

print(ggplot(junction_scores,aes(x=deltapsi,y=WB_scores_comp,color=modified))+
        geom_point()+labs(title="Junction Scores, Max"))

print(ggplot(junction_scores,aes(x=deltapsi,y=log2(pep_len),color=modified))+
        geom_point()+labs(title="Junction Scores, Max"))

print(ggplot(junction_scores,aes(x=modified,y=log2(pep_len),color=modified))+
        geom_violin()+labs(title="Junction Scores, Max"))

print(ggplot(junction_scores,aes(x=modified,y=SB_scores,color=modified))+
        geom_violin()+labs(title="Junction Scores, Max"))

print(ggplot(junction_scores,aes(x=modified,y=SB_scores_norm,color=modified))+
        geom_violin()+labs(title="Junction Scores, Max"))

print(ggplot(junction_scores,aes(x=modified,y=WB_scores,color=modified))+
        geom_violin()+labs(title="Junction Scores, Max"))

print(ggplot(junction_scores,aes(x=modified,y=WB_scores_norm,color=modified))+
        geom_violin()+labs(title="Junction Scores, Max"))

#

junction_scores_non_ann <- junction_scores %>% dplyr::filter(verdict != "annotated") 
print(ggplot(junction_scores_non_ann,aes(x=deltapsi,y=SB_scores,color=verdict))+
        geom_point() + labs(title="Junction Scores, Non-Annotated, Max"))

junction_scores_non_ann <- junction_scores %>% dplyr::filter(verdict != "annotated") 
print(ggplot(junction_scores_non_ann,aes(x=deltapsi,y=WB_scores,color=verdict))+
        geom_point()+labs(title="Junction Scores, Non-Annotated, Max"))

junction_scores_non_ann <- junction_scores %>% dplyr::filter(verdict != "annotated") 
print(ggplot(junction_scores_non_ann,aes(x=deltapsi,y=SB_scores_norm,color=verdict))+
        geom_point() + labs(title="Junction Scores, Non-Annotated, Max"))

junction_scores_non_ann <- junction_scores %>% dplyr::filter(verdict != "annotated") 
print(ggplot(junction_scores_non_ann,aes(x=deltapsi,y=WB_scores_norm,color=verdict))+
        geom_point()+labs(title="Junction Scores, Non-Annotated, Max"))

junction_scores_non_ann <- junction_scores %>% dplyr::filter(verdict != "annotated") 
print(ggplot(junction_scores_non_ann,aes(x=deltapsi,y=SB_scores_comp,color=verdict))+
        geom_point() + labs(title="Junction Scores, Non-Annotated, Max"))

junction_scores_non_ann <- junction_scores %>% dplyr::filter(verdict != "annotated") 
print(ggplot(junction_scores_non_ann,aes(x=deltapsi,y=WB_scores_comp,color=verdict))+
        geom_point()+labs(title="Junction Scores, Non-Annotated, Max"))

junction_scores_non_ann <- junction_scores %>% dplyr::filter(verdict != "annotated") 
print(ggplot(junction_scores_non_ann,aes(x=deltapsi,y=log2(pep_len),color=verdict))+
        geom_point()+labs(title="Junction Scores, Non-Annotated, Max"))

junction_scores_non_ann <- junction_scores %>% dplyr::filter(verdict != "annotated") 
print(ggplot(junction_scores_non_ann,aes(x=deltapsi,y=SB_scores,color=modified))+
        geom_point() + labs(title="Junction Scores, Non-Annotated, Max"))

junction_scores_non_ann <- junction_scores %>% dplyr::filter(verdict != "annotated") 
print(ggplot(junction_scores_non_ann,aes(x=deltapsi,y=WB_scores,color=modified))+
        geom_point()+labs(title="Junction Scores, Non-Annotated, Max"))

junction_scores_non_ann <- junction_scores %>% dplyr::filter(verdict != "annotated") 
print(ggplot(junction_scores_non_ann,aes(x=deltapsi,y=SB_scores_norm,color=modified))+
        geom_point() + labs(title="Junction Scores, Non-Annotated, Max"))

junction_scores_non_ann <- junction_scores %>% dplyr::filter(verdict != "annotated") 
print(ggplot(junction_scores_non_ann,aes(x=deltapsi,y=WB_scores_norm,color=modified))+
        geom_point()+labs(title="Junction Scores, Non-Annotated, Max"))

junction_scores_non_ann <- junction_scores %>% dplyr::filter(verdict != "annotated") 
print(ggplot(junction_scores_non_ann,aes(x=deltapsi,y=SB_scores_comp,color=modified))+
        geom_point() + labs(title="Junction Scores, Non-Annotated, Max"))

junction_scores_non_ann <- junction_scores %>% dplyr::filter(verdict != "annotated") 
print(ggplot(junction_scores_non_ann,aes(x=deltapsi,y=WB_scores_comp,color=modified))+
        geom_point()+labs(title="Junction Scores, Non-Annotated, Max"))

junction_scores_non_ann <- junction_scores %>% dplyr::filter(verdict != "annotated") 
print(ggplot(junction_scores_non_ann,aes(x=deltapsi,y=log2(pep_len),color=modified))+
        geom_point()+labs(title="Junction Scores, Non-Annotated, Max"))

junction_scores_non_ann <- junction_scores %>% dplyr::filter(verdict != "annotated") 
print(ggplot(junction_scores_non_ann,aes(x=modified,y=log2(pep_len),color=modified))+
        geom_violin()+labs(title="Junction Scores, Non-Annotated, Max"))

junction_scores_non_ann <- junction_scores %>% dplyr::filter(verdict != "annotated") 
print(ggplot(junction_scores_non_ann,aes(x=modified,y=SB_scores,color=modified))+
        geom_violin()+labs(title="Junction Scores, Non-Annotated, Max"))
      
junction_scores_non_ann <- junction_scores %>% dplyr::filter(verdict != "annotated") 
print(ggplot(junction_scores_non_ann,aes(x=modified,y=WB_scores,color=modified))+
        geom_violin()+labs(title="Junction Scores, Non-Annotated, Max"))

junction_scores_non_ann <- junction_scores %>% dplyr::filter(verdict != "annotated") 
print(ggplot(junction_scores_non_ann,aes(x=modified,y=SB_scores_norm,color=modified))+
        geom_violin()+labs(title="Junction Scores, Non-Annotated, Max"))

junction_scores_non_ann <- junction_scores %>% dplyr::filter(verdict != "annotated") 
print(ggplot(junction_scores_non_ann,aes(x=modified,y=WB_scores_norm,color=modified))+
        geom_violin()+labs(title="Junction Scores, Non-Annotated, Max"))
#

print(ggplot(junction_counts,aes(x=deltapsi,y=SB_scores,color=verdict))+
        geom_point()+labs(title="Junction Counts, Max"))

print(ggplot(junction_counts,aes(x=deltapsi,y=WB_scores,color=verdict))+
        geom_point()+labs(title="Junction Counts, Max"))

print(ggplot(junction_counts,aes(x=deltapsi,y=SB_scores_norm,color=verdict))+
        geom_point()+labs(title="Junction Counts, Max"))

print(ggplot(junction_counts,aes(x=deltapsi,y=WB_scores_norm,color=verdict))+
        geom_point()+labs(title="Junction Counts, Max"))

print(ggplot(junction_counts,aes(x=deltapsi,y=SB_scores_comp,color=verdict))+
        geom_point()+labs(title="Junction Counts, Max"))

print(ggplot(junction_counts,aes(x=deltapsi,y=WB_scores_comp,color=verdict))+
        geom_point()+labs(title="Junction Counts, Max"))

print(ggplot(junction_counts,aes(x=deltapsi,y=log2(pep_len),color=verdict))+
        geom_point()+labs(title="Junction Counts, Max"))


print(ggplot(junction_counts,aes(x=deltapsi,y=SB_scores,color=modified))+
        geom_point()+labs(title="Junction Counts, Max"))

print(ggplot(junction_counts,aes(x=deltapsi,y=WB_scores,color=modified))+
        geom_point()+labs(title="Junction Counts, Max"))

print(ggplot(junction_counts,aes(x=deltapsi,y=SB_scores_norm,color=modified))+
        geom_point()+labs(title="Junction Counts, Max"))

print(ggplot(junction_counts,aes(x=deltapsi,y=WB_scores_norm,color=modified))+
        geom_point()+labs(title="Junction Counts, Max"))

print(ggplot(junction_counts,aes(x=deltapsi,y=SB_scores_comp,color=modified))+
        geom_point()+labs(title="Junction Counts, Max"))

print(ggplot(junction_counts,aes(x=deltapsi,y=WB_scores_comp,color=modified))+
        geom_point()+labs(title="Junction Counts, Max"))

print(ggplot(junction_counts,aes(x=deltapsi,y=log2(pep_len),color=modified))+
        geom_point()+labs(title="Junction Counts, Max"))

print(ggplot(junction_counts,aes(x=modified,y=SB_scores,color=modified))+
        geom_violin()+labs(title="Junction Counts, Max"))

print(ggplot(junction_counts,aes(x=modified,y=WB_scores_comp,color=modified))+
        geom_violin()+labs(title="Junction Counts, Max"))

print(ggplot(junction_counts,aes(x=modified,y=SB_scores_comp,color=modified))+
        geom_violin()+labs(title="Junction Counts, Max"))

print(ggplot(junction_counts,aes(x=modified,y=WB_scores_comp,color=modified))+
        geom_violin()+labs(title="Junction Counts, Max"))

#

junction_counts_non_ann <- junction_counts %>% dplyr::filter(verdict != "annotated") 
print(ggplot(junction_counts_non_ann,aes(x=deltapsi,y=SB_scores,color=verdict))+
        geom_point() + labs(title="Junction Counts, Non-Annotated, Max"))

junction_counts_non_ann <- junction_counts %>% dplyr::filter(verdict != "annotated") 
print(ggplot(junction_counts_non_ann,aes(x=deltapsi,y=WB_scores,color=verdict))+
        geom_point()+labs(title="Junction Counts, Non-Annotated, Max"))

junction_counts_non_ann <- junction_counts %>% dplyr::filter(verdict != "annotated") 
print(ggplot(junction_counts_non_ann,aes(x=deltapsi,y=SB_scores_norm,color=verdict))+
        geom_point() + labs(title="Junction Counts, Non-Annotated, Max"))

junction_counts_non_ann <- junction_counts %>% dplyr::filter(verdict != "annotated") 
print(ggplot(junction_counts_non_ann,aes(x=deltapsi,y=WB_scores_norm,color=verdict))+
        geom_point()+labs(title="Junction Counts, Non-Annotated, Max"))

junction_counts_non_ann <- junction_counts %>% dplyr::filter(verdict != "annotated") 
print(ggplot(junction_counts_non_ann,aes(x=deltapsi,y=SB_scores_comp,color=verdict))+
        geom_point() + labs(title="Junction Counts, Non-Annotated, Max"))

junction_counts_non_ann <- junction_counts %>% dplyr::filter(verdict != "annotated") 
print(ggplot(junction_counts_non_ann,aes(x=deltapsi,y=WB_scores_comp,color=verdict))+
        geom_point()+labs(title="Junction Counts, Non-Annotated, Max"))

junction_counts_non_ann <- junction_counts %>% dplyr::filter(verdict != "annotated") 
print(ggplot(junction_counts_non_ann,aes(x=deltapsi,y=log2(pep_len),color=verdict))+
        geom_point()+labs(title="Junction Counts, Non-Annotated, Max"))



junction_counts_non_ann <- junction_counts %>% dplyr::filter(verdict != "annotated") 
print(ggplot(junction_counts_non_ann,aes(x=deltapsi,y=SB_scores,color=modified))+
        geom_point() + labs(title="Junction Counts, Non-Annotated, Max"))

junction_counts_non_ann <- junction_counts %>% dplyr::filter(verdict != "annotated") 
print(ggplot(junction_counts_non_ann,aes(x=deltapsi,y=WB_scores,color=modified))+
        geom_point()+labs(title="Junction Counts, Non-Annotated, Max"))

junction_counts_non_ann <- junction_counts %>% dplyr::filter(verdict != "annotated") 
print(ggplot(junction_counts_non_ann,aes(x=deltapsi,y=SB_scores_norm,color=modified))+
        geom_point() + labs(title="Junction Counts, Non-Annotated, Max"))

junction_counts_non_ann <- junction_counts %>% dplyr::filter(verdict != "annotated") 
print(ggplot(junction_counts_non_ann,aes(x=deltapsi,y=WB_scores_norm,color=modified))+
        geom_point()+labs(title="Junction Counts, Non-Annotated, Max"))

junction_counts_non_ann <- junction_counts %>% dplyr::filter(verdict != "annotated") 
print(ggplot(junction_counts_non_ann,aes(x=deltapsi,y=SB_scores_comp,color=modified))+
        geom_point() + labs(title="Junction Counts, Non-Annotated, Max"))

junction_counts_non_ann <- junction_counts %>% dplyr::filter(verdict != "annotated") 
print(ggplot(junction_counts_non_ann,aes(x=deltapsi,y=WB_scores_comp,color=modified))+
        geom_point()+labs(title="Junction Counts, Non-Annotated, Max"))

junction_counts_non_ann <- junction_counts %>% dplyr::filter(verdict != "annotated") 
print(ggplot(junction_counts_non_ann,aes(x=deltapsi,y=log2(pep_len),color=modified))+
        geom_point()+labs(title="Junction Counts, Non-Annotated, Max"))

junction_counts_non_ann <- junction_counts %>% dplyr::filter(verdict != "annotated") 
print(ggplot(junction_counts_non_ann,aes(x=modified,y=SB_scores,color=modified))+
        geom_violin() + labs(title="Junction Counts, Non-Annotated, Max"))

junction_counts_non_ann <- junction_counts %>% dplyr::filter(verdict != "annotated") 
print(ggplot(junction_counts_non_ann,aes(x=modified,y=WB_scores,color=modified))+
        geom_violin()+labs(title="Junction Counts, Non-Annotated, Max"))

junction_counts_non_ann <- junction_counts %>% dplyr::filter(verdict != "annotated") 
print(ggplot(junction_counts_non_ann,aes(x=modified,y=SB_scores_norm,color=modified))+
        geom_violin() + labs(title="Junction Counts, Non-Annotated, Max"))

junction_counts_non_ann <- junction_counts %>% dplyr::filter(verdict != "annotated") 
print(ggplot(junction_counts_non_ann,aes(x=modified,y=WB_scores_norm,color=modified))+
        geom_violin()+labs(title="Junction Counts, Non-Annotated, Max"))

# mean scores and counts
print(ggplot(junction_scores_mean,aes(x=deltapsi,y=SB_scores,color=verdict))+
        geom_point()+labs(title="Junction Scores Mean, Max"))

print(ggplot(junction_scores_mean,aes(x=deltapsi,y=WB_scores,color=verdict))+
        geom_point()+labs(title="Junction Scores Mean, Max"))

print(ggplot(junction_scores_mean,aes(x=deltapsi,y=SB_scores_norm,color=verdict))+
        geom_point()+labs(title="Junction Scores Mean, Max"))

print(ggplot(junction_scores_mean,aes(x=deltapsi,y=WB_scores_norm,color=verdict))+
        geom_point()+labs(title="Junction Scores Mean, Max"))

print(ggplot(junction_scores_mean,aes(x=deltapsi,y=SB_scores_comp,color=verdict))+
        geom_point()+labs(title="Junction Scores Mean, Max"))

print(ggplot(junction_scores_mean,aes(x=deltapsi,y=WB_scores_comp,color=verdict))+
        geom_point()+labs(title="Junction Scores Mean, Max"))

print(ggplot(junction_scores_mean,aes(x=deltapsi,y=log2(pep_len),color=verdict))+
        geom_point()+labs(title="Junction Scores Mean, Max"))


print(ggplot(junction_scores_mean,aes(x=deltapsi,y=SB_scores,color=modified))+
        geom_point()+labs(title="Junction Scores Mean, Max"))

print(ggplot(junction_scores_mean,aes(x=deltapsi,y=WB_scores,color=modified))+
        geom_point()+labs(title="Junction Scores Mean, Max"))

print(ggplot(junction_scores_mean,aes(x=deltapsi,y=SB_scores_norm,color=modified))+
        geom_point()+labs(title="Junction Scores Mean, Max"))

print(ggplot(junction_scores_mean,aes(x=deltapsi,y=WB_scores_norm,color=modified))+
        geom_point()+labs(title="Junction Scores Mean, Max"))

print(ggplot(junction_scores_mean,aes(x=deltapsi,y=SB_scores_comp,color=modified))+
        geom_point()+labs(title="Junction Scores Mean, Max"))

print(ggplot(junction_scores_mean,aes(x=deltapsi,y=WB_scores_comp,color=modified))+
        geom_point()+labs(title="Junction Scores Mean, Max"))

print(ggplot(junction_scores_mean,aes(x=deltapsi,y=log2(pep_len),color=modified))+
        geom_point()+labs(title="Junction Scores Mean, Max"))


print(ggplot(junction_scores_mean,aes(x=modified,y=SB_scores,color=modified))+
        geom_violin()+labs(title="Junction Scores Mean, Max"))

print(ggplot(junction_scores_mean,aes(x=modified,y=WB_scores,color=modified))+
        geom_violin()+labs(title="Junction Scores Mean, Max"))

print(ggplot(junction_scores_mean,aes(x=modified,y=SB_scores_norm,color=modified))+
        geom_violin()+labs(title="Junction Scores Mean, Max"))

print(ggplot(junction_scores_mean,aes(x=modified,y=WB_scores_norm,color=modified))+
        geom_violin()+labs(title="Junction Scores Mean, Max"))

#

junction_scores_non_ann <- junction_scores_mean %>% dplyr::filter(verdict != "annotated") 
print(ggplot(junction_scores_non_ann,aes(x=deltapsi,y=SB_scores,color=verdict))+
  geom_point() + labs(title="Junction Scores Mean, Non-Annotated, Max"))

junction_scores_non_ann <- junction_scores_mean %>% dplyr::filter(verdict != "annotated") 
print(ggplot(junction_scores_non_ann,aes(x=deltapsi,y=WB_scores,color=verdict))+
  geom_point()+labs(title="Junction Scores Mean, Non-Annotated, Max"))

junction_scores_non_ann <- junction_scores_mean %>% dplyr::filter(verdict != "annotated") 
print(ggplot(junction_scores_non_ann,aes(x=deltapsi,y=SB_scores_norm,color=verdict))+
        geom_point() + labs(title="Junction Scores Mean, Non-Annotated, Max"))

junction_scores_non_ann <- junction_scores_mean %>% dplyr::filter(verdict != "annotated") 
print(ggplot(junction_scores_non_ann,aes(x=deltapsi,y=WB_scores_norm,color=verdict))+
        geom_point()+labs(title="Junction Scores Mean, Non-Annotated, Max"))

junction_scores_non_ann <- junction_scores_mean %>% dplyr::filter(verdict != "annotated") 
print(ggplot(junction_scores_non_ann,aes(x=deltapsi,y=SB_scores_comp,color=verdict))+
        geom_point() + labs(title="Junction Scores Mean, Non-Annotated, Max"))

junction_scores_non_ann <- junction_scores_mean %>% dplyr::filter(verdict != "annotated") 
print(ggplot(junction_scores_non_ann,aes(x=deltapsi,y=WB_scores_comp,color=verdict))+
        geom_point()+labs(title="Junction Scores Mean, Non-Annotated, Max"))

junction_scores_non_ann <- junction_scores_mean %>% dplyr::filter(verdict != "annotated") 
print(ggplot(junction_scores_non_ann,aes(x=deltapsi,y=log2(pep_len),color=verdict))+
  geom_point()+labs(title="Junction Scores Mean, Non-Annotated, Max"))


junction_scores_non_ann <- junction_scores_mean %>% dplyr::filter(verdict != "annotated") 
print(ggplot(junction_scores_non_ann,aes(x=deltapsi,y=SB_scores,color=modified))+
        geom_point() + labs(title="Junction Scores Mean, Non-Annotated, Max"))

junction_scores_non_ann <- junction_scores_mean %>% dplyr::filter(verdict != "annotated") 
print(ggplot(junction_scores_non_ann,aes(x=deltapsi,y=WB_scores,color=modified))+
        geom_point()+labs(title="Junction Scores Mean, Non-Annotated, Max"))

junction_scores_non_ann <- junction_scores_mean %>% dplyr::filter(verdict != "annotated") 
print(ggplot(junction_scores_non_ann,aes(x=deltapsi,y=SB_scores_norm,color=modified))+
        geom_point() + labs(title="Junction Scores Mean, Non-Annotated, Max"))

junction_scores_non_ann <- junction_scores_mean %>% dplyr::filter(verdict != "annotated") 
print(ggplot(junction_scores_non_ann,aes(x=deltapsi,y=WB_scores_norm,color=modified))+
        geom_point()+labs(title="Junction Scores Mean, Non-Annotated, Max"))

junction_scores_non_ann <- junction_scores_mean %>% dplyr::filter(verdict != "annotated") 
print(ggplot(junction_scores_non_ann,aes(x=deltapsi,y=SB_scores_comp,color=modified))+
        geom_point() + labs(title="Junction Scores Mean, Non-Annotated, Max"))

junction_scores_non_ann <- junction_scores_mean %>% dplyr::filter(verdict != "annotated") 
print(ggplot(junction_scores_non_ann,aes(x=deltapsi,y=WB_scores_comp,color=verdict))+
        geom_point()+labs(title="Junction Scores Mean, Non-Annotated, Max"))

junction_scores_non_ann <- junction_scores_mean %>% dplyr::filter(verdict != "annotated") 
print(ggplot(junction_scores_non_ann,aes(x=deltapsi,y=log2(pep_len),color=modified))+
        geom_point()+labs(title="Junction Scores Mean, Non-Annotated, Max"))


junction_scores_non_ann <- junction_scores_mean %>% dplyr::filter(verdict != "annotated") 
print(ggplot(junction_scores_non_ann,aes(x=modified,y=SB_scores_norm,color=modified))+
        geom_violin() + labs(title="Junction Scores Mean, Non-Annotated, Max"))

junction_scores_non_ann <- junction_scores_mean %>% dplyr::filter(verdict != "annotated") 
print(ggplot(junction_scores_non_ann,aes(x=modified,y=WB_scores_norm,color=modified))+
        geom_violin()+labs(title="Junction Scores Mean, Non-Annotated, Max"))

junction_scores_non_ann <- junction_scores_mean %>% dplyr::filter(verdict != "annotated") 
print(ggplot(junction_scores_non_ann,aes(x=modified,y=SB_scores,color=modified))+
        geom_violin() + labs(title="Junction Scores Mean, Non-Annotated, Max"))

junction_scores_non_ann <- junction_scores_mean %>% dplyr::filter(verdict != "annotated") 
print(ggplot(junction_scores_non_ann,aes(x=modified,y=WB_scores,color=modified))+
        geom_violin()+labs(title="Junction Scores Mean, Non-Annotated, Max"))

#

print(ggplot(junction_counts_mean,aes(x=deltapsi,y=SB_scores,color=verdict))+
        geom_point()+labs(title="Junction Counts Mean, Max"))

print(ggplot(junction_counts_mean,aes(x=deltapsi,y=WB_scores,color=verdict))+
        geom_point()+labs(title="Junction Counts Mean, Max"))

print(ggplot(junction_counts_mean,aes(x=deltapsi,y=SB_scores_norm,color=verdict))+
        geom_point()+labs(title="Junction Counts Mean, Max"))

print(ggplot(junction_counts_mean,aes(x=deltapsi,y=WB_scores_norm,color=verdict))+
        geom_point()+labs(title="Junction Counts Mean, Max"))

print(ggplot(junction_counts_mean,aes(x=deltapsi,y=SB_scores_comp,color=verdict))+
        geom_point()+labs(title="Junction Counts Mean, Max"))

print(ggplot(junction_counts_mean,aes(x=deltapsi,y=WB_scores_comp,color=verdict))+
        geom_point()+labs(title="Junction Counts Mean, Max"))

print(ggplot(junction_counts_mean,aes(x=deltapsi,y=log2(pep_len),color=verdict))+
        geom_point()+labs(title="Junction Counts Mean, Max"))

print(ggplot(junction_counts_mean,aes(x=deltapsi,y=SB_scores,color=modified))+
        geom_point()+labs(title="Junction Counts Mean, Max"))

print(ggplot(junction_counts_mean,aes(x=deltapsi,y=WB_scores,color=modified))+
        geom_point()+labs(title="Junction Counts Mean, Max"))

print(ggplot(junction_counts_mean,aes(x=deltapsi,y=SB_scores_norm,color=modified))+
        geom_point()+labs(title="Junction Counts Mean, Max"))

print(ggplot(junction_counts_mean,aes(x=deltapsi,y=WB_scores_norm,color=modified))+
        geom_point()+labs(title="Junction Counts Mean, Max"))

print(ggplot(junction_counts_mean,aes(x=deltapsi,y=SB_scores_comp,color=modified))+
        geom_point()+labs(title="Junction Counts Mean, Max"))

print(ggplot(junction_counts_mean,aes(x=deltapsi,y=WB_scores_comp,color=modified))+
        geom_point()+labs(title="Junction Counts Mean, Max"))

print(ggplot(junction_counts_mean,aes(x=deltapsi,y=log2(pep_len),color=modified))+
        geom_point()+labs(title="Junction Counts Mean, Max"))


print(ggplot(junction_counts_mean,aes(x=modified,y=SB_scores,color=modified))+
        geom_violin()+labs(title="Junction Counts Mean, Max"))

print(ggplot(junction_counts_mean,aes(x=modified,y=WB_scores,color=modified))+
        geom_violin()+labs(title="Junction Counts Mean, Max"))

print(ggplot(junction_counts_mean,aes(x=modified,y=SB_scores_norm,color=modified))+
        geom_violin()+labs(title="Junction Counts Mean, Max"))

print(ggplot(junction_counts_mean,aes(x=modified,y=WB_scores_norm,color=modified))+
        geom_violin()+labs(title="Junction Counts Mean, Max"))
#

junction_counts_non_ann <- junction_counts_mean %>% dplyr::filter(verdict != "annotated") 
print(ggplot(junction_counts_non_ann,aes(x=deltapsi,y=SB_scores,color=verdict))+
        geom_point() + labs(title="Junction Counts Mean, Non-Annotated, Max"))

junction_counts_non_ann <- junction_counts_mean %>% dplyr::filter(verdict != "annotated") 
print(ggplot(junction_counts_non_ann,aes(x=deltapsi,y=WB_scores,color=verdict))+
        geom_point()+labs(title="Junction Counts Mean, Non-Annotated, Max"))

junction_counts_non_ann <- junction_counts_mean %>% dplyr::filter(verdict != "annotated") 
print(ggplot(junction_counts_non_ann,aes(x=deltapsi,y=SB_scores_norm,color=verdict))+
        geom_point() + labs(title="Junction Counts Mean, Non-Annotated, Max"))

junction_counts_non_ann <- junction_counts_mean %>% dplyr::filter(verdict != "annotated") 
print(ggplot(junction_counts_non_ann,aes(x=deltapsi,y=WB_scores_norm,color=verdict))+
        geom_point()+labs(title="Junction Counts Mean, Non-Annotated, Max"))

junction_counts_non_ann <- junction_counts_mean %>% dplyr::filter(verdict != "annotated") 
print(ggplot(junction_counts_non_ann,aes(x=deltapsi,y=SB_scores_comp,color=verdict))+
        geom_point() + labs(title="Junction Counts Mean, Non-Annotated, Max"))

junction_counts_non_ann <- junction_counts_mean %>% dplyr::filter(verdict != "annotated") 
print(ggplot(junction_counts_non_ann,aes(x=deltapsi,y=WB_scores_comp,color=verdict))+
        geom_point()+labs(title="Junction Counts Mean, Non-Annotated, Max"))

junction_counts_non_ann <- junction_counts_mean %>% dplyr::filter(verdict != "annotated") 
print(ggplot(junction_counts_non_ann,aes(x=deltapsi,y=log2(pep_len),color=verdict))+
        geom_point()+labs(title="Junction Counts Mean, Non-Annotated, Max"))


junction_counts_non_ann <- junction_counts_mean %>% dplyr::filter(verdict != "annotated") 
print(ggplot(junction_counts_non_ann,aes(x=deltapsi,y=SB_scores,color=modified))+
        geom_point() + labs(title="Junction Counts Mean, Non-Annotated, Max"))

junction_counts_non_ann <- junction_counts_mean %>% dplyr::filter(verdict != "annotated") 
print(ggplot(junction_counts_non_ann,aes(x=deltapsi,y=WB_scores,color=modified))+
        geom_point()+labs(title="Junction Counts Mean, Non-Annotated, Max"))

junction_counts_non_ann <- junction_counts_mean %>% dplyr::filter(verdict != "annotated") 
print(ggplot(junction_counts_non_ann,aes(x=deltapsi,y=SB_scores_norm,color=modified))+
        geom_point() + labs(title="Junction Counts Mean, Non-Annotated, Max"))

junction_counts_non_ann <- junction_counts_mean %>% dplyr::filter(verdict != "annotated") 
print(ggplot(junction_counts_non_ann,aes(x=deltapsi,y=WB_scores_norm,color=modified))+
        geom_point()+labs(title="Junction Counts Mean, Non-Annotated, Max"))

junction_counts_non_ann <- junction_counts_mean %>% dplyr::filter(verdict != "annotated") 
print(ggplot(junction_counts_non_ann,aes(x=deltapsi,y=SB_scores_comp,color=modified))+
        geom_point() + labs(title="Junction Counts Mean, Non-Annotated, Max"))

junction_counts_non_ann <- junction_counts_mean %>% dplyr::filter(verdict != "annotated") 
print(ggplot(junction_counts_non_ann,aes(x=deltapsi,y=WB_scores_comp,color=modified))+
        geom_point()+labs(title="Junction Counts Mean, Non-Annotated, Max"))

junction_counts_non_ann <- junction_counts_mean %>% dplyr::filter(verdict != "annotated") 
print(ggplot(junction_counts_non_ann,aes(x=deltapsi,y=log2(pep_len),color=modified))+
        geom_point()+labs(title="Junction Counts Mean, Non-Annotated, Max"))


junction_counts_non_ann <- junction_counts_mean %>% dplyr::filter(verdict != "annotated") 
print(ggplot(junction_counts_non_ann,aes(x=modified,y=SB_scores,color=modified))+
        geom_violin() + labs(title="Junction Counts Mean, Non-Annotated, Max"))

junction_counts_non_ann <- junction_counts_mean %>% dplyr::filter(verdict != "annotated") 
print(ggplot(junction_counts_non_ann,aes(x=modified,y=WB_scores,color=modified))+
        geom_violin()+labs(title="Junction Counts Mean, Non-Annotated, Max"))

junction_counts_non_ann <- junction_counts_mean %>% dplyr::filter(verdict != "annotated") 
print(ggplot(junction_counts_non_ann,aes(x=modified,y=SB_scores_norm,color=modified))+
        geom_violin() + labs(title="Junction Counts Mean, Non-Annotated, Max"))

junction_counts_non_ann <- junction_counts_mean %>% dplyr::filter(verdict != "annotated") 
print(ggplot(junction_counts_non_ann,aes(x=modified,y=WB_scores_norm,color=modified))+
        geom_violin()+labs(title="Junction Counts Mean, Non-Annotated, Max"))

dev.off()

data_splicemutr_filt_novel <- data_splicemutr_filt %>% dplyr::filter(verdict == "novel annotated pair")
ggplot(data_splicemutr_filt_novel,aes(x=deltapsi,y=abund_tumor-abund_normal,color=verdict))+
  geom_point()

ggplot(data_splicemutr_filt_novel,aes(x=deltapsi,y=log10(pep_len),color=verdict))+
  geom_point()

ggplot(data_splicemutr_filt_novel,aes(x=log2(abund_tumor+1),y=log10(pep_len),color=verdict))+
  geom_point()+xlim(c(0,12))

ggplot(data_splicemutr_filt_novel,aes(x=log2(abund_normal+1),y=log10(pep_len),color=verdict))+
  geom_point()+xlim(c(0,12))

ggplot(data_splicemutr_filt,aes(x=log2(abund_tumor+1),y=log10(pep_len),color=verdict))+
  geom_point()+xlim(c(0,12))

png(sprintf(),width=700,height=700)
ggplot(data_splicemutr_filt,aes(x=log2(abund_normal+1),y=log10(pep_len),color=verdict))+
  geom_point()+xlim(c(0,12))
dev.off()

