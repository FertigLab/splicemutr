# This script is for analyzing splicemutr validation data and Daria's data

#------------------------------------------------------------------------------#
# loading libraries

print("loading libraries")

library(stringi)
library(stringr)
library(optparse)
library(data.table)
library(dplyr)
library(biomaRt)
library(ensembldb)
library(IRanges)
library(Homo.sapiens)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicFeatures)
library(polyester)
library(Biostrings)
library(ggplot2)

#------------------------------------------------------------------------------#
# internal functions

readRDS_global<-function(file){ 
  if(file.exists(file)){
    var<-str_split(basename(file),".rds")[[1]][1]
    assign(var,readRDS(file),envir=globalenv())
  }
}

load_ex_dat<-function(ex_dir){
  in_files<-c("binders.rds","gene_counts.rds","juncs_adjusted.rds","overlapping_junctions.rds","splicemutr_junc_counts.rds")
  if (dir.exists(ex_dir)){
    readRDS_global(sprintf("%s/%s",ex_dir,in_files[1]))
    readRDS_global(sprintf("%s/%s",ex_dir,in_files[2]))
    readRDS_global(sprintf("%s/%s",ex_dir,in_files[3]))
    readRDS_global(sprintf("%s/%s",ex_dir,in_files[4]))
    readRDS_global(sprintf("%s/%s",ex_dir,in_files[5]))
  }
}


#------------------------------------------------------------------------------#
# loading in and processing the validation fasta

validation_fasta<-readDNAStringSet("/media/theron/My_Passport/head_and_neck_DARIA/data/IEDB_benchmarking/HLA-A02-01/validation.fa")

names<-validation_fasta@ranges@NAMES
genes<-unname(vapply(names,function(x){
  str_split(x,"[|]")[[1]][6]
},character(1)))

genes_target<-unique(genes[seq(119)])
genes_random<-unique(genes[seq(119+1,length(genes))])


#------------------------------------------------------------------------------#
# loading in the validation epitopes file

epitope_dir<-"/media/theron/My_Passport/head_and_neck_DARIA/data/IEDB_benchmarking/HLA-A02-01/epitopes.txt"
epitopes<-read.table(epitope_dir,header = T)
epitopes<-unique(epitopes$epitope) # turning dataframe into a character vector
sizes<-vapply(epitopes,function(x){
  nchar(x)
},numeric(1))
epitopes_9<-epitopes[sizes==9]
epitopes_10<-epitopes[sizes==10]
epitopes_11<-epitopes[sizes==11]

#------------------------------------------------------------------------------#
# loading in the validation splicemutr data

validation_splicemutr<-unique(readRDS("/media/theron/My_Passport/head_and_neck_DARIA/data/IEDB_benchmarking/HLA-A02-01/simulated_RNAseq_data/introns/data_splicemutr.rds"))
validation_splicemutr_coding <- validation_splicemutr  %>% dplyr::filter(validation_splicemutr$error %in% c(4,3,2,"tx"))
validation_splicemutr_normal <- validation_splicemutr_coding %>% dplyr::filter(validation_splicemutr_coding$modified == "normal")

#------------------------------------------------------------------------------#
# loading in the validation extracted data, IC9 and the scores

val_ex_dat_dir<-"/media/theron/My_Passport/head_and_neck_DARIA/data/IEDB_benchmarking/HLA-A02-01/simulated_RNAseq_data/extracted_data/IC9"
load_ex_dat(val_ex_dat_dir)
scores_9<-"/media/theron/My_Passport/head_and_neck_DARIA/data/IEDB_benchmarking/HLA-A02-01/simulated_RNAseq_data/kmers/predictions/HLA-A02-01_kmers_9_NetMHC.txt.scores.rds"
scores_9<-readRDS(scores_9)

#------------------------------------------------------------------------------#
# loading in the validation extracted data, IC10

val_ex_dat_dir<-"/media/theron/My_Passport/head_and_neck_DARIA/data/IEDB_benchmarking/HLA-A02-01/simulated_RNAseq_data/extracted_data/IC10"
load_ex_dat(val_ex_dat_dir)
scores_10<-"/media/theron/My_Passport/head_and_neck_DARIA/data/IEDB_benchmarking/HLA-A02-01/simulated_RNAseq_data/kmers/predictions/HLA-A02-01_kmers_10_NetMHC.txt.scores.rds"
scores_10<-readRDS(scores_10)

#------------------------------------------------------------------------------#
# loading in the validation extracted data, IC11

val_ex_dat_dir<-"/media/theron/My_Passport/head_and_neck_DARIA/data/IEDB_benchmarking/HLA-A02-01/simulated_RNAseq_data/extracted_data/IC11"
load_ex_dat(val_ex_dat_dir)
scores_11<-"/media/theron/My_Passport/head_and_neck_DARIA/data/IEDB_benchmarking/HLA-A02-01/simulated_RNAseq_data/kmers/predictions/HLA-A02-01_kmers_11_NetMHC.txt.scores.rds"
scores_11<-readRDS(scores_11)

#------------------------------------------------------------------------------#
# looking to see if the validation targets were identified in the data

targets<-vapply(seq(length(binders)),function(x){ # the epiropes can be found in the targets
  binds<-c(binders[[x]][["SB"]],binders[[x]][["WB"]])
  any(names(binds) %in% epitopes_9)
},logical(1))

# IN
gene_target<-"MLANA"
epitope<-"AAGIGILTV"
rows<-which(validation_splicemutr_coding$gene == gene_target)
target<-vapply(rows,function(x){
  binds<-scores_9[[x]]
  any(epitope %in% names(binds))
},logical(1))
any(target)

# gene not in
gene_target<-"MC1R"
epitope<-"AIIDPLIYA"
rows<-which(validation_splicemutr_coding$gene == gene_target)
target<-vapply(rows,function(x){ # the epiropes can be found in the targets
  binds<-c(binders[[x]][["SB"]],binders[[x]][["WB"]])
  any(epitope %in% names(binds))
},logical(1))
any(target)

# IN
gene_target<-"TEK"
epitope<-"AIRIRTMKM"
rows<-which(validation_splicemutr_coding$gene == gene_target)
target<-vapply(rows,function(x){
  binds<-scores_9[[x]]
  any(epitope %in% names(binds))
},logical(1))
any(target)

# IN
gene_target<-"ACTB"
epitope<-"ALAPSTMKI"
rows<-which(validation_splicemutr_coding$gene == gene_target)
target<-vapply(rows,function(x){ # the epiropes can be found in the targets
  binds<-c(binders[[x]][["SB"]],binders[[x]][["WB"]])
  any(epitope %in% names(binds))
},logical(1))
any(target)

# IN
gene_target<-"ACP3"
epitope<-"ALDVYNGLL"
rows<-which(validation_splicemutr_coding$gene == gene_target)
target<-vapply(rows,function(x){ # the epiropes can be found in the targets
  binds<-c(binders[[x]][["SB"]],binders[[x]][["WB"]])
  any(epitope %in% names(binds))
},logical(1))
any(target)

# IN
gene_target<-"MYH9"
epitope<-"ALEAKIAQL"
rows<-which(validation_splicemutr_coding$gene == gene_target)
target<-vapply(rows,function(x){ # the epiropes can be found in the targets
  binds<-c(binders[[x]][["SB"]],binders[[x]][["WB"]])
  any(epitope %in% names(binds))
},logical(1))
any(target)

# IN
gene_target<-"PLTP"
epitope<-"ALFGALFLA"
rows<-which(validation_splicemutr_coding$gene == gene_target)
target<-vapply(rows,function(x){ # the epiropes can be found in the targets
  binds<-c(binders[[x]][["SB"]],binders[[x]][["WB"]])
  any(epitope %in% names(binds))
},logical(1))
any(target)

# IN
gene_target<-"EIF3D"
epitope<-"ALLAGSEYL"
rows<-which(validation_splicemutr_coding$gene == gene_target)
target<-vapply(rows,function(x){ # the epiropes can be found in the targets
  binds<-c(binders[[x]][["SB"]],binders[[x]][["WB"]])
  any(epitope %in% names(binds))
},logical(1))
any(target)

# IN
gene_target<-"LMNB1"
epitope<-"ALNSKDAAL"
rows<-which(validation_splicemutr_coding$gene == gene_target)
target<-vapply(rows,function(x){ # the epiropes can be found in the targets
  binds<-c(binders[[x]][["SB"]],binders[[x]][["WB"]])
  any(epitope %in% names(binds))
},logical(1))
any(target)

# IN
gene_target<-"ACTB"
epitope<-"ALPHAILRL"
rows<-which(validation_splicemutr_coding$gene == gene_target)
target<-vapply(rows,function(x){ # the epiropes can be found in the targets
  binds<-c(binders[[x]][["SB"]],binders[[x]][["WB"]])
  any(epitope %in% names(binds))
},logical(1))
any(target)

# IN
gene_target<-"DCT"
epitope<-"ALPYWNFAT"
rows<-which(validation_splicemutr_coding$gene == gene_target)
target<-vapply(rows,function(x){ # the epiropes can be found in the targets
  binds<-c(binders[[x]][["SB"]],binders[[x]][["WB"]])
  any(epitope %in% names(binds))
},logical(1))
any(target)

# IN
gene_target<-"EMC7"
epitope<-"ALWGFFPVL"
rows<-which(validation_splicemutr_coding$gene == gene_target)
target<-vapply(rows,function(x){ # the epiropes can be found in the targets
  binds<-c(binders[[x]][["SB"]],binders[[x]][["WB"]])
  any(epitope %in% names(binds))
},logical(1))
any(target)

# IN
gene_target<-"MCF2"
epitope<-"AMLDLLKSV"
rows<-which(validation_splicemutr_coding$gene == gene_target)
target<-vapply(rows,function(x){ # the epiropes can be found in the targets
  binds<-c(binders[[x]][["SB"]],binders[[x]][["WB"]])
  any(epitope %in% names(binds))
},logical(1))
any(target)

# IN
gene_target<-"ACTB"
epitope<-"AMYVAIQAV"
rows<-which(validation_splicemutr_coding$gene == gene_target)
target<-vapply(rows,function(x){ # the epiropes can be found in the targets
  binds<-c(binders[[x]][["SB"]],binders[[x]][["WB"]])
  any(epitope %in% names(binds))
},logical(1))
any(target)

# IN
gene_target<-"CEACAM5"
epitope<-"ATVGIMIGV"
rows<-which(validation_splicemutr_coding$gene == gene_target)
target<-vapply(rows,function(x){ # the epiropes can be found in the targets
  binds<-c(binders[[x]][["SB"]],binders[[x]][["WB"]])
  any(epitope %in% names(binds))
},logical(1))
any(target)

# OUT
gene_target<-"CEACAM5"
epitope<-"AYVHMVTHF"
rows<-which(validation_splicemutr_coding$gene == gene_target)
target<-vapply(rows,function(x){
  binds<-scores_9[[x]]
  any(epitope %in% names(binds))
},logical(1))
any(target)

# IN
gene_target<-"MAGEA2"
epitope<-"FLWGPRALI"
rows<-which(validation_splicemutr_coding$gene == gene_target)
target<-vapply(rows,function(x){ # the epiropes can be found in the targets
  binds<-c(binders[[x]][["SB"]],binders[[x]][["WB"]])
  any(epitope %in% names(binds))
},logical(1))
any(target)

# epitope 10s

# OUT
gene_target<-"MAGEA2"
epitope<-"ALDGGNKHFL"
rows<-which(validation_splicemutr_coding$gene == gene_target)
target<-vapply(rows,function(x){ # the epiropes can be found in the targets
  binds<-c(binders[[x]][["SB"]],binders[[x]][["WB"]])
  any(epitope %in% names(binds))
},logical(1))
any(target)

# OUT
gene_target<-"MYH9"
epitope<-"ALKTELEDTL"
rows<-which(validation_splicemutr_coding$gene == gene_target)
target<-vapply(rows,function(x){ # the epiropes can be found in the targets
  binds<-c(binders[[x]][["SB"]],binders[[x]][["WB"]])
  any(epitope %in% names(binds))
},logical(1))
any(target)

# IN
gene_target<-"PMEL"
epitope<-"AMLGTHTMEV"
rows<-which(validation_splicemutr_coding$gene == gene_target)
target<-vapply(rows,function(x){ # the epiropes can be found in the targets
  binds<-c(binders[[x]][["SB"]],binders[[x]][["WB"]])
  any(epitope %in% names(binds))
},logical(1))
any(target)

# IN
gene_target<-"FASN"
epitope<-"FLFDGSPTYV"
rows<-which(validation_splicemutr_coding$gene == gene_target)
target<-vapply(rows,function(x){ # the epiropes can be found in the targets
  binds<-c(binders[[x]][["SB"]],binders[[x]][["WB"]])
  any(epitope %in% names(binds))
},logical(1))
any(target)

# epitope 11s

# IN
gene_target<-"PLP1"
epitope<-"ALFCGCGHEAL"
rows<-which(validation_splicemutr_coding$gene == gene_target)
target<-vapply(rows,function(x){
  binds<-scores_11[[x]]
  any(epitope %in% names(binds))
},logical(1))
any(target)

#------------------------------------------------------------------------------#
# tx in fasta per gene


genes_from_val<-c("MLANA","MC1R","TEK","ACTB","ACP3","MYH9","PLTP","MAGEA1",'MC1R',"EIF3D","LMNB1",
                  "ACTB","DCT","MAGEA3","EMC7","MCF2","TYR","ACTB","CHCHD2","CEACAM5","CEACAM5",
                  "MAGEA2","MAGEA2","MYH9","PMEL","FASN","H1-4","PLP1")
genes_from_val<-unique(genes_from_val)
genes_from_val_counts<-vapply(genes_from_val,function(x){
  length(which(genes == x))
},numeric(1))
names(genes_from_val_counts)<-genes_from_val

genes_target[!genes_target %in% validation_splicemutr_coding$gene]


#------------------------------------------------------------------------------#
# looking at polyester sims differential splicing

countmat<-readRDS("/media/theron/My_Passport/head_and_neck_DARIA/data/IEDB_benchmarking/HLA-A02-01/simulated_RNAseq_data/countmat.rds")

countmat_ggplot<-data.frame(c(countmat[1:5,1],countmat[1:5,7]))

countmat_ggplot$isoform<-NA
countmat_ggplot$isoform[1:5]<-1:5
countmat_ggplot$isoform[6:10]<-1:5
countmat_ggplot$group<-NA
countmat_ggplot$group[1:5]<-"S1"
countmat_ggplot$group[6:10]<-"S4"

colnames(countmat_ggplot)<-c("reads","isoform","group")

ggplot(data=countmat_ggplot, aes(x=isoform, y=log2(reads), group=group)) +
  geom_line(aes(color=group))+
  geom_point(aes(color=group))+
  xlab("ACP3 Isoforms")+theme(axis.text=element_text(size=12),
                              axis.title=element_text(size=14,face="bold"))
write.csv(log2(countmat),
            file="/media/theron/My_Passport/head_and_neck_DARIA/data/IEDB_benchmarking/HLA-A02-01/simulated_RNAseq_data/countmat.csv")

#------------------------------------------------------------------------------#
# looking at sample-based metrics: number of overlapping junctions per sample

ex_dat_dir<-"/media/theron/My_Passport/head_and_neck_DARIA/data/IEDB_benchmarking/HLA-A02-01/simulated_RNAseq_data/extracted_data/IC9"
load_ex_dat(ex_dat_dir)
overlapping_junctions_dat<-vapply(seq(length(overlapping_junctions)),function(x){
  return(sum(overlapping_junctions[[x]]))
},numeric(1))
over_dat<-data.frame(overlapping_junctions_dat)
over_dat$sample<-seq(12)
colnames(over_dat)<-c("overlaps","samples")
ggplot(data=over_dat, aes(x=samples, y=overlaps, group=1)) +
  geom_line()+
  geom_point()

#------------------------------------------------------------------------------#
# looking at sample-based metrics: number of overlapping junctions per sample as overlaps change

ex_dat_dir<-"/media/theron/My_Passport/head_and_neck_DARIA/data/IEDB_benchmarking/simulated_RNAseq_data_overlaps/ex_dat/overlapping_junctions.rds"
overlapping_junctions<-readRDS(ex_dat_dir)
overlapping_junctions_dat<-vapply(seq(length(overlapping_junctions)),function(x){
  return(sum(overlapping_junctions[[x]]))
},numeric(1))
over_dat<-data.frame(overlapping_junctions_dat)
over_dat$samples<-seq(24)
percentages<-sort(rep(seq(1,8),3))
over_dat$group<-percentages
over_dat_summ<-lapply(seq(length(unique(over_dat$group))),function(x){
  o<-c(mean(log2(over_dat$overlapping_junctions_dat[over_dat$group==x])),sd(log2(over_dat$overlapping_junctions_dat[over_dat$group==x])))
  names(o)<-c("mean","sd")
  return(o)
})
names(over_dat_summ)<-as.character(seq(8))
over_dat_summ<-data.frame(t(data.frame(over_dat_summ)))
over_dat_summ$percentage<-c(10,20,30,40,50,60,70,80)
over_dat_summ$percentage<-as.character(100-over_dat_summ$percentage)
colnames(over_dat)<-c("overlaps","samples","group")
ggplot(data=over_dat_summ, aes(x=percentage, y=mean, group=1)) +
  geom_line()+
  geom_point()+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.1)+
  ylab("log2(overlapping junction count)")+
  xlab("Percentage of Total Isoforms with Expression")



#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#-------------------------------------DARIAS DATA------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# loading in and processing the splicemutr data

splicemutr_file<-"/media/theron/My_Passport/head_and_neck_DARIA/data/leafcutter_10_21_2020/introns/filenames.txt"
data_files<-read.table(splicemutr_file)
for (i in seq(nrow(data_files))){
  file<-data_files[i,]
  if (file.exists(file)){
    if (i == 1){
      peps<-readRDS(file)
      peps<-unique(peps)
    } else {
      peps_fill<-readRDS(file)
      peps_fill<-unique(peps_fill)
      peps<-rbind(peps,peps_fill)
    }
  }
}
peps<-unique(peps)
data_splicemutr<-peps

#------------------------------------------------------------------------------#
# creating metadata plots for all data

# From first frame, translatable based on junction type
codes<-c("0","1","2","3","4","5","6","tx")

peps_small_3<-peps[(peps$error == "0"),]
no_code_seq<-table(peps_small_3$verdict)
no_code_seq<-c(no_code_seq[1],no_code_seq[1],
               no_code_seq[1],no_code_seq[1],no_code_seq[1])
no_code_seq[2:5]<-0

peps_small_4<-peps[(peps$error == "1"),]
stop_before_start<-table(peps_small_4$verdict)

peps_small_5<-peps[(peps$error == "2"),]
multiple_stop<-table(peps_small_5$verdict)

peps_small_6<-peps[(peps$error == "3"),]
stop_not_end<-table(peps_small_6$verdict)
stop_not_end<-c(stop_not_end,stop_not_end[4])
stop_not_end[4]<-0

peps_small_7<-peps[(peps$error == "4"),]
start_not_beginning<-table(peps_small_7$verdict)

peps_small_7<-peps[(peps$error == "5"),]
starts_no_stop<-table(peps_small_7$verdict)

peps_small_7<-peps[(peps$error == "6"),]
no_start_but_stops<-table(peps_small_7$verdict)
no_start_but_stops<-c(no_start_but_stops,no_start_but_stops[4])
no_start_but_stops[4]<-0

peps_small_7<-peps[(peps$error == "tx"),]
translatable<-table(peps_small_7$verdict)

total<-table(peps$verdict)

# creating grouped bar plots 
type_num<-5
condition_num<-8
errors<-c(rep("translatable",type_num), rep("no coding sequence",type_num),
          rep("stop before start",type_num), rep("multiple stop",type_num), 
          rep("sequence aft. stop",type_num), rep("start not at beginning",type_num),
          rep("starts but not stop",type_num),rep("stops but not start",type_num))
vals<-data.frame()
vals<-rbind(vals,translatable)
vals<-rbind(vals,no_code_seq)
vals<-rbind(vals,stop_before_start)
vals<-rbind(vals,multiple_stop)
vals<-rbind(vals,stop_not_end)
vals<-rbind(vals,start_not_beginning)
vals<-rbind(vals,starts_no_stop)
vals<-rbind(vals,no_start_but_stops)

cond<-c("annotated","cryptic fiveprime","cryptic threeprime","cryptic unanchored","novel annotated pair")
colnames(vals)<-cond
tots<-apply(vals,2,sum)

condition<-rep(c(sprintf("%s (%s)","annotated",tots[1]), sprintf("%s (%s)", "cryptic fiveprime", tots[2]),
                 sprintf("%s (%s)","cryptic threeprime",tots[3]), 
                 sprintf("%s (%s)","cryptic unanchored",tots[4]),sprintf("%s (%s)","novel annotated pair",tots[5])),condition_num)
values<-c((translatable/tots)*100, (no_code_seq/tots)*100, (stop_before_start/tots)*100, 
          (multiple_stop/tots)*100, 
          (stop_not_end/tots)*100, 
          (start_not_beginning/tots)*100,
          (starts_no_stop/tots)*100,
          (no_start_but_stops/tots)*100)
vals<-c(translatable, no_code_seq, stop_before_start, multiple_stop,stop_not_end, start_not_beginning,starts_no_stop,no_start_but_stops)
data<-data.frame(errors,condition,values,vals)
# Grouped
ggplot(data, aes(fill=errors, y=values, x=condition)) + 
  geom_bar(position="dodge", stat="identity")+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))+
  ylab("Percentage of Condition Total")

#------------------------------------------------------------------------------------------------------------------------------------------------#
#creating a divergent plot

# psi analysis < 0
codes<-c("0","1","2","3","4","5","6","tx")
peps_small<-peps[peps$deltapsi < 0,]
# From first frame, translatable based on junction type

peps_small_3<-peps_small[(peps_small$error == "0"),]
no_code_seq_2<-table(peps_small_3$verdict)
no_code_seq_2[1:5]<-0

peps_small_4<-peps_small[(peps_small$error == "1"),]
stop_before_start_2<-table(peps_small_4$verdict)

peps_small_5<-peps_small[(peps_small$error == "2"),]
multiple_stop_2<-table(peps_small_5$verdict)

peps_small_6<-peps_small[(peps_small$error == "3"),]
stop_not_end_2<-table(peps_small_6$verdict)
stop_not_end_2<-c(stop_not_end_2,stop_not_end_2[4])
stop_not_end_2[4]<-0

peps_small_7<-peps_small[(peps_small$error == "4"),]
start_not_beginning_2<-table(peps_small_7$verdict)

peps_small_7<-peps_small[(peps_small$error == "5"),]
starts_no_stop_2<-table(peps_small_7$verdict)

peps_small_7<-peps_small[(peps_small$error == "6"),]
no_start_but_stops_2<-table(peps_small_7$verdict)
no_start_but_stops_2<-c(no_start_but_stops_2,no_start_but_stops_2[4])
no_start_but_stops_2[4]<-0

peps_small_7<-peps_small[(peps_small$error == "tx"),]
translatable_2<-table(peps_small_7$verdict)

total_2<-table(peps_small$verdict)


# psi analysis > 0
peps_small<-peps[peps$deltapsi > 0,]
# From first frame, translatable based on junction type

peps_small_3<-peps_small[(peps_small$error == "0"),]
no_code_seq_1<-table(peps_small_3$verdict)
no_code_seq_1<-c(no_code_seq_1[1],no_code_seq_1[1],
                 no_code_seq_1[1],no_code_seq_1[1],no_code_seq_1[1])
no_code_seq_1[2:5]<-0

peps_small_4<-peps_small[(peps_small$error == "1"),]
stop_before_start_1<-table(peps_small_4$verdict)

peps_small_5<-peps_small[(peps_small$error == "2"),]
multiple_stop_1<-table(peps_small_5$verdict)

peps_small_6<-peps_small[(peps_small$error == "3"),]
stop_not_end_1<-table(peps_small_6$verdict)
stop_not_end_1<-c(stop_not_end_1,stop_not_end_1[4])
stop_not_end_1[4]<-0

peps_small_7<-peps_small[(peps_small$error == "4"),]
start_not_beginning_1<-table(peps_small_7$verdict)
start_not_beginning_1<-c(start_not_beginning_1,start_not_beginning_1[4])
start_not_beginning_1[4]<-0

peps_small_7<-peps_small[(peps_small$error == "5"),]
starts_no_stop_1<-table(peps_small_7$verdict)

peps_small_7<-peps_small[(peps_small$error == "6"),]
no_start_but_stops_1<-table(peps_small_7$verdict)
no_start_but_stops_1<-c(no_start_but_stops_1,no_start_but_stops_1[4])
no_start_but_stops_1[4]<-0

peps_small_7<-peps_small[(peps_small$error == "tx"),]
translatable_1<-table(peps_small_7$verdict)

total_1<-table(peps_small$verdict)

# getting differential values
translatable<-translatable_1-translatable_2
no_code_seq<-no_code_seq_1-no_code_seq_2
stop_before_start<-stop_before_start_1-stop_before_start_2
multiple_stop<-multiple_stop_1-multiple_stop_2
stop_not_end<-stop_not_end_1-stop_not_end_2
start_not_beginning<-start_not_beginning_1-start_not_beginning_2
total<-total_1-total_2

# creating grouped bar plots 
type_num<-5
condition_num<-8
errors<-c(rep("translatable",type_num), rep("no coding sequence",type_num),
          rep("stop before start",type_num), rep("multiple stop",type_num), 
          rep("sequence aft. stop",type_num), rep("start not at beginning",type_num),
          rep("starts but not stop",type_num),rep("stops but not start",type_num))
vals<-data.frame()
vals<-rbind(vals,translatable)
vals<-rbind(vals,no_code_seq)
vals<-rbind(vals,stop_before_start)
vals<-rbind(vals,multiple_stop)
vals<-rbind(vals,stop_not_end)
vals<-rbind(vals,start_not_beginning)
vals<-rbind(vals,starts_no_stop)
vals<-rbind(vals,no_start_but_stops)

cond<-c("annotated","cryptic fiveprime","cryptic threeprime","cryptic unanchored","novel annotated pair")
colnames(vals)<-cond
tots<-apply(vals,2,function(x){
  sum(abs(x))
})
condition<-rep(c(sprintf("%s (%s)","annotated",tots[1]), sprintf("%s (%s)", "cryptic fiveprime", tots[2]),
                 sprintf("%s (%s)","cryptic threeprime",tots[3]), 
                 sprintf("%s (%s)","cryptic unanchored",tots[4]),sprintf("%s (%s)","novel annotated pair",tots[5])),condition_num)
values<-c((translatable/tots)*100, (no_code_seq/tots)*100, (stop_before_start/tots)*100, 
          (multiple_stop/tots)*100, 
          (stop_not_end/tots)*100, 
          (start_not_beginning/tots)*100,
          (starts_no_stop/tots)*100,
          (no_start_but_stops/tots)*100)
vals<-c(translatable, no_code_seq, stop_before_start, multiple_stop,stop_not_end, start_not_beginning,starts_no_stop,no_start_but_stops)
data<-data.frame(errors,condition,values,vals)
# Grouped
ggplot(data, aes(fill=errors, y=values, x=condition)) + 
  geom_bar(position="dodge", stat="identity")+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))+
  ylab("Percentage of Condition Total")

#------------------------------------------------------------------------------#
# loading binders

data_splicemutr_filt<- data_splicemutr %>% dplyr::filter(data_splicemutr$error %in% c(4,3,2,"tx"))

binder_dir<-"/media/theron/My_Passport/head_and_neck_DARIA/data/leafcutter_10_21_2020/extracted_data"
if (dir.exists(binder_dir)){
  binder_files<-list.dirs(binder_dir)
  for (x in seq(2,length(binder_files))){
    print(basename(binder_files[x]))
    file<-sprintf("%s/%s",binder_files[x],"binders.rds")
    binder_list<-readRDS(file)
    SB<-vapply(seq(length(binder_list)),function(y){
      SB<-binder_list[[y]][["SB"]]
      return(length(SB))
    },numeric(1))
    WB<-vapply(seq(length(binder_list)),function(y){
      WB<-binder_list[[y]][["WB"]]
      return(length(WB))
    },numeric(1))
    allele_kmer_name<-basename(binder_files[x])
    data_splicemutr_filt[,sprintf("%s_%s",allele_kmer_name,"SB")]<-SB
    data_splicemutr_filt[,sprintf("%s_%s",allele_kmer_name,"WB")]<-WB
  }
}
data_colnames<-sort(colnames(data_splicemutr_filt)[17:40])[11:24]

#------------------------------------------------------------------------------#
# looking at binding with top scoring differentially expressed genes
# binders start at col 17

load("/media/theron/My_Passport/head_and_neck_DARIA/data/leafcutter_10_21_2020/data.Rdata")
cluster_order<-order(clusters$FDR,decreasing = F)
top_clusters<-clusters[cluster_order[seq(100)],]
clus<-top_clusters$clusterID
data_splicemutr_top_clus<-data_splicemutr_filt %>% dplyr::filter(data_splicemutr_filt$cluster %in% clus)
data_splicemutr_top_clus<-unique(data_splicemutr_top_clus[,c(seq(1,4),seq(6,15),seq(26,ncol(data_splicemutr_top_clus)))])
genes<-unique((data_splicemutr_top_clus$gene))
binding_summary<-lapply(seq(length(genes)),function(x){
  attach(data_splicemutr_top_clus)
  curr_dat<-data_splicemutr_top_clus %>% dplyr::filter(gene==genes[x])
  curr_dat<-data.frame(curr_dat[,"peptide"],curr_dat[,seq(15,ncol(curr_dat))])
  apply(curr_dat[,seq(2,ncol(curr_dat))],2,mean)
})
names(binding_summary)<-gene
binding_summary<-data.frame(binding_summary)
binding_summary[is.na(binding_summary)]<-0

pheatmap::pheatmap(binding_summary,display_numbers=T,cluster_rows=F,cluster_cols=F)

binding_summary_tumor<-lapply(seq(length(gene)),function(x){
  attach(data_splicemutr_top_clus)
  curr_dat<-data_splicemutr_top_clus %>% dplyr::filter(gene==genes[x] & deltapsi>0)
  curr_dat<-data.frame(curr_dat[,"peptide"],curr_dat[,seq(15,ncol(curr_dat))])
  apply(curr_dat[,seq(2,ncol(curr_dat))],2,mean)
})
names(binding_summary_tumor)<-gene
binding_summary_tumor<-data.frame(binding_summary_tumor)
binding_summary_tumor[is.na(binding_summary_tumor)]<-0


binding_summary_normal<-lapply(seq(length(gene)),function(x){
  attach(data_splicemutr_top_clus)
  curr_dat<-data_splicemutr_top_clus %>% dplyr::filter(gene==genes[x] & deltapsi<0)
  curr_dat<-data.frame(curr_dat[,"peptide"],curr_dat[,seq(15,ncol(curr_dat))])
  apply(curr_dat[,seq(2,ncol(curr_dat))],2,mean)
})
names(binding_summary_normal)<-gene
binding_summary_normal<-data.frame(binding_summary_normal)
binding_summary_normal[is.na(binding_summary_normal)]<-0

diff_binding_summary<-binding_summary_tumor-binding_summary_normal

pheatmap::pheatmap(diff_binding_summary,display_numbers=T,cluster_rows=F,cluster_cols=F)

#------------------------------------------------------------------------------#
# plotting SB and WB counts per junction, error type, and HLA type:

peps<-data_splicemutr_filt
codes<-c("2","3","4","tx")
codes_named<-c("multiple_stop","stop_not_end","start_not_beginning","translatable")
cond<-c("annotated","cryptic_fiveprime","cryptic_threeprime","cryptic_unanchored","novel annotated pair")

for (i in seq(1,14,by=2)){
  
  HLA<-data_colnames[i]
  
  junc_counts<-lapply(seq(length(codes)),function(x){
    peps_small<-peps %>% dplyr::filter(peps$error == codes[x])
    cond_vals<-vapply(seq(length(cond)),function(y){
      peps_smaller<-peps_small %>% dplyr::filter(peps_small$verdict == cond[y])
      nrow(peps_smaller)
    },numeric(1))
    names(cond_vals)<-cond
    return(cond_vals)
  })
  names(junc_counts)<-codes_named
  junc_counts<-data.frame(junc_counts)
  
  
  SB_vals<-lapply(seq(length(codes)),function(x){
    peps_small<-peps %>% dplyr::filter(peps$error == codes[x])
    cond_vals<-vapply(seq(length(cond)),function(y){
      peps_smaller<-peps_small %>% dplyr::filter(peps_small$verdict == cond[y])
      sum(peps_smaller[,HLA])
    },numeric(1))
    names(cond_vals)<-cond
    return(cond_vals)
  })
  names(SB_vals)<-codes_named
  SB_vals<-(data.frame(SB_vals)+1)/(junc_counts+1)
  SB_vals<-SB_vals
  pheatmap::pheatmap(SB_vals,display_numbers=T,cluster_rows=F,cluster_cols=F,main=HLA,fontsize = 14)
  
  HLA<-data_colnames[i+1]
  WB_vals<-lapply(seq(length(codes)),function(x){
    peps_small<-peps %>% dplyr::filter(peps$error == codes[x])
    cond_vals<-vapply(seq(length(cond)),function(y){
      peps_smaller<-peps_small %>% dplyr::filter(peps_small$verdict == cond[y])
      sum(peps_smaller[,HLA])
    },numeric(1))
    names(cond_vals)<-cond
    return(cond_vals)
  })
  names(WB_vals)<-codes_named
  WB_vals<-(data.frame(WB_vals)+1)/(junc_counts+1)
  WB_vals<-WB_vals
  pheatmap::pheatmap(WB_vals,display_numbers=T,cluster_rows=F,cluster_cols=F,main=HLA,fontsize = 14)
}

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# differential immunogenicity--------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#


peps<-data_splicemutr_filt
codes<-c("2","3","4","tx")
codes_named<-c("multiple_stop","stop_not_end","start_not_beginning","translatable")
cond<-c("annotated","cryptic_fiveprime","cryptic_threeprime","cryptic_unanchored","novel annotated pair")
for (i in seq(1,14,by=2)){
  
  HLA<-data_colnames[i]
  
  junc_counts_tumor<-lapply(seq(length(codes)),function(x){
    peps_small<-peps %>% dplyr::filter(peps$error == codes[x] & peps$deltapsi > 0)
    cond_vals<-vapply(seq(length(cond)),function(y){
      peps_smaller<-peps_small %>% dplyr::filter(peps_small$verdict == cond[y])
      nrow(peps_smaller)
    },numeric(1))
    names(cond_vals)<-cond
    return(cond_vals)
  })
  names(junc_counts_tumor)<-codes_named
  junc_counts_tumor<-data.frame(junc_counts_tumor)
  
  junc_counts_normal<-lapply(seq(length(codes)),function(x){
    peps_small<-peps %>% dplyr::filter(peps$error == codes[x] & peps$deltapsi < 0)
    cond_vals<-vapply(seq(length(cond)),function(y){
      peps_smaller<-peps_small %>% dplyr::filter(peps_small$verdict == cond[y])
      nrow(peps_smaller)
    },numeric(1))
    names(cond_vals)<-cond
    return(cond_vals)
  })
  names(junc_counts_normal)<-codes_named
  junc_counts_normal<-data.frame(junc_counts_normal)
  
  SB_vals_tumor<-lapply(seq(length(codes)),function(x){
    peps_small<-peps %>% dplyr::filter(peps$error == codes[x] & peps$deltapsi > 0)
    cond_vals<-vapply(seq(length(cond)),function(y){
      peps_smaller<-peps_small %>% dplyr::filter(peps_small$verdict == cond[y])
      sum(peps_smaller[,HLA])
    },numeric(1))
    names(cond_vals)<-cond
    return(cond_vals)
  })
  names(SB_vals_tumor)<-codes_named
  SB_vals_tumor<-(data.frame(SB_vals_tumor)+1)/(junc_counts_tumor+1)
  
  SB_vals_normal<-lapply(seq(length(codes)),function(x){
    peps_small<-peps %>% dplyr::filter(peps$error == codes[x] & peps$deltapsi < 0)
    cond_vals<-vapply(seq(length(cond)),function(y){
      peps_smaller<-peps_small %>% dplyr::filter(peps_small$verdict == cond[y])
      sum(peps_smaller[,HLA])
    },numeric(1))
    names(cond_vals)<-cond
    return(cond_vals)
  })
  names(SB_vals_normal)<-codes_named
  SB_vals_normal<-(data.frame(SB_vals_normal)+1)/(junc_counts_normal+1)
  
  SB_vals<-SB_vals_tumor-SB_vals_normal
  pheatmap::pheatmap(SB_vals,display_numbers=T,cluster_rows=F,cluster_cols=F,main=HLA,fontsize = 14)
  
  HLA<-data_colnames[i+1]
  
  WB_vals_tumor<-lapply(seq(length(codes)),function(x){
    peps_small<-peps %>% dplyr::filter(peps$error == codes[x] & peps$deltapsi > 0)
    cond_vals<-vapply(seq(length(cond)),function(y){
      peps_smaller<-peps_small %>% dplyr::filter(peps_small$verdict == cond[y])
      sum(peps_smaller[,HLA])
    },numeric(1))
    names(cond_vals)<-cond
    return(cond_vals)
  })
  names(WB_vals_tumor)<-codes_named
  WB_vals_tumor<-(data.frame(WB_vals_tumor)+1)/(junc_counts_tumor+1)
  
  WB_vals_normal<-lapply(seq(length(codes)),function(x){
    peps_small<-peps %>% dplyr::filter(peps$error == codes[x] & peps$deltapsi < 0)
    cond_vals<-vapply(seq(length(cond)),function(y){
      peps_smaller<-peps_small %>% dplyr::filter(peps_small$verdict == cond[y])
      sum(peps_smaller[,HLA])
    },numeric(1))
    names(cond_vals)<-cond
    return(cond_vals)
  })
  names(WB_vals_normal)<-codes_named
  WB_vals_normal<-(data.frame(WB_vals_normal)+1)/(junc_counts_normal+1)
  
  WB_vals<-WB_vals_tumor-WB_vals_normal
  pheatmap::pheatmap(WB_vals,display_numbers=T,cluster_rows=F,cluster_cols=F,main=HLA,fontsize = 14)

}

#------------------------------------------------------------------------------#
# plotting SB and WB counts per junction, error type, and HLA type:


HLA<-data_colnames[3]

# From first frame, translatable based on junction type
codes<-c("2","3","4","tx")
codes_named<-c("multiple_stop","stop_not_end","start_not_beginning","translatable")
cond<-c("annotated","cryptic_fiveprime","cryptic_threeprime","cryptic_unanchored","novel annotated pair")

SB_vals_tumor<-lapply(seq(length(codes)),function(x){
  peps_small<-peps %>% dplyr::filter(peps$error == codes[x] & peps$deltapsi > 0)
  cond_vals<-vapply(seq(length(cond)),function(y){
    peps_smaller<-peps_small %>% dplyr::filter(peps_small$verdict == cond[y])
    sum(peps_smaller[,HLA])
  },numeric(1))
  names(cond_vals)<-cond
  return(cond_vals)
})
names(SB_vals_tumor)<-codes_named
SB_vals_tumor<-data.frame(SB_vals_tumor)

SB_vals_normal<-lapply(seq(length(codes)),function(x){
  peps_small<-peps %>% dplyr::filter(peps$error == codes[x] & peps$deltapsi < 0)
  cond_vals<-vapply(seq(length(cond)),function(y){
    peps_smaller<-peps_small %>% dplyr::filter(peps_small$verdict == cond[y])
    sum(peps_smaller[,HLA])
  },numeric(1))
  names(cond_vals)<-cond
  return(cond_vals)
})
names(SB_vals_normal)<-codes_named
SB_vals_normal<-data.frame(SB_vals_normal)

SB_vals<-log10(SB_vals_tumor+1)-log10(SB_vals_normal+1)
pheatmap::pheatmap(SB_vals,display_numbers=T,cluster_rows=F,cluster_cols=F,main=HLA)

HLA<-data_colnames[4]
WB_vals_tumor<-lapply(seq(length(codes)),function(x){
  peps_small<-peps %>% dplyr::filter(peps$error == codes[x] & peps$deltapsi > 0)
  cond_vals<-vapply(seq(length(cond)),function(y){
    peps_smaller<-peps_small %>% dplyr::filter(peps_small$verdict == cond[y])
    sum(peps_smaller[,HLA])
  },numeric(1))
  names(cond_vals)<-cond
  return(cond_vals)
})
names(WB_vals_tumor)<-codes_named
WB_vals_tumor<-data.frame(WB_vals_tumor)

WB_vals_normal<-lapply(seq(length(codes)),function(x){
  peps_small<-peps %>% dplyr::filter(peps$error == codes[x] & peps$deltapsi < 0)
  cond_vals<-vapply(seq(length(cond)),function(y){
    peps_smaller<-peps_small %>% dplyr::filter(peps_small$verdict == cond[y])
    sum(peps_smaller[,HLA])
  },numeric(1))
  names(cond_vals)<-cond
  return(cond_vals)
})
names(WB_vals_normal)<-codes_named
WB_vals_normal<-data.frame(WB_vals_normal)

WB_vals<-log10(WB_vals_tumor+1)-log10(WB_vals_normal+1)
pheatmap::pheatmap(WB_vals,display_numbers=T,cluster_rows=F,cluster_cols=F,main=HLA)

#------------------------------------------------------------------------------#
# plotting SB and WB counts per junction, error type, and HLA type:


HLA<-data_colnames[5]

# From first frame, translatable based on junction type
codes<-c("2","3","4","tx")
codes_named<-c("multiple_stop","stop_not_end","start_not_beginning","translatable")
cond<-c("annotated","cryptic_fiveprime","cryptic_threeprime","cryptic_unanchored","novel annotated pair")

SB_vals_tumor<-lapply(seq(length(codes)),function(x){
  peps_small<-peps %>% dplyr::filter(peps$error == codes[x] & peps$deltapsi > 0)
  cond_vals<-vapply(seq(length(cond)),function(y){
    peps_smaller<-peps_small %>% dplyr::filter(peps_small$verdict == cond[y])
    sum(peps_smaller[,HLA])
  },numeric(1))
  names(cond_vals)<-cond
  return(cond_vals)
})
names(SB_vals_tumor)<-codes_named
SB_vals_tumor<-data.frame(SB_vals_tumor)

SB_vals_normal<-lapply(seq(length(codes)),function(x){
  peps_small<-peps %>% dplyr::filter(peps$error == codes[x] & peps$deltapsi < 0)
  cond_vals<-vapply(seq(length(cond)),function(y){
    peps_smaller<-peps_small %>% dplyr::filter(peps_small$verdict == cond[y])
    sum(peps_smaller[,HLA])
  },numeric(1))
  names(cond_vals)<-cond
  return(cond_vals)
})
names(SB_vals_normal)<-codes_named
SB_vals_normal<-data.frame(SB_vals_normal)

SB_vals<-log10(SB_vals_tumor+1)-log10(SB_vals_normal+1)
pheatmap::pheatmap(SB_vals,display_numbers=T,cluster_rows=F,cluster_cols=F,main=HLA)

HLA<-data_colnames[6]

WB_vals_tumor<-lapply(seq(length(codes)),function(x){
  peps_small<-peps %>% dplyr::filter(peps$error == codes[x] & peps$deltapsi > 0)
  cond_vals<-vapply(seq(length(cond)),function(y){
    peps_smaller<-peps_small %>% dplyr::filter(peps_small$verdict == cond[y])
    sum(peps_smaller[,HLA])
  },numeric(1))
  names(cond_vals)<-cond
  return(cond_vals)
})
names(WB_vals_tumor)<-codes_named
WB_vals_tumor<-data.frame(WB_vals_tumor)

WB_vals_normal<-lapply(seq(length(codes)),function(x){
  peps_small<-peps %>% dplyr::filter(peps$error == codes[x] & peps$deltapsi < 0)
  cond_vals<-vapply(seq(length(cond)),function(y){
    peps_smaller<-peps_small %>% dplyr::filter(peps_small$verdict == cond[y])
    sum(peps_smaller[,HLA])
  },numeric(1))
  names(cond_vals)<-cond
  return(cond_vals)
})
names(WB_vals_normal)<-codes_named
WB_vals_normal<-data.frame(WB_vals_normal)

WB_vals<-log10(WB_vals_tumor+1)-log10(WB_vals_normal+1)
pheatmap::pheatmap(WB_vals,display_numbers=T,cluster_rows=F,cluster_cols=F,main=HLA)

#------------------------------------------------------------------------------#
# plotting SB and WB counts per junction, error type, and HLA type:


HLA<-data_colnames[7]

# From first frame, translatable based on junction type
codes<-c("2","3","4","tx")
codes_named<-c("multiple_stop","stop_not_end","start_not_beginning","translatable")
cond<-c("annotated","cryptic_fiveprime","cryptic_threeprime","cryptic_unanchored","novel annotated pair")

SB_vals_tumor<-lapply(seq(length(codes)),function(x){
  peps_small<-peps %>% dplyr::filter(peps$error == codes[x] & peps$deltapsi > 0)
  cond_vals<-vapply(seq(length(cond)),function(y){
    peps_smaller<-peps_small %>% dplyr::filter(peps_small$verdict == cond[y])
    sum(peps_smaller[,HLA])
  },numeric(1))
  names(cond_vals)<-cond
  return(cond_vals)
})
names(SB_vals_tumor)<-codes_named
SB_vals_tumor<-data.frame(SB_vals_tumor)

SB_vals_normal<-lapply(seq(length(codes)),function(x){
  peps_small<-peps %>% dplyr::filter(peps$error == codes[x] & peps$deltapsi < 0)
  cond_vals<-vapply(seq(length(cond)),function(y){
    peps_smaller<-peps_small %>% dplyr::filter(peps_small$verdict == cond[y])
    sum(peps_smaller[,HLA])
  },numeric(1))
  names(cond_vals)<-cond
  return(cond_vals)
})
names(SB_vals_normal)<-codes_named
SB_vals_normal<-data.frame(SB_vals_normal)

SB_vals<-log10(SB_vals_tumor+1)-log10(SB_vals_normal+1)
pheatmap::pheatmap(SB_vals,display_numbers=T,cluster_rows=F,cluster_cols=F,main=HLA)

HLA<-data_colnames[8]

WB_vals_tumor<-lapply(seq(length(codes)),function(x){
  peps_small<-peps %>% dplyr::filter(peps$error == codes[x] & peps$deltapsi > 0)
  cond_vals<-vapply(seq(length(cond)),function(y){
    peps_smaller<-peps_small %>% dplyr::filter(peps_small$verdict == cond[y])
    sum(peps_smaller[,HLA])
  },numeric(1))
  names(cond_vals)<-cond
  return(cond_vals)
})
names(WB_vals_tumor)<-codes_named
WB_vals_tumor<-data.frame(WB_vals_tumor)

WB_vals_normal<-lapply(seq(length(codes)),function(x){
  peps_small<-peps %>% dplyr::filter(peps$error == codes[x] & peps$deltapsi < 0)
  cond_vals<-vapply(seq(length(cond)),function(y){
    peps_smaller<-peps_small %>% dplyr::filter(peps_small$verdict == cond[y])
    sum(peps_smaller[,HLA])
  },numeric(1))
  names(cond_vals)<-cond
  return(cond_vals)
})
names(WB_vals_normal)<-codes_named
WB_vals_normal<-data.frame(WB_vals_normal)

WB_vals<-log10(WB_vals_tumor+1)-log10(WB_vals_normal+1)
pheatmap::pheatmap(WB_vals,display_numbers=T,cluster_rows=F,cluster_cols=F,main=HLA)

#------------------------------------------------------------------------------#
# plotting SB and WB counts per junction, error type, and HLA type:


HLA<-data_colnames[9]

# From first frame, translatable based on junction type
codes<-c("2","3","4","tx")
codes_named<-c("multiple_stop","stop_not_end","start_not_beginning","translatable")
cond<-c("annotated","cryptic_fiveprime","cryptic_threeprime","cryptic_unanchored","novel annotated pair")

SB_vals_tumor<-lapply(seq(length(codes)),function(x){
  peps_small<-peps %>% dplyr::filter(peps$error == codes[x] & peps$deltapsi > 0)
  cond_vals<-vapply(seq(length(cond)),function(y){
    peps_smaller<-peps_small %>% dplyr::filter(peps_small$verdict == cond[y])
    sum(peps_smaller[,HLA])
  },numeric(1))
  names(cond_vals)<-cond
  return(cond_vals)
})
names(SB_vals_tumor)<-codes_named
SB_vals_tumor<-data.frame(SB_vals_tumor)

SB_vals_normal<-lapply(seq(length(codes)),function(x){
  peps_small<-peps %>% dplyr::filter(peps$error == codes[x] & peps$deltapsi < 0)
  cond_vals<-vapply(seq(length(cond)),function(y){
    peps_smaller<-peps_small %>% dplyr::filter(peps_small$verdict == cond[y])
    sum(peps_smaller[,HLA])
  },numeric(1))
  names(cond_vals)<-cond
  return(cond_vals)
})
names(SB_vals_normal)<-codes_named
SB_vals_normal<-data.frame(SB_vals_normal)

SB_vals<-log10(SB_vals_tumor+1)-log10(SB_vals_normal+1)
pheatmap::pheatmap(SB_vals,display_numbers=T,cluster_rows=F,cluster_cols=F,main=HLA)

HLA<-data_colnames[10]

WB_vals_tumor<-lapply(seq(length(codes)),function(x){
  peps_small<-peps %>% dplyr::filter(peps$error == codes[x] & peps$deltapsi > 0)
  cond_vals<-vapply(seq(length(cond)),function(y){
    peps_smaller<-peps_small %>% dplyr::filter(peps_small$verdict == cond[y])
    sum(peps_smaller[,HLA])
  },numeric(1))
  names(cond_vals)<-cond
  return(cond_vals)
})
names(WB_vals_tumor)<-codes_named
WB_vals_tumor<-data.frame(WB_vals_tumor)

WB_vals_normal<-lapply(seq(length(codes)),function(x){
  peps_small<-peps %>% dplyr::filter(peps$error == codes[x] & peps$deltapsi < 0)
  cond_vals<-vapply(seq(length(cond)),function(y){
    peps_smaller<-peps_small %>% dplyr::filter(peps_small$verdict == cond[y])
    sum(peps_smaller[,HLA])
  },numeric(1))
  names(cond_vals)<-cond
  return(cond_vals)
})
names(WB_vals_normal)<-codes_named
WB_vals_normal<-data.frame(WB_vals_normal)

WB_vals<-log10(WB_vals_tumor+1)-log10(WB_vals_normal+1)
pheatmap::pheatmap(WB_vals,display_numbers=T,cluster_rows=F,cluster_cols=F,main=HLA)
#------------------------------------------------------------------------------#
# plotting SB and WB counts per junction, error type, and HLA type:


HLA<-data_colnames[11]

# From first frame, translatable based on junction type
codes<-c("2","3","4","tx")
codes_named<-c("multiple_stop","stop_not_end","start_not_beginning","translatable")
cond<-c("annotated","cryptic_fiveprime","cryptic_threeprime","cryptic_unanchored","novel annotated pair")

SB_vals_tumor<-lapply(seq(length(codes)),function(x){
  peps_small<-peps %>% dplyr::filter(peps$error == codes[x] & peps$deltapsi > 0)
  cond_vals<-vapply(seq(length(cond)),function(y){
    peps_smaller<-peps_small %>% dplyr::filter(peps_small$verdict == cond[y])
    sum(peps_smaller[,HLA])
  },numeric(1))
  names(cond_vals)<-cond
  return(cond_vals)
})
names(SB_vals_tumor)<-codes_named
SB_vals_tumor<-data.frame(SB_vals_tumor)

SB_vals_normal<-lapply(seq(length(codes)),function(x){
  peps_small<-peps %>% dplyr::filter(peps$error == codes[x] & peps$deltapsi < 0)
  cond_vals<-vapply(seq(length(cond)),function(y){
    peps_smaller<-peps_small %>% dplyr::filter(peps_small$verdict == cond[y])
    sum(peps_smaller[,HLA])
  },numeric(1))
  names(cond_vals)<-cond
  return(cond_vals)
})
names(SB_vals_normal)<-codes_named
SB_vals_normal<-data.frame(SB_vals_normal)

SB_vals<-log10(SB_vals_tumor+1)-log10(SB_vals_normal+1)
pheatmap::pheatmap(SB_vals,display_numbers=T,cluster_rows=F,cluster_cols=F,main=HLA)

HLA<-data_colnames[12]

WB_vals_tumor<-lapply(seq(length(codes)),function(x){
  peps_small<-peps %>% dplyr::filter(peps$error == codes[x] & peps$deltapsi > 0)
  cond_vals<-vapply(seq(length(cond)),function(y){
    peps_smaller<-peps_small %>% dplyr::filter(peps_small$verdict == cond[y])
    sum(peps_smaller[,HLA])
  },numeric(1))
  names(cond_vals)<-cond
  return(cond_vals)
})
names(WB_vals_tumor)<-codes_named
WB_vals_tumor<-data.frame(WB_vals_tumor)

WB_vals_normal<-lapply(seq(length(codes)),function(x){
  peps_small<-peps %>% dplyr::filter(peps$error == codes[x] & peps$deltapsi < 0)
  cond_vals<-vapply(seq(length(cond)),function(y){
    peps_smaller<-peps_small %>% dplyr::filter(peps_small$verdict == cond[y])
    sum(peps_smaller[,HLA])
  },numeric(1))
  names(cond_vals)<-cond
  return(cond_vals)
})
names(WB_vals_normal)<-codes_named
WB_vals_normal<-data.frame(WB_vals_normal)

WB_vals<-log10(WB_vals_tumor+1)-log10(WB_vals_normal+1)
pheatmap::pheatmap(WB_vals,display_numbers=T,cluster_rows=F,cluster_cols=F,main=HLA)

#------------------------------------------------------------------------------#
# plotting SB and WB counts per junction, error type, and HLA type:


HLA<-data_colnames[13]

# From first frame, translatable based on junction type
codes<-c("2","3","4","tx")
codes_named<-c("multiple_stop","stop_not_end","start_not_beginning","translatable")
cond<-c("annotated","cryptic_fiveprime","cryptic_threeprime","cryptic_unanchored","novel annotated pair")

SB_vals_tumor<-lapply(seq(length(codes)),function(x){
  peps_small<-peps %>% dplyr::filter(peps$error == codes[x] & peps$deltapsi > 0)
  cond_vals<-vapply(seq(length(cond)),function(y){
    peps_smaller<-peps_small %>% dplyr::filter(peps_small$verdict == cond[y])
    sum(peps_smaller[,HLA])
  },numeric(1))
  names(cond_vals)<-cond
  return(cond_vals)
})
names(SB_vals_tumor)<-codes_named
SB_vals_tumor<-data.frame(SB_vals_tumor)

SB_vals_normal<-lapply(seq(length(codes)),function(x){
  peps_small<-peps %>% dplyr::filter(peps$error == codes[x] & peps$deltapsi < 0)
  cond_vals<-vapply(seq(length(cond)),function(y){
    peps_smaller<-peps_small %>% dplyr::filter(peps_small$verdict == cond[y])
    sum(peps_smaller[,HLA])
  },numeric(1))
  names(cond_vals)<-cond
  return(cond_vals)
})
names(SB_vals_normal)<-codes_named
SB_vals_normal<-data.frame(SB_vals_normal)

SB_vals<-log10(SB_vals_tumor+1)-log10(SB_vals_normal+1)
pheatmap::pheatmap(SB_vals,display_numbers=T,cluster_rows=F,cluster_cols=F,main=HLA)

HLA<-data_colnames[14]

WB_vals_tumor<-lapply(seq(length(codes)),function(x){
  peps_small<-peps %>% dplyr::filter(peps$error == codes[x] & peps$deltapsi > 0)
  cond_vals<-vapply(seq(length(cond)),function(y){
    peps_smaller<-peps_small %>% dplyr::filter(peps_small$verdict == cond[y])
    sum(peps_smaller[,HLA])
  },numeric(1))
  names(cond_vals)<-cond
  return(cond_vals)
})
names(WB_vals_tumor)<-codes_named
WB_vals_tumor<-data.frame(WB_vals_tumor)

WB_vals_normal<-lapply(seq(length(codes)),function(x){
  peps_small<-peps %>% dplyr::filter(peps$error == codes[x] & peps$deltapsi < 0)
  cond_vals<-vapply(seq(length(cond)),function(y){
    peps_smaller<-peps_small %>% dplyr::filter(peps_small$verdict == cond[y])
    sum(peps_smaller[,HLA])
  },numeric(1))
  names(cond_vals)<-cond
  return(cond_vals)
})
names(WB_vals_normal)<-codes_named
WB_vals_normal<-data.frame(WB_vals_normal)

WB_vals<-log10(WB_vals_tumor+1)-log10(WB_vals_normal+1)
pheatmap::pheatmap(WB_vals,display_numbers=T,cluster_rows=F,cluster_cols=F,main=HLA)



