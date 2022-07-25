# The purpose of this script is to process the immunogenicity with respect to the genotype per patient

# input the alleles and parse them
# organize the alleles such that they are in the file format
# parse out those rows that correspond to tumor patients
# per patient, read in the respective mhc alleles and possible pairs of alleles as apppropriate
  # A,B,C; DRA:DRB,DPA:DPB,DQA:DQB
# create all potential genotypes: 1 from A,B,C; 1 from D*A:D*B pairing
# within this genotype determine those rows that are shared amongst all to be immunogenic
# do this for all genotypes
# take the union of these immunogenic rows per pairing as the set of potential immunogenic transcripts for the patient


#------------------------------------------------------------------------------#
# loading libraries

library(stringr)
library(ggplot2)

#------------------------------------------------------------------------------#
# internal functions

arcas_to_file<-function(hla){
  hla_breaks<-str_locate_all(hla,":")[[1]]
  if (length(hla_breaks)>2){
    hla_small<-substr(hla,1,hla_breaks[2]-1)
  } else {
    hla_small<-hla
  }
  hla_small<-str_replace(hla_small,"[*]","")
  hla_small<-str_replace(hla_small,"[:]","-")
  hla_small<-sprintf("HLA-%s",hla_small)
}

separate_alleles<-function(alleles_tum){
  
  ABC<-c("HLA-A","HLA-B","HLA-C")
  DPA<-"HLA-DPA"
  DPB<-"HLA-DPB"
  DQA<-"HLA-DQA"
  DQB<-"HLA-DQB"
  DRA<-"HLA-DRA"
  DRB<-"HLA-DRB"
  
  ABCs<-alleles_tum[vapply(alleles_tum,function(allele){
    if (any(str_detect(allele,ABC))){return(T)}else{return(F)}
  },logical(1))]
  if (length(ABCs)==0){ABCs=""}
  DPAs<-alleles_tum[vapply(alleles_tum,function(allele){
    if (any(str_detect(allele,DPA))){return(T)}else{return(F)}
  },logical(1))]
  if (length(DPAs)==0){DPAs=""}
  DPBs<-alleles_tum[vapply(alleles_tum,function(allele){
    if (any(str_detect(allele,DPB))){return(T)}else{return(F)}
  },logical(1))]
  if (length(DPBs)==0){DPBs=""}
  DQAs<-alleles_tum[vapply(alleles_tum,function(allele){
    if (any(str_detect(allele,DQA))){return(T)}else{return(F)}
  },logical(1))]
  if (length(DQAs)==0){DQAs=""}
  DQBs<-alleles_tum[vapply(alleles_tum,function(allele){
    if (any(str_detect(allele,DQB))){return(T)}else{return(F)}
  },logical(1))]
  if (length(DQBs)==0){DQBs=""}
  DRAs<-alleles_tum[vapply(alleles_tum,function(allele){
    if (any(str_detect(allele,DRA))){return(T)}else{return(F)}
  },logical(1))]
  if (length(DRAs)==0){DRAs=""}
  DRBs<-alleles_tum[vapply(alleles_tum,function(allele){
    if (any(str_detect(allele,DRB))){return(T)}else{return(F)}
  },logical(1))]
  if (length(DRBs)==0){DRBs=""}
  alleles_split<-list(ABCs,DPAs,DPBs,DQAs,DQBs,DRAs,DRBs)
  names(alleles_split)<-c("ABC","DPA","DPB","DQA","DQB","DRA","DBA")
  return(alleles_split)
}

#------------------------------------------------------------------------------#
# reading in the alleles

allele_file<-"/media/theron/My_Passport/head_and_neck_DARIA/data/leafcutter_10_21_2020/leafcutter_introns_split/alleles/alleles.txt"
alleles<-read.table(allele_file,header=F,sep="\t")
colnames(alleles)<-"allele_set"

#------------------------------------------------------------------------------#
# organize the alleles such that they are in the file format

alleles_rep<-c("A","B","C","DPA1","DPB1","DQA1","DQB1","DRA","DRB1")
alleles_split<-str_split(alleles$allele_set,"[-]")
alleles_clean<-lapply(seq(length(alleles_split)),function(val){
  types<-str_split(str_replace_all(gsub("\\[|\\]|", "", alleles_split[[val]])," ",""),"[,]")
  names(types)<-alleles_rep
  types<-vapply(unlist(types),function(allele){
    arcas_to_file(allele)
  },character(1))
  return(unname(types))
})

#------------------------------------------------------------------------------#
# reading in the patients files

key_file<-"/media/theron/My_Passport/head_and_neck_DARIA/data/leafcutter_10_21_2020/keys.txt"
keys<-read.table(key_file,header=F)
colnames(keys)<-c("original","ID","type")

# reading in the alleles labels
samples_file<-"/media/theron/My_Passport/head_and_neck_DARIA/data/leafcutter_10_21_2020/leafcutter_introns_split/alleles/samples.txt"
samples<-read.table(samples_file)
colnames(samples)<-"samples"

#------------------------------------------------------------------------------#
# determining which allele sets are Tumor or Normal

samples_not_in<-samples$sample[which(!(samples$sample %in% keys$original))]
tumor_or_normal<-vapply(samples$samples,function(sample){
  if (sample %in% samples_not_in){
    return("0")
  } else {
    return(keys$type[which(keys$original==sample)])
  }
},character(1))

samples$type<-tumor_or_normal

tumors<-which(samples$type=="T")
normals<-which(samples$type=="N")

#------------------------------------------------------------------------------#
# extracting those mhc sets that are tumors

tumor_alleles_clean<-alleles_clean[tumors]
IDS<-vapply(samples$samples[tumors],function(samp){
  keys$ID[which(keys$original==samp)]
},character(1))
names(tumor_alleles_clean)<-IDS

#------------------------------------------------------------------------------#
# per patient, read in the appropriate immunogenic rows per allele

mhc_I_dir<-"/media/theron/My_Passport/head_and_neck_DARIA/data/leafcutter_10_21_2020/leafcutter_introns_split/predictions_1/rows_per_allele"
mhc_II_dir<-"/media/theron/My_Passport/head_and_neck_DARIA/data/leafcutter_10_21_2020/leafcutter_introns_split/predictions_2/rows_per_allele"

ABC<-c("HLA-A","HLA-B","HLA-C")
DPA<-"HLA-DPA"
DPB<-"HLA-DPB"
DQA<-"HLA-DQA"
DQB<-"HLA-DQB"
DRA<-"HLA-DRA"
DRB<-"HLA-DRB"

immunogenic_rows<-lapply(names(tumor_alleles_clean),function(ID){
  alleles_tum<-tumor_alleles_clean[[ID]]
  alleles_tum_split<-separate_alleles(alleles_tum)
  alleles_pairings<-unique(expand.grid(alleles_tum_split[[1]],
                                alleles_tum_split[[2]],
                                alleles_tum_split[[3]],
                                alleles_tum_split[[4]],
                                alleles_tum_split[[5]],
                                alleles_tum_split[[6]],
                                alleles_tum_split[[7]]))
  ABC_rows<-c()
  # DPA_rows<-c()
  DPB_rows<-c()
  # DQA_rows<-c()
  DQB_rows<-c()
  # DRA_rows<-c()
  DRB_rows<-c()
  for (i in alleles_tum_split[[1]]){
    for (j in seq(9,13)){
      ABC_file<-sprintf("%s/%s_tx_dict_%d.pickle.rows.txt",mhc_I_dir,i,j)
      if (file.exists(ABC_file) & file.info(ABC_file)$size != 0){
        ABC_rows<-c(ABC_rows,read.table(ABC_file)[,1])
      }
    }
  }
  for (i in alleles_tum_split[[3]]){
    for (j in seq(13,25)){
      # DPA_file<-sprintf("%s/%s_tx_dict_%d.pickle.rows.txt",mhc_II_dir,alleles_pairings[i,2],j)
      # if (file.exists(DPA_file)){
      #   DPA_rows<-c(DPA_rows,read.table(DPA_file)[,1])
      # }
      DPB_file<-sprintf("%s/%s_tx_dict_%d.pickle.rows.txt",mhc_II_dir,alleles_pairings[i,3],j)
      if (file.exists(DPB_file)){
        DPB_rows<-c(DPB_rows,read.table(DPB_file)[,1])
      }
    }
  }
  for (i in alleles_tum_split[[3]]){
    for (j in seq(13,25)){
      # DQA_file<-sprintf("%s/%s_tx_dict_%d.pickle.rows.txt",mhc_II_dir,alleles_pairings[i,4],j)
      # if (file.exists(DQA_file)){
      #   DQA_rows<-c(DQA_rows,read.table(DQA_file)[,1])
      # }
      DQB_file<-sprintf("%s/%s_tx_dict_%d.pickle.rows.txt",mhc_II_dir,alleles_pairings[i,5],j)
      if (file.exists(DQB_file)){
        DQB_rows<-c(DQB_rows,read.table(DQB_file)[,1])
      }
    }
  }
  for (i in alleles_tum_split[[3]]){
    for (j in seq(13,25)){
      # DRA_file<-sprintf("%s/%s_tx_dict_%d.pickle.rows.txt",mhc_II_dir,alleles_pairings[i,6],j)
      # if (file.exists(DRA_file)){
      #   DRA_rows<-c(DRA_rows,read.table(DRA_file)[,1])
      # }
      DRB_file<-sprintf("%s/%s_tx_dict_%d.pickle.rows.txt",mhc_II_dir,alleles_pairings[i,7],j)
      if (file.exists(DRB_file)){
        DRB_rows<-c(DRB_rows,read.table(DRB_file)[,1])
      }
    }
  }
  ABC_rows<-unique(ABC_rows)
  # DPA_rows<-unique(DPA_rows)
  DPB_rows<-unique(DPB_rows)
  # DQA_rows<-unique(DQA_rows)
  DQB_rows<-unique(DQB_rows)
  # DRA_rows<-unique(DRA_rows)
  DRB_rows<-unique(DRB_rows)
  if (length(DPB_rows)>0){
    rows_1<-ABC_rows[which(ABC_rows %in% DPB_rows)]
  } else {rows_1<-ABC_rows}
  if (length(DQB_rows)>0){
    rows_2<-ABC_rows[which(ABC_rows %in% DQB_rows)]
  } else {rows_2<-ABC_rows}
  if (length(DRB_rows)>0){
    rows_3<-ABC_rows[which(ABC_rows %in% DRB_rows)]
  } else {rows_3<-ABC_rows}
  return(unique(c(rows_1,rows_2,rows_3)))
})

saveRDS(immunogenic_rows,
        file="/media/theron/My_Passport/head_and_neck_DARIA/data/leafcutter_10_21_2020/leafcutter_introns_split/predictions_1/immunogenic_rows.rds")

#------------------------------------------------------------------------------#
# read in the conglomerate data_splicemutr information
immunogenic_rows<-readRDS("/media/theron/My_Passport/head_and_neck_DARIA/data/leafcutter_10_21_2020/leafcutter_introns_split/predictions_1/immunogenic_rows.rds")
splice_file<-"/media/theron/My_Passport/head_and_neck_DARIA/data/leafcutter_10_21_2020/leafcutter_introns_split/data_splicemutr_dup.txt"
splicemutr_dat<-read.table(splice_file,skip = 1)

all_rows<-as.data.frame(unlist(immunogenic_rows))
colnames(all_rows)<-"rows"
all_rows_tab<-table(unlist(immunogenic_rows))

png("/media/theron/My_Passport/head_and_neck_DARIA/data/leafcutter_10_21_2020/leafcutter_introns_split/immun_hist_2.png",width=700,height=700)
  ggplot(all_rows,aes(x=rows))+geom_histogram(bins=200000)
  #barplot(all_rows_tab)
dev.off()
