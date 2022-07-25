# This script processes the mhc alleles json file

#------------------------------------------------------------------------------#
# loading libraries

library(stringr)
library(ggplot2)

#------------------------------------------------------------------------------#
# internal functions

arcas_to_mhcnuggets<-function(hla){
  hla_breaks<-str_locate_all(hla,":")[[1]]
  if (length(hla_breaks)>2){
    hla_small<-substr(hla,1,hla_breaks[2]-1)
  } else {
    hla_small<-hla
  }
  hla_small<-str_replace(hla_small,"[*]","")
  hla_small<-sprintf("HLA-%s",hla_small)
}


#------------------------------------------------------------------------------#
# reading in the allele file

allele_file<-"/media/theron/My_Passport/head_and_neck_DARIA/data/leafcutter_10_21_2020/leafcutter_introns_split/alleles/alleles.txt"
alleles<-read.table(allele_file,header=F,sep="\t")
colnames(alleles)<-"allele_set"

alleles_rep<-c("A","B","C","DPA1","DPB1","DQA1","DQB1","DRA","DRB1")
alleles_split<-str_split(alleles$allele_set,"[-]")
alleles_clean<-lapply(seq(length(alleles_split)),function(val){
  types<-str_split(str_replace_all(gsub("\\[|\\]", "", alleles_split[[val]])," ",""),"[,]")
  names(types)<-alleles_rep
  types<-vapply(unlist(types),function(allele){
    arcas_to_mhcnuggets(allele)
  },character(1))
  return(unname(types))
})
alleles_clean_orig<-lapply(seq(length(alleles_split)),function(val){
  types<-str_split(str_replace_all(gsub("\\[|\\]|", "", alleles_split[[val]])," ",""),"[,]")
  names(types)<-alleles_rep
  return(unname(types))
})

alleles_unique<-unique(unlist(alleles_clean))
alleles_unique_orig<-unique(unlist(alleles_clean_orig))

#------------------------------------------------------------------------------#
# checking to see which alleles can be processed by mhcnuggets currently

mhc_I<-read.table("/media/theron/My_Passport/splicemutr/MHC_I_alleles.txt")
mhc_II<-read.table("/media/theron/My_Passport/splicemutr/MHC_II_alleles.txt")
mhc_tot<-rbind(mhc_I,mhc_II)
colnames(mhc_tot)<-"mhc"

num_to_process<-length(which(alleles_unique %in% mhc_tot$mhc))
mhc_I_to_process<-alleles_unique[which(alleles_unique %in% mhc_I$V1)]
mhc_II_to_process<-alleles_unique[which(alleles_unique %in% mhc_II$V1)]

allele_dir<-"/media/theron/My_Passport/head_and_neck_DARIA/data/leafcutter_10_21_2020/leafcutter_introns_split/alleles"

write.table(mhc_I_to_process,
            file=sprintf("%s/%s",allele_dir,"mhc_I_to_process.txt"),
            quote=F,col.names=F,row.names=F)

write.table(mhc_II_to_process,
            file=sprintf("%s/%s",allele_dir,"mhc_II_to_process.txt"),
            quote=F,col.names=F,row.names=F)

#------------------------------------------------------------------------------#
# reading in the key tumor file

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

#------------------------------------------------------------------------------#
# generating alleles metrics

alleles_unique_counts<-vapply(alleles_unique,function(allele){
  logs<-vapply(seq(length(alleles_clean)),function(num){
    if (allele %in% alleles_clean[[num]]){
      return(T)
    } else {
      return(F)
    }
  },logical(1))
  length(which(logs==T))
},numeric(1))
alleles_unique_counts<-data.frame(sort(alleles_unique_counts,decreasing = T))
alleles_unique_counts$HLA<-rownames(alleles_unique_counts)
colnames(alleles_unique_counts)<-c("counts","HLA")
ggplot(alleles_unique_counts[1:20,],aes(x=reorder(HLA ,-counts),y=counts))+
  geom_bar(stat="identity")+
  geom_text(stat='identity', aes(label=counts), vjust=-1)+
  theme(axis.text.x = element_text(angle = 270, vjust = 0.5)) +
  xlab("") + labs(title=sprintf("All Samples (%s)",length(tumors)+length(normals)))

# alleles_unique_orig_counts<-vapply(alleles_unique_orig,function(allele){
#   logs<-vapply(seq(length(alleles_clean_orig)),function(num){
#     if (allele %in% alleles_clean_orig[[num]]){
#       return(T)
#     } else {
#       return(F)
#     }
#   },logical(1))
#   length(which(logs==T))
# },numeric(1))
# alleles_unique_orig_counts<-data.frame(sort(alleles_unique_orig_counts,decreasing = T))
# alleles_unique_orig_counts$HLA<-rownames(alleles_unique_orig_counts)
# colnames(alleles_unique_orig_counts)<-c("counts","HLA")
# ggplot(alleles_unique_orig_counts[1:20,],aes(x=reorder(HLA ,-counts),y=counts))+
#   geom_bar(stat="identity")+
#   geom_text(stat='identity', aes(label=counts), vjust=-1)+
#   theme(axis.text.x = element_text(angle = 270, vjust = 0.5)) +
#   xlab("")

# generating tumor alleles metrics

tumors<-which(samples$type=="T")
normals<-which(samples$type=="N")

alleles_clean_tum<-alleles_clean[tumors]
alleles_clean_norm<-alleles_clean[normals]

alleles_unique_counts_tumor<-vapply(alleles_unique,function(allele){
  logs<-vapply(seq(length(alleles_clean_tum)),function(num){
    if (allele %in% alleles_clean_tum[[num]]){
      return(T)
    } else {
      return(F)
    }
  },logical(1))
  length(which(logs==T))
},numeric(1))
alleles_unique_counts_tumor<-data.frame(sort(alleles_unique_counts_tumor,decreasing = T))
alleles_unique_counts_tumor$HLA<-rownames(alleles_unique_counts_tumor)
colnames(alleles_unique_counts_tumor)<-c("counts","HLA")
ggplot(alleles_unique_counts_tumor[1:20,],aes(x=reorder(HLA ,-counts),y=counts))+
  geom_bar(stat="identity")+
  geom_text(stat='identity', aes(label=counts), vjust=-1)+
  theme(axis.text.x = element_text(angle = 270, vjust = 0.5)) +
  xlab("") + labs(title=sprintf("Tumor Samples (%s)",length(tumors)))

alleles_unique_counts_normal<-vapply(alleles_unique,function(allele){
  logs<-vapply(seq(length(alleles_clean_norm)),function(num){
    if (allele %in% alleles_clean_norm[[num]]){
      return(T)
    } else {
      return(F)
    }
  },logical(1))
  length(which(logs==T))
},numeric(1))
alleles_unique_counts_normal<-data.frame(sort(alleles_unique_counts_normal,decreasing = T))
alleles_unique_counts_normal$HLA<-rownames(alleles_unique_counts_normal)
colnames(alleles_unique_counts_normal)<-c("counts","HLA")
ggplot(alleles_unique_counts_normal[1:20,],aes(x=reorder(HLA ,-counts),y=counts))+
  geom_bar(stat="identity")+
  geom_text(stat='identity', aes(label=counts), vjust=-1)+
  theme(axis.text.x = element_text(angle = 270, vjust = 0.5)) +
  xlab("") + labs(title=sprintf("Normal Samples (%s)",length(normals)))

alleles_clean_orig_tum<-alleles_clean_orig[tumors]
alleles_clean_orig_norm<-alleles_clean_orig[normals]

alleles_unique_orig_counts_tum<-vapply(alleles_unique_orig,function(allele){
  logs<-vapply(seq(length(alleles_clean_orig_tum)),function(num){
    if (allele %in% alleles_clean_orig_tum[[num]]){
      return(T)
    } else {
      return(F)
    }
  },logical(1))
  length(which(logs==T))
},numeric(1))

alleles_unique_orig_counts_norm<-vapply(alleles_unique_orig,function(allele){
  logs<-vapply(seq(length(alleles_clean_orig_norm)),function(num){
    if (allele %in% alleles_clean_orig_norm[[num]]){
      return(T)
    } else {
      return(F)
    }
  },logical(1))
  length(which(logs==T))
},numeric(1))

