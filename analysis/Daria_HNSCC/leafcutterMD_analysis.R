library(stringr)

# loading in the outlier analysis results used in Daria's paper
outlier_dir<-"/media/theron/My_Passport/data/leafcutter_08_24_2020"
out_calls<-read.csv(sprintf("%s/%s",outlier_dir,"outlier_calls.csv"))
out_calls_meta<-read.csv(sprintf("%s/%s",outlier_dir,"outlier_calls_meta.csv"))
out_calls_novel<-read.csv(sprintf("%s/%s",outlier_dir,"outlier_calls_novel.csv"))

key<-read.csv(sprintf("%s/%s",outlier_dir,"2015-09-15_DiscoveryCohort-Final_ID-RNA_nolab.csv"))
nme<-colnames(out_calls)[5:51]
nme_split<-vapply(str_split(nme,"\\."),function(x){
  fir<-key[,2] == x[2]
  if (any(fir == TRUE)){
    as.character(key[fir,1])
  } else {
    as.character(key[key[,3] == x[2],1])
  }
},character(1))
for (i in str_split(nme,"\\.")){
  print(sprintf("%s %s",as.character(any((key[,2]==i[2])==TRUE)),as.character(any((key[,3]==i[2])==TRUE))))
}

# loading in the outlier analysis results done using leafcutterMD
files<-read.table(sprintf("%s/%s",outlier_dir,"files.txt"))
file_split<-str_split(c(files[,1]),".junc")
file<-vapply(seq_len(length(file_split)),function(x){
  file<-file_split[[x]][1]
  dir<-sprintf("%s/%s",outlier_dir,file)
  if (dir.exists(dir)){
    dir
  } else {
    "Z"
  }
},character(1))
file<-file[file!="Z"]

cluster_pvals<-"leafcutter_outlier_clusterPvals.txt"
eff_size<-"leafcutter_outlier_effSize.txt"
pvals<-"leafcutter_outlier_pVals.txt"

outputs<-vapply(file,function(x){
  cluster_pvals_dat<-read.table(sprintf("%s/%s",x,cluster_pvals),header=TRUE)
  nme<-vapply(str_split(colnames(cluster_pvals_dat),"X"),function(x){
    x[2]
  },character(1))
  colnames(cluster_pvals_dat)<-nme
  eff_size_dat<-read.table(sprintf("%s/%s",x,eff_size),header=TRUE)
  nme<-vapply(str_split(colnames(eff_size_dat),"X"),function(x){
    x[2]
  },character(1))
  colnames(eff_size_dat)<-nme
  pvals_dat<-read.table(sprintf("%s/%s",x,pvals),header=TRUE)3
  nme<-vapply(str_split(colnames(pvals_dat),"X"),function(x){
    x[2]
  },character(1))
  colnames(pvals_dat)<-nme
  
  junctions_MD<-vapply(row.names((pvals_dat)),function(x){
    split_junc<-str_split(x,":")
    junc<-sprintf("%s:%s-%s",split_junc[[1]][1],split_junc[[1]][2],split_junc[[1]][3])
  },character(1))
  junctions_outlier<-out_calls[5:(nrow(out_calls)-3),1]
  
  length(which((junctions_outlier %in% junctions_MD) == TRUE))},numeric(1))

cluster_pvals_dat<-read.table(sprintf("%s/%s",file[1],cluster_pvals),header=TRUE)
nme<-vapply(str_split(colnames(cluster_pvals_dat),"X"),function(x){
  x[2]
},character(1))
colnames(cluster_pvals_dat)<-nme
eff_size_dat<-read.table(sprintf("%s/%s",file[1],eff_size),header=TRUE)
nme<-vapply(str_split(colnames(eff_size_dat),"X"),function(x){
  x[2]
},character(1))
colnames(eff_size_dat)<-nme
pvals_dat<-read.table(sprintf("%s/%s",file[1],pvals),header=TRUE)3
nme<-vapply(str_split(colnames(pvals_dat),"X"),function(x){
  x[2]
},character(1))
colnames(pvals_dat)<-nme

junctions_MD<-vapply(row.names((pvals_dat)),function(x){
  split_junc<-str_split(x,":")
  junc<-sprintf("%s:%s-%s",split_junc[[1]][1],split_junc[[1]][2],split_junc[[1]][3])
},character(1))
junctions_outlier<-out_calls[5:(nrow(out_calls)-3),1]

junctions_outlier %in% junctions_MD

