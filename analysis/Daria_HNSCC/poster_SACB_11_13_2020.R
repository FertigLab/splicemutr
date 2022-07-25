# this script is to explore data for the presentation taking place 10/29/2020

library(stringr)
library(stringi)
library(data.table)
library(Repitools)
library(optparse)
library(ggplot2)
library(ggpubr)

#------------------------------------------------------------------------------------------------------------------------------------------------#
countKmers <- function(seqs_vec, K, step = 1){
  # seqs_vec: vector of sequences (i.e. data_canon[,2])
  # K: kmer value
  # step: step size for kmer counting - defaults to 1
  
  # returns: a list of integer vectors with kmer 
  # character vectors as names
  # the integers represent the counts of each kmer 
  # sequence (the names, i.e. a count table)
  
  
  # uses lapply to iterate over each seq in seqs_vec
  # lapply will automatically return a list with index
  # in seqs_vec as list access index in returned list
  
  kmers <- lapply(seqs_vec, function(x){
    table(stri_sub(str = x,
                   from = seq(1, nchar(x) - K + 1,
                              by = step),
                   length = K))
  }
  )
  
  return(kmers) # return is not necessary in R, but I like to make it explicit
}

getUniqKmers <- function(kmers){
  # kmers: take kmer object returned by countKmers function
  # this is a list of integer vectors with kmer 
  # character vectors as names
  # the integers represent the counts of each kmer 
  # sequence (the names, i.e. a count table)
  
  # returns: a vector of character vectors representing
  # a unique list of kmers from the union of all sequences
  
  return(unique(unlist(lapply(kmers, names))))
  
}

split_kmers<-function(kmers,save_dir){
  # This script splits the kmer files into blocks of 2000000 or less and saves them to the out directory. 
  # This helps with computation efficiency
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
    write.table(kmers_test,file=sprintf(save_dir,x),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
    return(nrow(kmers_test))
  },numeric(1))
}

freq_to_scores<-function(mhc_dir,kmers_mhc,st = 1){
  # The purpose of this function is to take the kmer frequency per peptide list and transform 
  # that to a list of IC50 scores per kmer per peptide. The input mhc_dir is the directory and file
  # that contains the list of mhc alleles and associated scores for the unique list of kmers. kmers_mhc
  # is the list of kmer frequencies per peptide.
  
  filenames<-read.csv(mhc_dir,  header = FALSE)
  num_file<-nrow(filenames)
  start<-Sys.time()
  end<-start
  for (fi in seq(st,num_file)){
    print(sprintf("%s out of %s alleles: %s",fi,num_file,end-start))
    score<-read.csv(as.character(str_replace(filenames[fi,],":","-")))
    score<-score[!(score$peptide == "peptide"),]
    score[,2]<-as.numeric(score[,2])
    score_dt<-data.table(score)
    setkey(score_dt,peptide)
    num<-length(kmers_mhc)
    kmer_mhc_dat<-vector(mode = "list", length = num)
    kmer_mhc_dat<-lapply(kmers_mhc,function(x){
      a<-data.frame(x)
      a[3]<-score_dt[as.character(names(x))]$ic50
      freqs<-unique(a[,2])
      for (j in 1:length(freqs)){
        if (freqs[j] > 1){
          temp<-a[a[,2] == freqs[j],]
          for (k in 1:(freqs[j]-1)){a<-rbind(a,temp)}
        }
      }
      b<-a[,3]
      names(b)<-a[,1]
      return(b)
    })
    
    f<-sprintf("%s.scores.rds",filenames[fi,])
    f<-str_replace_all(f,":","-")
    saveRDS(kmer_mhc_dat,file=f)
    end<-Sys.time()
  }
}

orf<-function(seq,i){
  # finds the largest open reading frame a sequence
  # inputs:
  #   seq: the DNA sequence to search Overrepresented_sequences
  #   i: the frame to use
  # output:
  #   out: the largest orf, 0 if there is a stop before the first start, 1 if there is no start or no stop
  
  out<-0 # no coding sequence
  if ( all(unique(strsplit(seq,"")[[1]]) %in% c("A","T","C","G")) ) {
    seq<-DNAStringSet(substr(seq,i,nchar(seq)-((nchar(seq)-i+1) %% 3)))
    trans<-translate(seq,no.init.codon = FALSE)
  } else {
    trans<-seq
  }
  a<-strsplit(as.character(trans),"")
  stop<-which((a[[1]] %in% "*") == TRUE)
  start<-which((a[[1]] %in% "M") == TRUE)
  if (length(start)>0 & length(stop)>0) {
    if (stop[1] < start[1]){ # stop before first start
      out<-1
      return(out)
    } else if (length(stop)>1){ # more than one stop codon
      out<-2
    } else if (stop[1] != length(a[[1]])){ # stop not at end
      out<-3
    } else if (start[1] != 1){ # start not at beginning
      out<-4
    } else {
      out<-substr(as.character(trans),start[1],stop[1])
      return(out)
    }
  }
  return(out)
}

#------------------------------------------------------------------------------------------------------------------------------------------------#
# looking at specific target genes

data_dir<-"/media/theron/My_Passport/head_and_neck_DARIA/data/leafcutter_10_21_2020/data.Rdata"
load(data_dir)

bp<- ggplot(intron_summary, aes(x="", y=n, fill=Results))+
  geom_bar(width = 1, stat = "identity")
pie <- bp + coord_polar("y", start=0) +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5)) + 
  scale_fill_brewer(palette="Greens")
pie

genes_of_interest<-c("PIK3R1","NPHP3","TTC39A","CD44","PLS3")

out_dir<-"/media/theron/My_Passport/head_and_neck_DARIA/data/leafcutter_10_21_2020/kmers"

# loading in the peps files
introns_dir<-"/media/theron/My_Passport/head_and_neck_DARIA/data/leafcutter_10_21_2020/introns"
iter<-1
for (i in seq(30)){
  file<-sprintf("%s/peptides%s.rds",introns_dir,i)
  if (file.exists(file)){
    if (iter == 1){
      peps<-readRDS(file)
    } else {
      peps_fill<-readRDS(file)
      peps<-rbind(peps,peps_fill)
    }
    iter<-iter+1
  }
}
peps<-unique(peps)

peps_targets<-peps[peps$gene %in% genes_of_interest,]
peps_targets<-unique(peps_targets)
peps_full<-peps[!(peps_targets$peptide == 1 | peps_targets$peptide == 0),]
peps_no_stop<-sapply(peps_full$peptide,function(x){return(substr(x,1,nchar(x)-1))},USE.NAMES = FALSE)
peps_no_stop_frame<-data.frame(as.character(peps_no_stop))
colnames(peps_no_stop_frame)<-"peps"

for (i in 14:25){
  print(i)
  kmer_length=i
  peps_no_stop<-sapply(peps_no_stop,function(x){
    Z<-c()
    len<-nchar(x)-kmer_length
    if(len < 0){Z[1:abs(len)]<-"Z"}
    paste(x,paste(Z,collapse=""),sep="")
  },USE.NAMES = FALSE)
  
  kmers<-countKmers(peps_no_stop,kmer_length,1)
  kmers_unique<-getUniqKmers(kmers)
  kmers_frame<-data.frame(kmers_unique)
  
  out<-sprintf("%s/%s",out_dir,"kmers_%s.txt")
  out_rds<-sprintf("%s/kmers_%s.rds",out_dir,kmer_length)
  write.table(kmers_unique,file=sprintf(out,kmer_length),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
  saveRDS(kmers,file=out_rds)
}

#------------------------------------------------------------------------------------------------------------------------------------------------#
# looking at all genes

data_dir<-"/media/theron/My_Passport/head_and_neck_DARIA/data/leafcutter_10_21_2020/data.Rdata"
load(data_dir)

out_dir<-"/media/theron/My_Passport/head_and_neck_DARIA/data/leafcutter_10_21_2020/kmers"

# loading in the peps files
introns_dir<-"/media/theron/My_Passport/head_and_neck_DARIA/data/leafcutter_10_21_2020/introns"
iter<-1
for (i in seq(30)){
  file<-sprintf("%s/peptides%s.rds",introns_dir,i)
  if (file.exists(file)){
    if (iter == 1){
      peps<-readRDS(file)
    } else {
      peps_fill<-readRDS(file)
      peps<-rbind(peps,peps_fill)
    }
    iter<-iter+1
  }
}
peps<-unique(peps)

peps_full<-peps[!(peps$peptide == 1 | peps$peptide == 0),]
peps_no_stop<-sapply(peps_full$peptide,function(x){return(substr(x,1,nchar(x)-1))},USE.NAMES = FALSE)
peps_no_stop_frame<-data.frame(as.character(peps_no_stop))
colnames(peps_no_stop_frame)<-"peps"

for (i in 22:25){
  print(i)
  kmer_length=i
  peps_no_stop<-sapply(peps_no_stop,function(x){
    Z<-c()
    len<-nchar(x)-kmer_length
    if(len < 0){Z[1:abs(len)]<-"Z"}
    paste(x,paste(Z,collapse=""),sep="")
  },USE.NAMES = FALSE)
  
  kmers<-countKmers(peps_no_stop,kmer_length,1)
  kmers_unique<-getUniqKmers(kmers)
  kmers_frame<-data.frame(kmers_unique)
  
  out<-sprintf("%s/%s",out_dir,"kmers_%s.txt")
  out_rds<-sprintf("%s/kmers_%s.rds",out_dir,kmer_length)
  write.table(kmers_unique,file=sprintf(out,kmer_length),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
  saveRDS(kmers,file=out_rds)
}

#------------------------------------------------------------------------------------------------------------------------------------------------#

peps$peptide<-NA
for (i in 1:nrow(peps)) {
  print(sprintf("%d rows complete out of %d total",i,nrow(peps)))
  peps[i,"peptide"]<-orf(peps$sequence[i],peps$frame[i])
}
saveRDS(peps,file="/media/theron/My_Passport/head_and_neck_DARIA/data/leafcutter_10_21_2020/introns/peps.rds")

#------------------------------------------------------------------------------------------------------------------------------------------------#
# creating metadata plots for all data

# From first frame, translatable based on junction type
codes<-c("0","1","2","3","4")
peps_small<-peps[peps$frame == 1,]
peps_small_2<-peps_small[!(peps_small$peptide %in% codes),]
translatable<-table(peps_small_2$verdict)

peps_small_3<-peps_small[(peps_small$peptide == "0"),]
no_code_seq<-table(peps_small_3$verdict)

peps_small_4<-peps_small[(peps_small$peptide == "1"),]
stop_before_start<-table(peps_small_4$verdict)

peps_small_5<-peps_small[(peps_small$peptide == "2"),]
multiple_stop<-table(peps_small_5$verdict)

peps_small_6<-peps_small[(peps_small$peptide == "3"),]
stop_not_end<-table(peps_small_6$verdict)

peps_small_7<-peps_small[(peps_small$peptide == "4"),]
start_not_beginning<-table(peps_small_7$verdict)

total<-table(peps_small$verdict)

# creating grouped bar plots 
type_num<-5
condition_num<-6
errors<-c(rep("translatable",type_num), rep("no coding sequence",type_num),
          rep("stop before start",type_num), rep("multiple stop",type_num), 
          rep("sequence aft. stop",type_num), rep("start not at beginning",type_num))
vals<-data.frame()
vals<-rbind(vals,translatable)
vals<-rbind(vals,no_code_seq)
vals<-rbind(vals,stop_before_start)
vals<-rbind(vals,multiple_stop)
vals<-rbind(vals,stop_not_end)
vals<-rbind(vals,start_not_beginning)
cond<-c("annotated","cryptic fiveprime","cryptic threeprime","cryptic unanchored","novel annotated pair")
colnames(vals)<-cond
tots<-apply(vals,2,sum)

condition<-rep(c(sprintf("%s (%s)","annotated",tots[1]), sprintf("%s (%s)", "cryptic fiveprime", tots[2]),
                 sprintf("%s (%s)","cryptic threeprime",tots[3]), 
                 sprintf("%s (%s)","cryptic unanchored",tots[4]),sprintf("%s (%s)","novel annotated pair",tots[5])),condition_num)
values<-c((translatable/tots)*100, (no_code_seq/tots)*100, (stop_before_start/tots)*100, 
          (multiple_stop/tots)*100, 
          (stop_not_end/tots)*100, 
          (start_not_beginning/tots)*100)
vals<-c(translatable, no_code_seq, stop_before_start, multiple_stop,stop_not_end, start_not_beginning)
data<-data.frame(errors,condition,values,vals)
# Grouped
ggplot(data, aes(fill=errors, y=values, x=condition)) + 
  geom_bar(position="dodge", stat="identity")

# creating grouped bar plots just total
type_num<-5
condition_num<-1
errors<-c(rep("total",type_num))
condition<-rep(c("annotated","cryptic fiveprime","cryptic threeprime","cryptic unanchored","novel annotated pair"),condition_num)
values<-c(total)
data<-data.frame(errors,condition,values)
# Grouped
ggplot(data, aes(fill=errors, y=values, x=condition)) + 
  geom_bar(position="dodge", stat="identity", fill="black", colour="black") +
ylim(0,45000)

#------------------------------------------------------------------------------------------------------------------------------------------------#
  
# psi analysis > 0
  codes<-c("0","1","2","3","4")
peps_small<-peps[peps$frame == 1,]
peps_small<-peps_small[peps_small$deltapsi > 0,]
peps_small_2<-peps_small[!(peps_small$peptide %in% codes),]
translatable<-table(peps_small_2$verdict)

peps_small_3<-peps_small[(peps_small$peptide == "0"),]
no_code_seq<-table(peps_small_3$verdict)

peps_small_4<-peps_small[(peps_small$peptide == "1"),]
stop_before_start<-table(peps_small_4$verdict)

peps_small_5<-peps_small[(peps_small$peptide == "2"),]
multiple_stop<-table(peps_small_5$verdict)

peps_small_6<-peps_small[(peps_small$peptide == "3"),]
stop_not_end<-table(peps_small_6$verdict)
stop_not_end<-c(stop_not_end,stop_not_end[4])
stop_not_end[4]<-0

peps_small_7<-peps_small[(peps_small$peptide == "4"),]
start_not_beginning<-table(peps_small_7$verdict)
start_not_beginning<-c(start_not_beginning,start_not_beginning[4])
stop_not_end[4]<-0

total<-table(peps_small$verdict)

# creating grouped bar plots 
type_num<-5
condition_num<-6
errors<-c(rep("translatable",type_num), rep("no coding sequence",type_num),
          rep("stop before start",type_num), rep("multiple stop",type_num), 
          rep("sequence aft. stop",type_num), rep("start not at beginning",type_num))
vals<-data.frame()
vals<-rbind(vals,translatable)
vals<-rbind(vals,no_code_seq)
vals<-rbind(vals,stop_before_start)
vals<-rbind(vals,multiple_stop)
vals<-rbind(vals,stop_not_end)
vals<-rbind(vals,start_not_beginning)
cond<-c("annotated","cryptic fiveprime","cryptic threeprime","cryptic unanchored","novel annotated pair")
colnames(vals)<-cond
tots<-apply(vals,2,sum)

condition<-rep(c(sprintf("%s (%s)","annotated",tots[1]), sprintf("%s (%s)", "cryptic fiveprime", tots[2]),
                 sprintf("%s (%s)","cryptic threeprime",tots[3]), 
                 sprintf("%s (%s)","cryptic unanchored",tots[4]),sprintf("%s (%s)","novel annotated pair",tots[5])),condition_num)
values<-c((translatable/tots)*100, (no_code_seq/tots)*100, (stop_before_start/tots)*100, 
          (multiple_stop/tots)*100, 
          (stop_not_end/tots)*100, 
          (start_not_beginning/tots)*100)
vals<-c(translatable, no_code_seq, stop_before_start, multiple_stop,stop_not_end, start_not_beginning)
data<-data.frame(errors,condition,values,vals)
# Grouped
ggplot(data, aes(fill=errors, y=values, x=condition)) + 
  geom_bar(position="dodge", stat="identity")

# creating grouped bar plots just total
type_num<-5
condition_num<-1
errors<-c(rep("total",type_num))
condition<-rep(c("annotated","cryptic fiveprime","cryptic threeprime","cryptic unanchored","novel annotated pair"),condition_num)
values<-c(total)
data<-data.frame(errors,condition,values)
# Grouped
ggplot(data, aes(fill=errors, y=values, x=condition)) + 
  geom_bar(position="dodge", stat="identity", fill= "black",colour="black") +
  ylim(0,25000)


# psi analysis < 0
codes<-c("0","1","2","3","4")
peps_small<-peps[peps$frame == 1,]
peps_small<-peps_small[peps_small$deltapsi < 0,]
peps_small_2<-peps_small[!(peps_small$peptide %in% codes),]
translatable<-table(peps_small_2$verdict)

peps_small_3<-peps_small[(peps_small$peptide == "0"),]
no_code_seq<-table(peps_small_3$verdict)

peps_small_4<-peps_small[(peps_small$peptide == "1"),]
stop_before_start<-table(peps_small_4$verdict)

peps_small_5<-peps_small[(peps_small$peptide == "2"),]
multiple_stop<-table(peps_small_5$verdict)

peps_small_6<-peps_small[(peps_small$peptide == "3"),]
stop_not_end<-table(peps_small_6$verdict)

peps_small_7<-peps_small[(peps_small$peptide == "4"),]
start_not_beginning<-table(peps_small_7$verdict)

total<-table(peps_small$verdict)

# creating grouped bar plots 
type_num<-5
condition_num<-6
errors<-c(rep("translatable",type_num), rep("no coding sequence",type_num),
          rep("stop before start",type_num), rep("multiple stop",type_num), 
          rep("sequence aft. stop",type_num), rep("start not at beginning",type_num))
vals<-data.frame()
vals<-rbind(vals,translatable)
vals<-rbind(vals,no_code_seq)
vals<-rbind(vals,stop_before_start)
vals<-rbind(vals,multiple_stop)
vals<-rbind(vals,stop_not_end)
vals<-rbind(vals,start_not_beginning)
cond<-c("annotated","cryptic fiveprime","cryptic threeprime","cryptic unanchored","novel annotated pair")
colnames(vals)<-cond
tots<-apply(vals,2,sum)

condition<-rep(c(sprintf("%s (%s)","annotated",tots[1]), sprintf("%s (%s)", "cryptic fiveprime", tots[2]),
                 sprintf("%s (%s)","cryptic threeprime",tots[3]), 
                 sprintf("%s (%s)","cryptic unanchored",tots[4]),sprintf("%s (%s)","novel annotated pair",tots[5])),condition_num)
values<-c((translatable/tots)*100, (no_code_seq/tots)*100, (stop_before_start/tots)*100, 
          (multiple_stop/tots)*100, 
          (stop_not_end/tots)*100, 
          (start_not_beginning/tots)*100)
vals<-c(translatable, no_code_seq, stop_before_start, multiple_stop,stop_not_end, start_not_beginning)
data<-data.frame(errors,condition,values,vals)
# Grouped
ggplot(data, aes(fill=errors, y=values, x=condition)) + 
  geom_bar(position="dodge", stat="identity")

# creating grouped bar plots just total
type_num<-5
condition_num<-1
errors<-c(rep("total",type_num))
condition<-rep(c("annotated","cryptic fiveprime","cryptic threeprime","cryptic unanchored","novel annotated pair"),condition_num)
values<-c(total)
data<-data.frame(errors,condition,values)
# Grouped
ggplot(data, aes(fill=errors, y=values, x=condition)) + 
  geom_bar(position="dodge", stat="identity", fill= "black",colour="black") +
  ylim(0,25000)

#------------------------------------------------------------------------------------------------------------------------------------------------#
#creating a divergent plot

# psi analysis > 0
codes<-c("0","1","2","3","4")
peps_small<-peps[peps$frame == 1,]
peps_small<-peps_small[peps_small$deltapsi > 0,]
peps_small_2<-peps_small[!(peps_small$peptide %in% codes),]
translatable_1<-table(peps_small_2$verdict)

peps_small_3<-peps_small[(peps_small$peptide == "0"),]
no_code_seq_1<-table(peps_small_3$verdict)

peps_small_4<-peps_small[(peps_small$peptide == "1"),]
stop_before_start_1<-table(peps_small_4$verdict)

peps_small_5<-peps_small[(peps_small$peptide == "2"),]
multiple_stop_1<-table(peps_small_5$verdict)

peps_small_6<-peps_small[(peps_small$peptide == "3"),]
stop_not_end_1<-table(peps_small_6$verdict)
stop_not_end_1<-c(stop_not_end_1,stop_not_end[4])
stop_not_end_1[4]<-0

peps_small_7<-peps_small[(peps_small$peptide == "4"),]
start_not_beginning_1<-table(peps_small_7$verdict)
start_not_beginning_1<-c(start_not_beginning_1,start_not_beginning[4])
stop_not_end_1[4]<-0

total_1<-table(peps_small$verdict)

# psi analysis < 0
codes<-c("0","1","2","3","4")
peps_small<-peps[peps$frame == 1,]
peps_small<-peps_small[peps_small$deltapsi < 0,]
peps_small_2<-peps_small[!(peps_small$peptide %in% codes),]
translatable_2<-table(peps_small_2$verdict)

peps_small_3<-peps_small[(peps_small$peptide == "0"),]
no_code_seq_2<-table(peps_small_3$verdict)

peps_small_4<-peps_small[(peps_small$peptide == "1"),]
stop_before_start_2<-table(peps_small_4$verdict)

peps_small_5<-peps_small[(peps_small$peptide == "2"),]
multiple_stop_2<-table(peps_small_5$verdict)

peps_small_6<-peps_small[(peps_small$peptide == "3"),]
stop_not_end_2<-table(peps_small_6$verdict)

peps_small_7<-peps_small[(peps_small$peptide == "4"),]
start_not_beginning_2<-table(peps_small_7$verdict)

total_2<-table(peps_small$verdict)

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
condition_num<-6
errors<-c(rep("translatable",type_num), rep("no coding sequence",type_num),
          rep("stop before start",type_num), rep("multiple stop",type_num), 
          rep("sequence aft. stop",type_num), rep("start not at beginning",type_num))
vals<-data.frame()
vals<-rbind(vals,translatable)
vals<-rbind(vals,no_code_seq)
vals<-rbind(vals,stop_before_start)
vals<-rbind(vals,multiple_stop)
vals<-rbind(vals,stop_not_end)
vals<-rbind(vals,start_not_beginning)
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
          (start_not_beginning/tots)*100)
vals<-c(translatable, no_code_seq, stop_before_start, multiple_stop,stop_not_end, start_not_beginning)
data<-data.frame(errors,condition,values,vals)
# Grouped
ggplot(data, aes(fill=errors, y=values, x=condition)) + 
  geom_bar(position="dodge", stat="identity")


# creating grouped bar plots just total
type_num<-5
condition_num<-1
errors<-c(rep("total",type_num))
condition<-rep(c("annotated","cryptic fiveprime","cryptic threeprime","cryptic unanchored","novel annotated pair"),condition_num)
values<-c(total)
data<-data.frame(errors,condition,values)
# Grouped
ggplot(data, aes(fill=condition, y=values, x=errors)) + 
  geom_bar(position="dodge", stat="identity") 

#------------------------------------------------------------------------------------------------------------------------------------------------#

# looking at CD44
peps_small<-peps[peps$frame == 1,]
peps_small<-peps_small[peps_small$gene== "CD44",]
peps_small_2<-peps_small[!(peps_small$peptide %in% codes),]
table(peps_small_2$verdict)

peps_small_3<-peps_small[(peps_small$peptide == "0"),]
table(peps_small_3$verdict)

peps_small_4<-peps_small[(peps_small$peptide == "1"),]
table(peps_small_4$verdict)

peps_small_5<-peps_small[(peps_small$peptide == "2"),]
table(peps_small_5$verdict)

peps_small_6<-peps_small[(peps_small$peptide == "3"),]
table(peps_small_6$verdict)

peps_small_7<-peps_small[(peps_small$peptide == "4"),]
table(peps_small_7$verdict)

table(peps_small$verdict)

#NPHP3
peps_small<-peps[peps$frame == 1,]
peps_small<-peps_small[peps_small$gene== "NPHP3",]
peps_small_2<-peps_small[!(peps_small$peptide %in% codes),]
table(peps_small_2$verdict)

peps_small_3<-peps_small[(peps_small$peptide == "0"),]
table(peps_small_3$verdict)

peps_small_4<-peps_small[(peps_small$peptide == "1"),]
table(peps_small_4$verdict)

peps_small_5<-peps_small[(peps_small$peptide == "2"),]
table(peps_small_5$verdict)

peps_small_6<-peps_small[(peps_small$peptide == "3"),]
table(peps_small_6$verdict)

peps_small_7<-peps_small[(peps_small$peptide == "4"),]
table(peps_small_7$verdict)

table(peps_small$verdict)

#------------------------------------------------------------------------------------------------------------------------------------------------#
# annotating peps with deltapsi value

#------------------------------------------------------------------------------------------------------------------------------------------------#

filenames<-"/media/theron/My_Passport/head_and_neck_DARIA/data/leafcutter_10_21_2020/kmers/predictions/filenames.txt"
kmers<-readRDS("/media/theron/My_Passport/head_and_neck_DARIA/data/leafcutter_10_21_2020/kmers/kmers_9.rds")
freq_to_scores(filenames,kmers,1)

#------------------------------------------------------------------------------------------------------------------------------------------------#
# looking at translateable reference cds transcripts
ref_peps<-readRDS("/media/theron/My_Passport/head_and_neck_DARIA/data/leafcutter_10_21_2020/ref/peptides/peptides.rds")
ref_peps_tx<-vapply(ref_peps$peptide,function(x){
  if (nchar(x)>1){
    return("tx")
  } else {
    return(x)
  }
},character(1))
cds_breakdown<-table(ref_peps_tx)
cds_breakdown<-data.frame(cds_breakdown)
colnames(cds_breakdown)<-c("type","count")
percent<-vapply(cds_breakdown$count, function(x){
  (x/sum(cds_breakdown$count))*100
},numeric(1))
cds_breakdown$percent<-percent

bp<- ggplot(cds_breakdown, aes(x="", y=percent, fill=type))+
  geom_bar(width = 1, stat = "identity")
pie <- bp + coord_polar("y", start=0) +
  geom_text(aes(label = type), position = position_stack(vjust = 0.5))
pie

target_peps<-readRDS("/media/theron/My_Passport/head_and_neck_DARIA/data/leafcutter_10_21_2020/introns/peps.rds")
target_peps<-target_peps[target_peps$frame == "1",]
test_1<-all(unique(target_peps$tx_id) %in% ref_peps$tx)
test_2<-which(unique(target_peps$tx_id) %in% ref_peps$tx[which(ref_peps_tx == "tx")])
test_3<-which(target_peps$tx_id %in% ref_peps$tx[which(ref_peps_tx == "tx")])

stat_1<-length(test_3)/nrow(target_peps) * 100

translatable_target_peps<-target_peps[target_peps$tx_id %in% ref_peps$tx[which(ref_peps_tx == "tx")],]

#------------------------------------------------------------------------------------------------------------------------------------------------#
# From first frame, translatable based on junction type

codes<-c("0","1","2","3","4","5","6")
peps_small<-translatable_target_peps
peps_small_2<-peps_small[!(peps_small$peptide %in% codes),]
translatable<-table(peps_small_2$verdict)

peps_small_3<-peps_small[(peps_small$peptide == "0"),]
no_code_seq<-table(peps_small_3$verdict)

peps_small_4<-peps_small[(peps_small$peptide == "1"),]
stop_before_start<-table(peps_small_4$verdict)

peps_small_5<-peps_small[(peps_small$peptide == "2"),]
multiple_stop<-table(peps_small_5$verdict)

peps_small_6<-peps_small[(peps_small$peptide == "3"),]
stop_not_end<-table(peps_small_6$verdict)
stop_not_end<-c(stop_not_end,stop_not_end[4])
stop_not_end[4]<-0

peps_small_7<-peps_small[(peps_small$peptide == "4"),]
start_not_beginning<-table(peps_small_7$verdict)


total<-table(peps_small$verdict)
# creating grouped bar plots 
type_num<-5
condition_num<-6
errors<-c(rep("translatable",type_num), rep("no coding sequence",type_num),
          rep("stop before start",type_num), rep("multiple stop",type_num), 
          rep("sequence aft. stop",type_num), rep("start not at beginning",type_num))
vals<-data.frame()
vals<-rbind(vals,translatable)
vals<-rbind(vals,no_code_seq)
vals<-rbind(vals,stop_before_start)
vals<-rbind(vals,multiple_stop)
vals<-rbind(vals,stop_not_end)
vals<-rbind(vals,start_not_beginning)

cond<-c("annotated","cryptic fiveprime","cryptic threeprime","cryptic unanchored","novel annotated pair")
colnames(vals)<-cond
tots<-apply(vals,2,sum)

condition<-rep(c(sprintf("%s (%s)","annotated",tots[1]), sprintf("%s (%s)", "cryptic fiveprime", tots[2]),
                 sprintf("%s (%s)","cryptic threeprime",tots[3]), 
                 sprintf("%s (%s)","cryptic unanchored",tots[4]),sprintf("%s (%s)","novel annotated pair",tots[5])),condition_num)
values<-c((translatable/tots)*100, (no_code_seq/tots)*100, (stop_before_start/tots)*100, 
          (multiple_stop/tots)*100, 
          (stop_not_end/tots)*100, 
          (start_not_beginning/tots)*100)
vals<-c(translatable, no_code_seq, stop_before_start, multiple_stop,stop_not_end, 
        start_not_beginning)
data<-data.frame(errors,condition,values,vals)
# Grouped
ggplot(data, aes(fill=errors, y=values, x=condition)) + 
  geom_bar(position="dodge", stat="identity") + ylab("Percentage of Total") + 
  xlab("Leafcutter Verdict") + theme(text = element_text(size=14),axis.text = element_text(size=14))

# creating grouped bar plots just total
type_num<-5
condition_num<-1
errors<-c(rep("total",type_num))
condition<-rep(c("annotated","cryptic fiveprime","cryptic threeprime","cryptic unanchored","novel annotated pair"),condition_num)
values<-c(total)
data<-data.frame(errors,condition,values)
# Grouped
ggplot(data, aes(fill=errors, y=values, x=condition)) + 
  geom_bar(position="dodge", stat="identity", fill="black", colour="black") +
  ylim(0,45000) + theme(text = element_text(size=16),axis.text = element_text(size=16))

#------------------------------------------------------------------------------------------------------------------------------------------------#
#creating a divergent plot

# psi analysis > 0
codes<-c("0","1","2","3","4","5","6")
peps_small<-translatable_target_peps
peps_small<-peps_small[peps_small$deltapsi > 0,]
peps_small_2<-peps_small[!(peps_small$peptide %in% codes),]
translatable_1<-table(peps_small_2$verdict)

peps_small_3<-peps_small[(peps_small$peptide == "0"),]
no_code_seq_1<-table(peps_small_3$verdict)

peps_small_4<-peps_small[(peps_small$peptide == "1"),]
stop_before_start_1<-table(peps_small_4$verdict)

peps_small_5<-peps_small[(peps_small$peptide == "2"),]
multiple_stop_1<-table(peps_small_5$verdict)

peps_small_6<-peps_small[(peps_small$peptide == "3"),]
stop_not_end_1<-table(peps_small_6$verdict)
stop_not_end_1<-c(stop_not_end_1,stop_not_end_1[4])
stop_not_end_1[4]<-0

peps_small_7<-peps_small[(peps_small$peptide == "4"),]
start_not_beginning_1<-table(peps_small_7$verdict)
start_not_beginning_1<-c(start_not_beginning_1,start_not_beginning_1[4])
start_not_beginning_1[4]<-0

total_1<-table(peps_small$verdict)

# psi analysis < 0
codes<-c("0","1","2","3","4","5","6")
peps_small<-translatable_target_peps
peps_small<-peps_small[peps_small$deltapsi < 0,]
peps_small_2<-peps_small[!(peps_small$peptide %in% codes),]
translatable_2<-table(peps_small_2$verdict)

peps_small_3<-peps_small[(peps_small$peptide == "0"),]
no_code_seq_2<-table(peps_small_3$verdict)

peps_small_4<-peps_small[(peps_small$peptide == "1"),]
stop_before_start_2<-table(peps_small_4$verdict)

peps_small_5<-peps_small[(peps_small$peptide == "2"),]
multiple_stop_2<-table(peps_small_5$verdict)

peps_small_6<-peps_small[(peps_small$peptide == "3"),]
stop_not_end_2<-table(peps_small_6$verdict)
stop_not_end_2<-c(stop_not_end_2,stop_not_end_2[4])
stop_not_end_2[4]<-0

peps_small_7<-peps_small[(peps_small$peptide == "4"),]
start_not_beginning_2<-table(peps_small_7$verdict)

total_2<-table(peps_small$verdict)

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
condition_num<-6
errors<-c(rep("translatable",type_num), rep("no coding sequence",type_num),
          rep("stop before start",type_num), rep("multiple stop",type_num), 
          rep("sequence aft. stop",type_num), rep("start not at beginning",type_num))
vals<-data.frame()
vals<-rbind(vals,translatable)
vals<-rbind(vals,no_code_seq)
vals<-rbind(vals,stop_before_start)
vals<-rbind(vals,multiple_stop)
vals<-rbind(vals,stop_not_end)
vals<-rbind(vals,start_not_beginning)
cond<-c("annotated","cryptic fiveprime","cryptic threeprime","cryptic unanchored","novel annotated pair")
colnames(vals)<-cond
tots_diff<-apply(vals,2,function(x){
  sum(abs(x))
})

condition<-rep(c(sprintf("%s (%s)","annotated",tots[1]), sprintf("%s (%s)", "cryptic fiveprime", tots[2]),
                 sprintf("%s (%s)","cryptic threeprime",tots[3]), 
                 sprintf("%s (%s)","cryptic unanchored",tots[4]),sprintf("%s (%s)","novel annotated pair",tots[5])),condition_num)
values<-c((translatable/tots)*100, (no_code_seq/tots)*100, (stop_before_start/tots)*100, 
          (multiple_stop/tots)*100, 
          (stop_not_end/tots)*100, 
          (start_not_beginning/tots)*100)
vals<-c(translatable, no_code_seq, stop_before_start, multiple_stop,stop_not_end, start_not_beginning)
data<-data.frame(errors,condition,values,vals)
# Grouped
ggplot(data, aes(fill=errors, y=values, x=condition)) + 
  geom_bar(position="dodge", stat="identity") + ylab("Percentage of Total") + 
  xlab("Leafcutter Verdict")

#------------------------------------------------------------------------------------------------------------------------------------------------#
# working with metric data

HLA_A02_01<-readRDS("/media/theron/My_Passport/head_and_neck_DARIA/data/leafcutter_10_21_2020/kmers/processed/HLA-A02-01_kmers_9_metrics_all.rds")
HLA_A02_01<-HLA_A02_01[!is.na(HLA_A02_01$KS),]
HLA_B07_02<-readRDS("/media/theron/My_Passport/head_and_neck_DARIA/data/leafcutter_10_21_2020/kmers/processed/HLA-B07-02_kmers_9_metrics_all.rds")
HLA_B07_02<-HLA_B07_02[!is.na(HLA_B07_02$KS),]
HLA_C07_01<-readRDS("/media/theron/My_Passport/head_and_neck_DARIA/data/leafcutter_10_21_2020/kmers/processed/HLA-C07-01_kmers_9_metrics_all.rds")
HLA_C07_01<-HLA_C07_01[!is.na(HLA_C07_01$KS),]

HLA_A02_01_test_1<-any(HLA_A02_01$SB_mod > 0)
HLA_A02_01_test_2<-any(HLA_A02_01$WB_mod > 0)
HLA_A02_01_test_3<-any(HLA_A02_01$SB_mod > 0)
HLA_A02_01_test_4<-any(HLA_A02_01$WB_mod > 0)

HLA_B07_02_test_1<-any(HLA_B07_02$SB_mod > 0)
HLA_B07_02_test_2<-any(HLA_B07_02$WB_mod > 0)
HLA_B07_02_test_3<-any(HLA_B07_02$SB_mod > 0)
HLA_B07_02_test_4<-any(HLA_B07_02$WB_mod > 0)

HLA_C07_01_test_1<-any(HLA_C07_01$SB_mod > 0)
HLA_C07_01_test_2<-any(HLA_C07_01$WB_mod > 0)
HLA_C07_01_test_3<-any(HLA_C07_01$SB_mod > 0)
HLA_C07_01_test_4<-any(HLA_C07_01$WB_mod > 0)


ggplot(ToothGrowth, aes(x=dose, y=len)) + 
  geom_violin() + geom_boxplot(width=0.1)

#------------------------------------------------------------------------------------------------------------------------------------------------#
# calculating peptide lengths

pep_len<-vapply(HLA_A02_01$peptide,function(x){nchar(x)},numeric(1))
HLA_A02_01$pep_len<-pep_len
HLA_A02_01$type<-"normal"
HLA_A02_01[HLA_A02_01$deltapsi>0,"type"]<-"tumor"

pep_len<-vapply(HLA_B07_02$peptide,function(x){nchar(x)},numeric(1))
HLA_B07_02$pep_len<-pep_len
HLA_B07_02$type<-"normal"
HLA_B07_02[HLA_B07_02$deltapsi>0,"type"]<-"tumor"

pep_len<-vapply(HLA_C07_01$peptide,function(x){nchar(x)},numeric(1))
HLA_C07_01$pep_len<-pep_len
HLA_C07_01$type<-"normal"
HLA_C07_01[HLA_C07_01$deltapsi>0,"type"]<-"tumor"

#------------------------------------------------------------------------------------------------------------------------------------------------#
# generating violin plots
#
ggplot(HLA_A02_01, aes(x=type, y=KS)) + 
  geom_violin() + geom_boxplot(width=0.1)

ggplot(HLA_A02_01, aes(x=type, y=perc)) + 
  geom_violin() + geom_boxplot(width=0.1)

ggplot(HLA_A02_01, aes(x=type, y=SB_mod)) + 
  geom_violin() + geom_boxplot(width=0.1)

ggplot(HLA_A02_01, aes(x=type, y=WB_mod)) + 
  geom_violin() + geom_boxplot(width=0.1)

ggplot(HLA_A02_01, aes(x=type, y=pep_len)) + 
  geom_violin() + geom_boxplot(width=0.1)


# 
ggplot(HLA_B07_02, aes(x=type, y=KS)) + 
  geom_violin() + geom_boxplot(width=0.1)

ggplot(HLA_B07_02, aes(x=type, y=perc)) + 
  geom_violin() + geom_boxplot(width=0.1)

ggplot(HLA_B07_02, aes(x=type, y=SB_mod)) + 
  geom_violin() + geom_boxplot(width=0.1)

ggplot(HLA_B07_02, aes(x=type, y=WB_mod)) + 
  geom_violin() + geom_boxplot(width=0.1)

ggplot(HLA_B07_02, aes(x=type, y=pep_len)) + 
  geom_violin() + geom_boxplot(width=0.1)

# 
ggplot(HLA_C07_01, aes(x=type, y=KS)) + 
  geom_violin() + geom_boxplot(width=0.1)

ggplot(HLA_C07_01, aes(x=type, y=perc)) + 
  geom_violin() + geom_boxplot(width=0.1)

ggplot(HLA_C07_01, aes(x=type, y=SB_mod)) + 
  geom_violin() + geom_boxplot(width=0.1)

ggplot(HLA_C07_01, aes(x=type, y=WB_mod)) + 
  geom_violin() + geom_boxplot(width=0.1)

ggplot(HLA_C07_01, aes(x=type, y=pep_len)) + 
  geom_violin() + geom_boxplot(width=0.1)

# changed and add_reg vs normal cds
changed<-(HLA_A02_01$modified == "changed" | HLA_A02_01$modified == "add_reg")
HLA_A02_01_violin<-ggplot(HLA_A02_01[changed,], aes(x=type, y=KS)) + 
  geom_violin() + geom_boxplot(width=0.1)
HLA_A02_01_tumor_out<-boxplot.stats(HLA_A02_01[HLA_A02_01$type[changed]=="tumor","KS"])$out
HLA_A02_01_normal_out<-boxplot.stats(HLA_A02_01[HLA_A02_01$type[changed]=="normal","KS"])$out
ggarrange(a, b, c, d, 
          labels = c("A", "B", "C","D"),
          ncol = 2, nrow = 2)

ggplot(HLA_A02_01[changed,], aes(x=type, y=perc)) + 
  geom_violin() + geom_boxplot(width=0.1)

ggplot(HLA_A02_01[changed,], aes(x=type, y=SB_mod)) + 
  geom_violin() + geom_boxplot(width=0.1)

ggplot(HLA_A02_01[changed,], aes(x=type, y=WB_mod)) + 
  geom_violin() + geom_boxplot(width=0.1)

ggplot(HLA_A02_01[changed,], aes(x=type, y=pep_len)) + 
  geom_violin() + geom_boxplot(width=0.1)


# 
changed<-(HLA_B07_02$modified == "changed" | HLA_B07_02$modified == "add_reg")
HLA_B07_02_violin<-ggplot(HLA_B07_02[changed,], aes(x=type, y=KS)) + 
  geom_violin() + geom_boxplot(width=0.1)
HLA_B07_02_violin

HLA_B07_02_tumor_out<-boxplot.stats(HLA_B07_02[HLA_B07_02$type[changed]=="tumor","KS"])$out
tumor_out<-HLA_B07_02[HLA_B07_02$type == "tumor","KS"] %in% HLA_B07_02_tumor_out
HLA_B07_02_normal_out<-boxplot.stats(HLA_B07_02[HLA_B07_02$type[changed]=="normal","KS"])$out
normal_out<-HLA_B07_02[HLA_B07_02$type == "normal","KS"] %in% HLA_B07_02_normal_out

HLA_B07_02_cong<-rbind(HLA_B07_02[normal_out,],HLA_B07_02[tumor_out,])
a<-ggplot(HLA_B07_02_cong[HLA_B07_02_cong$type=="normal",], aes(x = SB_mod, y = KS)) + geom_point() + labs(y="KS") + labs(x="Strong Binder Count") + labs(title="Normal Associated") + geom_smooth(method="gam",formula = y~s(x,bs="cs")) + 
  theme(text = element_text(size=20),axis.text = element_text(size=20))
b<-ggplot(HLA_B07_02_cong[HLA_B07_02_cong$type=="normal",], aes(x = WB_mod, y = KS)) + geom_point() + labs(y="KS") + labs(x="Weak Binder Count") + labs(title="Normal Associated") + geom_smooth(method="gam",formula = y~s(x,bs="cr")) +
  theme(text = element_text(size=20),axis.text = element_text(size=20))
c<-ggplot(HLA_B07_02_cong[HLA_B07_02_cong$type=="tumor",], aes(x = SB_mod, y = KS)) + geom_point() + labs(y="KS") + labs(x="Strong Binder Count") + labs(title="Tumor Associated") + geom_smooth(method="gam",formula = y~s(x,bs="cr")) + 
  theme(text = element_text(size=20),axis.text = element_text(size=20))
d<-ggplot(HLA_B07_02_cong[HLA_B07_02_cong$type=="tumor",], aes(x = WB_mod, y = KS)) + geom_point() + labs(y="KS") + labs(x="Weak Binder Count") + labs(title="Tumor Associated") + geom_smooth(method="gam",formula = y~s(x,bs="cr")) + 
  theme(text = element_text(size=20),axis.text = element_text(size=20))

ggarrange(a, b, c, d, 
          labels = c("A", "B", "C","D"),
          ncol = 2, nrow = 2)
  
ggplot(HLA_B07_02[changed,], aes(x=type, y=perc)) + 
  geom_violin() + geom_boxplot(width=0.1)

ggplot(HLA_B07_02[changed,], aes(x=type, y=SB_mod)) + 
  geom_violin() + geom_boxplot(width=0.1)

ggplot(HLA_B07_02[changed,], aes(x=type, y=WB_mod)) + 
  geom_violin() + geom_boxplot(width=0.1)

ggplot(HLA_B07_02[changed,], aes(x=type, y=pep_len)) + 
  geom_violin() + geom_boxplot(width=0.1)

# 
changed<-(HLA_C07_01$modified == "changed" | HLA_C07_01$modified == "add_reg")
HLA_C07_01_violin<-ggplot(HLA_C07_01[changed,], aes(x=type, y=KS)) + 
  geom_violin() + geom_boxplot(width=0.1)

HLA_C07_01_tumor_out<-boxplot.stats(HLA_C07_01[HLA_C07_01$type[changed]=="tumor","KS"])$out
tumor_out<-HLA_C07_01[HLA_C07_01$type == "tumor","KS"] %in% HLA_C07_01_tumor_out
HLA_C07_01_normal_out<-boxplot.stats(HLA_C07_01[HLA_C07_01$type[changed]=="normal","KS"])$out
normal_out<-HLA_C07_01[HLA_C07_01$type == "normal","KS"] %in% HLA_C07_01_normal_out

HLA_C07_01_cong<-rbind(HLA_C07_01[normal_out,],HLA_C07_01[tumor_out,])
a<-ggplot(HLA_C07_01_cong[HLA_C07_01_cong$type=="normal",], aes(x = SB_mod, y = KS)) + geom_point() + labs(y="KS") + labs(x="Strong Binder Count") + labs(title="Normal Associated")
b<-ggplot(HLA_C07_01_cong[HLA_C07_01_cong$type=="normal",], aes(x = WB_mod, y = KS)) + geom_point() + labs(y="KS") + labs(x="Weak Binder Count") + labs(title="Normal Associated")
c<-ggplot(HLA_C07_01_cong[HLA_C07_01_cong$type=="tumor",], aes(x = SB_mod, y = KS)) + geom_point() + labs(y="KS") + labs(x="Strong Binder Count") + labs(title="Tumor Associated")
d<-ggplot(HLA_C07_01_cong[HLA_C07_01_cong$type=="tumor",], aes(x = WB_mod, y = KS)) + geom_point() + labs(y="KS") + labs(x="Weak Binder Count") + labs(title="Tumor Associated")

ggarrange(a, b, c, d, 
          labels = c("E", "F", "G","H"),
          ncol = 2, nrow = 2)


ggplot(HLA_C07_01[changed,], aes(x=type, y=perc)) + 
  geom_violin() + geom_boxplot(width=0.1)

ggplot(HLA_C07_01[changed,], aes(x=type, y=SB_mod)) + 
  geom_violin() + geom_boxplot(width=0.1)

ggplot(HLA_C07_01[changed,], aes(x=type, y=WB_mod)) + 
  geom_violin() + geom_boxplot(width=0.1)

ggplot(HLA_C07_01[changed,], aes(x=type, y=pep_len)) + 
  geom_violin() + geom_boxplot(width=0.1)

# combining alleles and plotting violin
HLA_A02_01$HLA<-"HLA-A02:01"
HLA_B07_02$HLA<-"HLA-B07:02"
HLA_C07_01$HLA<-"HLA-C07:01"
HLA_comb<-rbind(HLA_A02_01,HLA_B07_02,HLA_C07_01)
changed<-(HLA_comb$modified == "changed" | HLA_comb$modified == "add_reg")
HLA_comb_violin<-ggplot(HLA_comb[changed,], aes(x=type, y=KS, fill=HLA)) + 
  geom_violin(position=position_dodge(1)) + 
  geom_boxplot(width=0.1,position=position_dodge(1)) + 
  labs(y="Kolmogorov-Smirnov of Reference and Modified Kmer Distribution") + labs(x="") + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))
HLA_comb_violin

# normal cds
changed<-!(HLA_A02_01$modified == "changed" | HLA_A02_01$modified == "add_reg")
ggplot(HLA_A02_01[changed,], aes(x=type, y=KS)) + 
  geom_violin() + geom_boxplot(width=0.1)

ggplot(HLA_A02_01[changed,], aes(x=type, y=perc)) + 
  geom_violin() + geom_boxplot(width=0.1)

ggplot(HLA_A02_01[changed,], aes(x=type, y=SB_mod)) + 
  geom_violin() + geom_boxplot(width=0.1)

ggplot(HLA_A02_01[changed,], aes(x=type, y=WB_mod)) + 
  geom_violin() + geom_boxplot(width=0.1)

ggplot(HLA_A02_01[changed,], aes(x=type, y=pep_len)) + 
  geom_violin() + geom_boxplot(width=0.1)


# 
changed<-!(HLA_B07_02$modified == "changed" | HLA_B07_02$modified == "add_reg")
ggplot(HLA_B07_02[changed,], aes(x=type, y=KS)) + 
  geom_violin() + geom_boxplot(width=0.1)

ggplot(HLA_B07_02[changed,], aes(x=type, y=perc)) + 
  geom_violin() + geom_boxplot(width=0.1)

ggplot(HLA_B07_02[changed,], aes(x=type, y=SB_mod)) + 
  geom_violin() + geom_boxplot(width=0.1)

ggplot(HLA_B07_02[changed,], aes(x=type, y=WB_mod)) + 
  geom_violin() + geom_boxplot(width=0.1)

ggplot(HLA_B07_02[changed,], aes(x=type, y=pep_len)) + 
  geom_violin() + geom_boxplot(width=0.1)

# 
changed<-!(HLA_C07_01$modified == "changed" | HLA_C07_01$modified == "add_reg")
ggplot(HLA_C07_01[changed,], aes(x=type, y=KS)) + 
  geom_violin() + geom_boxplot(width=0.1)

ggplot(HLA_C07_01[changed,], aes(x=type, y=perc)) + 
  geom_violin() + geom_boxplot(width=0.1)

ggplot(HLA_C07_01[changed,], aes(x=type, y=SB_mod)) + 
  geom_violin() + geom_boxplot(width=0.1)

ggplot(HLA_C07_01[changed,], aes(x=type, y=WB_mod)) + 
  geom_violin() + geom_boxplot(width=0.1)

ggplot(HLA_C07_01[changed,], aes(x=type, y=pep_len)) + 
  geom_violin() + geom_boxplot(width=0.1)

#------------------------------------------------------------------------------------------------------------------------------------------------#
# working with metric data all

#plotting peptide length HLA_A02_01
ggplot(HLA_A02_01[HLA_A02_01$pep_len < 10000,],aes(x=pep_len, y=KS)) + geom_point()

ggplot(HLA_A02_01[HLA_A02_01$pep_len < 10000,],aes(x=pep_len, y=perc)) + geom_point()

ggplot(HLA_A02_01[HLA_A02_01$pep_len < 10000,],aes(x=pep_len, y=SB_mod)) + geom_point()

ggplot(HLA_A02_01[HLA_A02_01$pep_len < 10000,],aes(x=pep_len, y=WB_mod)) + geom_point()

ggplot(HLA_A02_01[HLA_A02_01$pep_len < 10000,],aes(x=pep_len, y=SB_ref)) + geom_point()

ggplot(HLA_A02_01[HLA_A02_01$pep_len < 10000,],aes(x=pep_len, y=WB_ref)) + geom_point()


#plotting peptide length HLA_B07_02
ggplot(HLA_B07_02[HLA_B07_02$pep_len < 10000,],aes(x=pep_len, y=KS)) + geom_point()

ggplot(HLA_B07_02[HLA_B07_02$pep_len < 10000,],aes(x=pep_len, y=perc)) + geom_point()

ggplot(HLA_B07_02[HLA_B07_02$pep_len < 10000,],aes(x=pep_len, y=SB_mod)) + geom_point()

ggplot(HLA_B07_02[HLA_B07_02$pep_len < 10000,],aes(x=pep_len, y=WB_mod)) + geom_point()

ggplot(HLA_B07_02[HLA_B07_02$pep_len < 10000,],aes(x=pep_len, y=SB_ref)) + geom_point()

ggplot(HLA_B07_02[HLA_B07_02$pep_len < 10000,],aes(x=pep_len, y=WB_ref)) + geom_point()

#plotting peptide length HLA_C07_01
ggplot(HLA_C07_01[HLA_C07_01$pep_len < 10000,],aes(x=pep_len, y=KS)) + geom_point()

ggplot(HLA_C07_01[HLA_C07_01$pep_len < 10000,],aes(x=pep_len, y=perc)) + geom_point()

ggplot(HLA_C07_01[HLA_C07_01$pep_len < 10000,],aes(x=pep_len, y=SB_mod)) + geom_point()

ggplot(HLA_C07_01[HLA_C07_01$pep_len < 10000,],aes(x=pep_len, y=WB_mod)) + geom_point()

ggplot(HLA_C07_01[HLA_C07_01$pep_len < 10000,],aes(x=pep_len, y=SB_ref)) + geom_point()

ggplot(HLA_C07_01[HLA_C07_01$pep_len < 10000,],aes(x=pep_len, y=WB_ref)) + geom_point()

#------------------------------------------------------------------------------------------------------------------------------------------------#
# working with metric data psi > 0

HLA_A02_01_tumor<-HLA_A02_01[HLA_A02_01$deltapsi > 0,]
HLA_B07_02_tumor<-HLA_B07_02[HLA_B07_02$deltapsi > 0,]
HLA_C07_01_tumor<-HLA_C07_01[HLA_C07_01$deltapsi > 0,]

#plotting peptide length HLA_A02_01
ggplot(HLA_A02_01_tumor[HLA_A02_01_tumor$pep_len < 10000,],aes(x=pep_len, y=KS)) + geom_point()

ggplot(HLA_A02_01_tumor[HLA_A02_01_tumor$pep_len < 10000,],aes(x=pep_len, y=perc)) + geom_point()

ggplot(HLA_A02_01_tumor[HLA_A02_01_tumor$pep_len < 10000,],aes(x=pep_len, y=SB_mod)) + geom_point()

ggplot(HLA_A02_01_tumor[HLA_A02_01_tumor$pep_len < 10000,],aes(x=pep_len, y=WB_mod)) + geom_point()

ggplot(HLA_A02_01_tumor[HLA_A02_01_tumor$pep_len < 10000,],aes(x=pep_len, y=SB_ref)) + geom_point()

ggplot(HLA_A02_01_tumor[HLA_A02_01_tumor$pep_len < 10000,],aes(x=pep_len, y=WB_ref)) + geom_point()


#plotting peptide length HLA_B07_02
ggplot(HLA_B07_02_tumor[HLA_B07_02_tumor$pep_len < 10000,],aes(x=pep_len, y=KS)) + geom_point()

ggplot(HLA_B07_02_tumor[HLA_B07_02_tumor$pep_len < 10000,],aes(x=pep_len, y=perc)) + geom_point()

ggplot(HLA_B07_02_tumor[HLA_B07_02_tumor$pep_len < 10000,],aes(x=pep_len, y=SB_mod)) + geom_point()

ggplot(HLA_B07_02_tumor[HLA_B07_02_tumor$pep_len < 10000,],aes(x=pep_len, y=WB_mod)) + geom_point()

ggplot(HLA_B07_02_tumor[HLA_B07_02_tumor$pep_len < 10000,],aes(x=pep_len, y=SB_ref)) + geom_point()

ggplot(HLA_B07_02_tumor[HLA_B07_02_tumor$pep_len < 10000,],aes(x=pep_len, y=WB_ref)) + geom_point()

#plotting peptide length HLA_C07_01
ggplot(HLA_C07_01_tumor[HLA_C07_01_tumor$pep_len < 10000,],aes(x=pep_len, y=KS)) + geom_point()

ggplot(HLA_C07_01_tumor[HLA_C07_01_tumor$pep_len < 10000,],aes(x=pep_len, y=perc)) + geom_point()

ggplot(HLA_C07_01_tumor[HLA_C07_01_tumor$pep_len < 10000,],aes(x=pep_len, y=SB_mod)) + geom_point()

ggplot(HLA_C07_01_tumor[HLA_C07_01_tumor$pep_len < 10000,],aes(x=pep_len, y=WB_mod)) + geom_point()

ggplot(HLA_C07_01_tumor[HLA_C07_01_tumor$pep_len < 10000,],aes(x=pep_len, y=SB_ref)) + geom_point()

ggplot(HLA_C07_01_tumor[HLA_C07_01_tumor$pep_len < 10000,],aes(x=pep_len, y=WB_ref)) + geom_point()

#------------------------------------------------------------------------------------------------------------------------------------------------#
# working with metric data psi < 0

HLA_A02_01_normal<-HLA_A02_01[HLA_A02_01$deltapsi < 0,]
HLA_B07_02_normal<-HLA_B07_02[HLA_B07_02$deltapsi < 0,]
HLA_C07_01_normal<-HLA_C07_01[HLA_C07_01$deltapsi < 0,]

#plotting peptide length HLA_A02_01
ggplot(HLA_A02_01_normal[HLA_A02_01_normal$pep_len < 10000,],aes(x=pep_len, y=KS)) + geom_point()

ggplot(HLA_A02_01_normal[HLA_A02_01_normal$pep_len < 10000,],aes(x=pep_len, y=perc)) + geom_point()

ggplot(HLA_A02_01_normal[HLA_A02_01_normal$pep_len < 10000,],aes(x=pep_len, y=SB_mod)) + geom_point()

ggplot(HLA_A02_01_normal[HLA_A02_01_normal$pep_len < 10000,],aes(x=pep_len, y=WB_mod)) + geom_point()

ggplot(HLA_A02_01_normal[HLA_A02_01_normal$pep_len < 10000,],aes(x=pep_len, y=SB_ref)) + geom_point()

ggplot(HLA_A02_01_normal[HLA_A02_01_normal$pep_len < 10000,],aes(x=pep_len, y=WB_ref)) + geom_point()


#plotting peptide length HLA_B07_02
ggplot(HLA_B07_02_normal[HLA_B07_02_normal$pep_len < 10000,],aes(x=pep_len, y=KS)) + geom_point()

ggplot(HLA_B07_02_normal[HLA_B07_02_normal$pep_len < 10000,],aes(x=pep_len, y=perc)) + geom_point()

ggplot(HLA_B07_02_normal[HLA_B07_02_normal$pep_len < 10000,],aes(x=pep_len, y=SB_mod)) + geom_point()

ggplot(HLA_B07_02_normal[HLA_B07_02_normal$pep_len < 10000,],aes(x=pep_len, y=WB_mod)) + geom_point()

ggplot(HLA_B07_02_normal[HLA_B07_02_normal$pep_len < 10000,],aes(x=pep_len, y=SB_ref)) + geom_point()

ggplot(HLA_B07_02_normal[HLA_B07_02_normal$pep_len < 10000,],aes(x=pep_len, y=WB_ref)) + geom_point()

#plotting peptide length HLA_C07_01
ggplot(HLA_C07_01_normal[HLA_C07_01_normal$pep_len < 10000,],aes(x=pep_len, y=KS)) + geom_point()

ggplot(HLA_C07_01_normal[HLA_C07_01_normal$pep_len < 10000,],aes(x=pep_len, y=perc)) + geom_point()

ggplot(HLA_C07_01_normal[HLA_C07_01_normal$pep_len < 10000,],aes(x=pep_len, y=SB_mod)) + geom_point()

ggplot(HLA_C07_01_normal[HLA_C07_01_normal$pep_len < 10000,],aes(x=pep_len, y=WB_mod)) + geom_point()

ggplot(HLA_C07_01_normal[HLA_C07_01_normal$pep_len < 10000,],aes(x=pep_len, y=SB_ref)) + geom_point()

ggplot(HLA_C07_01_normal[HLA_C07_01_normal$pep_len < 10000,],aes(x=pep_len, y=WB_ref)) + geom_point()

#------------------------------------------------------------------------------------------------------------------------------------------------#
# working with metric data psi > 0 and cds changed

HLA_A02_01_tumor<-HLA_A02_01[HLA_A02_01$deltapsi > 0 & HLA_A02_01$modified == "changed",]
HLA_B07_02_tumor<-HLA_B07_02[HLA_B07_02$deltapsi > 0 & HLA_A02_01$modified == "changed",]
HLA_C07_01_tumor<-HLA_C07_01[HLA_C07_01$deltapsi > 0 & HLA_A02_01$modified == "changed",]

#plotting peptide length HLA_A02_01

plot(HLA_A02_01_tumor[HLA_A02_01_tumor$pep_len < 10000,13:19])
plot(HLA_A02_01_tumor[HLA_A02_01_tumor$pep_len < 10000,c(13:16,19)])

#plotting peptide length HLA_B07_02

plot(HLA_B07_02_tumor[HLA_B07_02_tumor$pep_len < 10000,13:19])
plot(HLA_B07_02_tumor[HLA_B07_02_tumor$pep_len < 10000,c(13:16,19)])

#plotting peptide length HLA_C07_01

plot(HLA_C07_01_tumor[HLA_C07_01_tumor$pep_len < 10000,13:19])
plot(HLA_C07_01_tumor[HLA_C07_01_tumor$pep_len < 10000,c(13:16,19)])

ggplot(ToothGrowth, aes(x=dose, y=len)) + 
  geom_violin() + geom_boxplot(width=0.1)

#------------------------------------------------------------------------------------------------------------------------------------------------#
# working with metric data psi > 0 and cds normal

HLA_A02_01_tumor<-HLA_A02_01[HLA_A02_01$deltapsi > 0 & HLA_A02_01$modified == "normal",]
HLA_B07_02_tumor<-HLA_B07_02[HLA_B07_02$deltapsi > 0 & HLA_A02_01$modified == "normal",]
HLA_C07_01_tumor<-HLA_C07_01[HLA_C07_01$deltapsi > 0 & HLA_A02_01$modified == "normal",]

#plotting peptide length HLA_A02_01

plot(HLA_A02_01_tumor[HLA_A02_01_tumor$pep_len < 10000,13:19])
plot(HLA_A02_01_tumor[HLA_A02_01_tumor$pep_len < 10000,c(13:16,19)])

plot(HLA_B07_02_tumor[HLA_B07_02_tumor$pep_len < 10000,13:19])
plot(HLA_B07_02_tumor[HLA_B07_02_tumor$pep_len < 10000,c(13:16,19)])


plot(HLA_C07_01_tumor[HLA_C07_01_tumor$pep_len < 10000,13:19])
plot(HLA_C07_01_tumor[HLA_C07_01_tumor$pep_len < 10000,c(13:16,19)])

#------------------------------------------------------------------------------------------------------------------------------------------------#
# working with metric data psi < 0 and cds changed

HLA_A02_01_tumor<-HLA_A02_01[HLA_A02_01$deltapsi < 0 & HLA_A02_01$modified == "changed",]
HLA_B07_02_tumor<-HLA_B07_02[HLA_B07_02$deltapsi < 0 & HLA_A02_01$modified == "changed",]
HLA_C07_01_tumor<-HLA_C07_01[HLA_C07_01$deltapsi < 0 & HLA_A02_01$modified == "changed",]

#plotting peptide length HLA_A02_01

plot(HLA_A02_01_tumor[HLA_A02_01_tumor$pep_len < 10000,13:19])

#plotting peptide length HLA_B07_02

plot(HLA_B07_02_tumor[HLA_B07_02_tumor$pep_len < 10000,13:19])

#plotting peptide length HLA_C07_01

plot(HLA_C07_01_tumor[HLA_C07_01_tumor$pep_len < 10000,13:19])

#------------------------------------------------------------------------------------------------------------------------------------------------#
# working with metric data psi < 0 and cds normal

HLA_A02_01_tumor<-HLA_A02_01[HLA_A02_01$deltapsi < 0 & HLA_A02_01$modified == "normal",]
HLA_B07_02_tumor<-HLA_B07_02[HLA_B07_02$deltapsi < 0 & HLA_A02_01$modified == "normal",]
HLA_C07_01_tumor<-HLA_C07_01[HLA_C07_01$deltapsi < 0 & HLA_A02_01$modified == "normal",]

#plotting peptide length HLA_A02_01

plot(HLA_A02_01_tumor[HLA_A02_01_tumor$pep_len < 10000,13:19])


plot(HLA_B07_02_tumor[HLA_B07_02_tumor$pep_len < 10000,13:19])


plot(HLA_C07_01_tumor[HLA_C07_01_tumor$pep_len < 10000,13:19])
