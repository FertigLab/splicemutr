#------------------------------------------------------------------------------#
# internal functions

orf<-function(seq,i){
  # finds the largest open reading frame a sequence
  # inputs:
  #   seq: the DNA sequence to search Overrepresented_sequences
  #   i: the frame to use
  # output:
  #   out: the largest orf, 0 if there is a stop before the first start, 1 if there is no start or no stop
  
  out<-0
  if ( all(unique(strsplit(seq,"")[[1]]) %in% c("A","T","C","G")) ) {
    seq<-DNAStringSet(substr(seq,i,nchar(seq)-((nchar(seq)-i+1) %% 3)))
    trans<-translate(seq,no.init.codon = FALSE)
  } else {
    trans<-seq
  }
  a<-strsplit(as.character(trans),"")
  stop<-which((a[[1]] %in% "*") == TRUE)
  start<-which((a[[1]] %in% "M") == TRUE)
  out<-list(out,trans,NA,NA)
  if (length(start)>0 & length(stop)>0) {
    if (stop[1] < start[1]){
      out<-1
      out<-list(out,trans,NA,paste(stop,collapse=","))
    } else if (length(stop)>1){
      out<-2
      trans_small<-substr(as.character(trans),start[1],stop[1])
      out<-list(out,trans,as.character(trans_small),paste(stop,collapse=","))
    } else if (stop[1] != length(a[[1]])){
      out<-3
      out<-list(out,trans,NA,paste(stop,collapse=","))
    } else if (start[1] != 1){
      out<-4
      out<-list(out,trans,NA,paste(stop,collapse=","))
    } else {
      out<-"tx"
      trans_small<-substr(as.character(trans),start[1],stop[1])
      return(list(out,trans,as.character(trans_small),paste(stop,collapse=",")))
    }
  } else if (length(start)>0 & !length(stop)>0){
    out<-5
    out<-list(out,trans,NA,paste(stop,collapse=","))
  } else if (!length(start)>0 & length(stop)>0){
    out<-6
    out<-list(out,trans,NA,paste(stop,collapse=","))
  }
  return(out)
}


#------------------------------------------------------------------------------------------------------------------------------------------------#
# looking at errors and how they relate to immunogenicity

HLA_A02_01<-readRDS("/media/theron/My_Passport/head_and_neck_DARIA/data/leafcutter_10_21_2020/kmers/processed/HLA-A02-01_kmers_9_metrics_all.rds")
# HLA_A02_01<-HLA_A02_01[!is.na(HLA_A02_01$KS),]
HLA_B07_02<-readRDS("/media/theron/My_Passport/head_and_neck_DARIA/data/leafcutter_10_21_2020/kmers/processed/HLA-B07-02_kmers_9_metrics_all.rds")
# HLA_B07_02<-HLA_B07_02[!is.na(HLA_B07_02$KS),]
HLA_C07_01<-readRDS("/media/theron/My_Passport/head_and_neck_DARIA/data/leafcutter_10_21_2020/kmers/processed/HLA-C07-01_kmers_9_metrics_all.rds")
# HLA_C07_01<-HLA_C07_01[!is.na(HLA_C07_01$KS),]

#------------------------------------------------------------------------------------------------------------------------------------------------#
# reprocessing HLA_data
data<-HLA_B07_02[,1:10]
HLA_B07_02_re<-data[1,]
HLA_B07_02_re$error<-NA
HLA_B07_02_re$original_pep<-NA
HLA_B07_02_re$peptide<-NA
HLA_B07_02_re$stops<-NA

iter<-1
for (i in 1:nrow(data)) {
  print(sprintf("%d rows complete out of %d total",i,nrow(data)))
  for (j in 1:3){
    orf_dat<-orf(data$sequence[i],j)
    HLA_A02_01_re[iter,]<-c(data[i,],orf_dat[[1]],orf_dat[[2]],orf_dat[[3]],paste(orf_dat[[4]]),j)
    iter<-iter+1
  }
}

#------------------------------------------------------------------------------------------------------------------------------------------------#
# From first frame, translatable based on junction type

vapply(HLA_B07_02$peptide,)

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