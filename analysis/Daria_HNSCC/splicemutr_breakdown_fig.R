# Theron Palmer
# Created: 11/19/2020
# Purpose: To create the plot used to characterize the junction types in splicemutr

#-------------------------------------------------------------------------------#
# loading in necessary libraries

library(plotly)
library(dplyr)
library(ggpubr)

#-------------------------------------------------------------------------------#
# internal functions

get_splicemutr_breakdown<-function(mod_dat,verdict,verdicts){
  # mod_dat is the modified data containing verdict and modified peptide information
  if (verdict == "all"){
    ann_dat<-data.frame(table(mod_dat$errors))
    colnames(ann_dat)<-c("error","freq")
  } else {
    verdict_data<-mod_dat[mod_dat$verdict == verdict,]
    ann_dat<-data.frame(table(verdict_data$errors))
    colnames(ann_dat)<-c("error","freq")
  }
  return(ann_dat)
}

get_tumor_normal_breakdown<-function(mod_dat,verdict,verdicts){
  all_dat<-data.frame()
  all_dat[seq(2),seq(6)]<-c(0,0,0,0,0,0)
  row.names(all_dat)<-c("N","T")
  colnames(all_dat)<-c("zero","one","two","three","four","five")
  if (verdict == "all"){
    all_dat[seq(2),seq(6)]<-vapply(seq(0,5),function(x){
      a<-length(which(mod_dat$errors[mod_dat$status == "N"] == x))
      b<-length(which(mod_dat$errors[mod_dat$status == "T"] == x))
      c(a,b)
    },numeric(2))
    all_dat$status<-c("Normal","Tumor")
  } else {
    verdict_dat<-mod_dat[mod_dat$verdict == verdict,]
    all_dat[seq(2),seq(6)]<-vapply(seq(0,5),function(x){
      a<-length(which(verdict_dat$errors[verdict_dat$status == "N"] == x))
      b<-length(which(verdict_dat$errors[verdict_dat$status == "T"] == x))
      c(a,b)
    },numeric(2))
    all_dat$status<-c("Normal","Tumor")
  }
  return(all_dat)
}

#-------------------------------------------------------------------------------#
# loading in peptide data

ref_dat<-readRDS("/media/theron/My_Passport/head_and_neck_DARIA/data/leafcutter_10_21_2020/ref/peptides/peptides.rds") # reading in the reference data
mod_dat<-readRDS("/media/theron/My_Passport/head_and_neck_DARIA/data/leafcutter_10_21_2020/introns/peps.rds")
mod_dat<-mod_dat[mod_dat$frame==1,]
errors<-c()
status<-c()
errors<-vapply(seq(length(mod_dat$peptide)),function(x){
  if(nchar(mod_dat$peptide[x]) > 1){
    return(5)
  } else {
    return(as.numeric(mod_dat$peptide[x]))
  }
},numeric(1))
status<-vapply(seq(length(mod_dat$deltapsi)),function(x){
  if(mod_dat$deltapsi[x] > 0){
    return("T")
  } else {
    return("N")
  }
},character(1))
mod_dat$errors<-errors
mod_dat$status<-status
verdicts<-unique(mod_dat$verdict)

annotated<-get_splicemutr_breakdown(mod_dat,"annotated",verdicts)
annotated_status<-get_tumor_normal_breakdown(mod_dat,"annotated",verdicts)

cryptic_fiveprime<-get_splicemutr_breakdown(mod_dat,"cryptic_fiveprime",verdicts)
cryptic_fiveprime_status<-get_tumor_normal_breakdown(mod_dat,"cryptic_fiveprime",verdicts)

cryptic_threeprime<-get_splicemutr_breakdown(mod_dat,"cryptic_threeprime",verdicts)
cryptic_threeprime_status<-get_tumor_normal_breakdown(mod_dat,"cryptic_threeprime",verdicts)

cryptic_unanchored<-get_splicemutr_breakdown(mod_dat,"cryptic_unanchored",verdicts)
cryptic_unanchored_status<-get_tumor_normal_breakdown(mod_dat,"cryptic_unanchored",verdicts)

novel_annotated_pair<-get_splicemutr_breakdown(mod_dat,"novel annotated pair",verdicts)
novel_annotated_pair_status<-get_tumor_normal_breakdown(mod_dat,"novel annotated pair",verdicts)

all<-get_splicemutr_breakdown(mod_dat,"all",verdicts)
all_status<-get_tumor_normal_breakdown(mod_dat,"all",verdicts)

data<-annotated
data<-cbind(data,cryptic_fiveprime$freq,cryptic_threeprime$freq,
            cryptic_unanchored$freq,novel_annotated_pair$freq,all$freq)
colnames(data)<-c("error","annotated","cryptic_fiveprime","cryptic_threeprime",
                  "cryptic_unanchored","novel_annotated_pair","all")
fig <- plot_ly(data, labels = ~error, values = ~annotated, type = 'pie',
               textinfo = "text + values", domain = list(row = 0, column = 0))
fig <- fig %>% add_pie(data, labels = ~error, values = ~cryptic_fiveprime, 
                       textinfo = "text + values", domain = list(row = 1, column = 0))
fig <- fig %>% add_pie(data, labels = ~error, values = ~cryptic_threeprime, 
                       textinfo = "text + values", domain = list(row = 2, column = 0))
fig <- fig %>% add_pie(data, labels = ~error, values = ~cryptic_unanchored, 
                       textinfo = "text + values", domain = list(row = 3, column = 0))
fig <- fig %>% add_pie(data, labels = ~error, values = ~novel_annotated_pair, 
                       textinfo = "text + values", domain = list(row = 4, column = 0))
fig <- fig %>% add_pie(data, labels = ~error, values = ~all,
                       textinfo = "text + values", domain = list(row = 5, column = 1))
fig <- fig %>% layout(showlegend = T,
                      grid=list(rows=5, columns=2))
fig

fig_status <- plot_ly(annotated_status, labels = ~status, values = ~zero, type = 'pie',
               textinfo = "text + values", domain = list(row = 0, column = 0))
fig_status <- fig_status %>% add_pie(annotated_status, labels = ~status, values = ~one, 
                           textinfo = "text + values", domain = list(row = 0, column = 1))
fig_status <- fig_status %>% add_pie(annotated_status, labels = ~status, values = ~two, 
                           textinfo = "text + values", domain = list(row = 0, column = 2))
fig_status <- fig_status %>% add_pie(annotated_status, labels = ~status, values = ~three, 
                           textinfo = "text + values", domain = list(row = 0, column = 3))
fig_status <- fig_status %>% add_pie(annotated_status, labels = ~status, values = ~four, 
                           textinfo = "text + values", domain = list(row = 0, column = 4))
fig_status <- fig_status %>% add_pie(annotated_status, labels = ~status, values = ~five, 
                           textinfo = "text + values", domain = list(row = 0, column = 5))

fig_status <- fig_status %>% add_pie(cryptic_fiveprime_status, labels = ~status, values = ~zero, 
                       textinfo = "text + values", domain = list(row = 1, column = 0))
fig_status <- fig_status %>% add_pie(cryptic_fiveprime_status, labels = ~status, values = ~one, 
                           textinfo = "text + values", domain = list(row = 1, column = 1))
fig_status <- fig_status %>% add_pie(cryptic_fiveprime_status, labels = ~status, values = ~two, 
                           textinfo = "text + values", domain = list(row = 1, column = 2))
fig_status <- fig_status %>% add_pie(cryptic_fiveprime_status, labels = ~status, values = ~three, 
                           textinfo = "text + values", domain = list(row = 1, column = 3))
fig_status <- fig_status %>% add_pie(cryptic_fiveprime_status, labels = ~status, values = ~four, 
                           textinfo = "text + values", domain = list(row = 1, column = 4))
fig_status <- fig_status %>% add_pie(cryptic_fiveprime_status, labels = ~status, values = ~five, 
                           textinfo = "text + values", domain = list(row = 1, column = 5))

fig_status <- fig_status %>% add_pie(cryptic_threeprime_status, labels = ~status, values = ~zero, 
                                     textinfo = "text + values", domain = list(row = 2, column = 0))
fig_status <- fig_status %>% add_pie(cryptic_threeprime_status, labels = ~status, values = ~one, 
                                     textinfo = "text + values", domain = list(row = 2, column = 1))
fig_status <- fig_status %>% add_pie(cryptic_threeprime_status, labels = ~status, values = ~two, 
                                     textinfo = "text + values", domain = list(row = 2, column = 2))
fig_status <- fig_status %>% add_pie(cryptic_threeprime_status, labels = ~status, values = ~three, 
                                     textinfo = "text + values", domain = list(row = 2, column = 3))
fig_status <- fig_status %>% add_pie(cryptic_threeprime_status, labels = ~status, values = ~four, 
                                     textinfo = "text + values", domain = list(row = 2, column = 4))
fig_status <- fig_status %>% add_pie(cryptic_threeprime_status, labels = ~status, values = ~five, 
                                     textinfo = "text + values", domain = list(row = 2, column = 5))

fig_status <- fig_status %>% add_pie(cryptic_unanchored_status, labels = ~status, values = ~zero, 
                                     textinfo = "text + values", domain = list(row = 3, column = 0))
fig_status <- fig_status %>% add_pie(cryptic_unanchored_status, labels = ~status, values = ~one, 
                                     textinfo = "text + values", domain = list(row = 3, column = 1))
fig_status <- fig_status %>% add_pie(cryptic_unanchored_status, labels = ~status, values = ~two, 
                                     textinfo = "text + values", domain = list(row = 3, column = 2))
fig_status <- fig_status %>% add_pie(cryptic_unanchored_status, labels = ~status, values = ~three, 
                                     textinfo = "text + values", domain = list(row = 3, column = 3))
fig_status <- fig_status %>% add_pie(cryptic_unanchored_status, labels = ~status, values = ~four, 
                                     textinfo = "text + values", domain = list(row = 3, column = 4))
fig_status <- fig_status %>% add_pie(cryptic_unanchored_status, labels = ~status, values = ~five, 
                                     textinfo = "text + values", domain = list(row = 3, column = 5))

fig_status <- fig_status %>% add_pie(novel_annotated_pair, labels = ~status, values = ~zero, 
                                     textinfo = "text + values", domain = list(row = 4, column = 0))
fig_status <- fig_status %>% add_pie(novel_annotated_pair, labels = ~status, values = ~one, 
                                     textinfo = "text + values", domain = list(row = 4, column = 1))
fig_status <- fig_status %>% add_pie(novel_annotated_pair, labels = ~status, values = ~two, 
                                     textinfo = "text + values", domain = list(row = 4, column = 2))
fig_status <- fig_status %>% add_pie(novel_annotated_pair, labels = ~status, values = ~three, 
                                     textinfo = "text + values", domain = list(row = 4, column = 3))
fig_status <- fig_status %>% add_pie(novel_annotated_pair, labels = ~status, values = ~four, 
                                     textinfo = "text + values", domain = list(row = 4, column = 4))
fig_status <- fig_status %>% add_pie(novel_annotated_pair, labels = ~status, values = ~five, 
                                     textinfo = "text + values", domain = list(row = 4, column = 5))

fig_status <- fig_status %>% layout(legend=list(orientation="h",x=0,y=0),grid=list(rows=5, columns=6))

fig_status

#-------------------------------------------------------------------------------#
# loading in peptide data
