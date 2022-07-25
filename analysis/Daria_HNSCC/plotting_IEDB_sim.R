#------------------------------------------------------------------------------#
# loading in libraries

library(ggplot2)

#------------------------------------------------------------------------------#
# reading in the filenames file 

filenames <- read.table("/media/theron/My_Passport/head_and_neck_DARIA/data/splicemutr_05_26_2021/split_peps/filenames.txt")

#------------------------------------------------------------------------------#
# conglomerating IEDB_peps

for (i in seq(nrow(filenames))){
  print(i)
  if (i == 1){
    sim_score <- read.table(filenames[i,])
  } else {
    sim_score_fill <- read.table(filenames[i,])
    sim_score <- rbind(sim_score,sim_score_fill)
  }
}
sim_score$V2 <- as.numeric(sim_score$V2)

#------------------------------------------------------------------------------#
# plotting sim_score histogram

ggplot(sim_score,aes(x=V2))+geom_histogram(binwidth=0.01)+
  scale_y_log10()+
  ylab("log10(count)")+
  xlab("R")+theme(text = element_text(size = 15))
ggsave("/media/theron/My_Passport/SARE_2021/project/spliemutr_images/R_vals_2.png",
       width=7, height = 7, units = "in",dpi=300)
ggplot(sim_score,aes(x=V2))+geom_histogram(binwidth=0.01)+
  xlab("R")+theme(text = element_text(size = 15))
ggsave("/media/theron/My_Passport/SARE_2021/project/spliemutr_images/R_vals_2.png",
       width=10, height = 10, units = "in",dpi=300)