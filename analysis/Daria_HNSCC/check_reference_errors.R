# Theron Palmer
# Created: 11/16/2020
# Purpose: To explore the errors found in the reference transcriptome

#-------------------------------------------------------------------------------#
# loading in the necessary libraries
library(GenomicFeatures)
library(rtracklayer)
library(BSgenome.Hsapiens.GENCODE.GRCh38.p13)

#-------------------------------------------------------------------------------#
# helper functions

get_AA<-function(mod_tx,ref_dat){
  return(ref_dat$peptide[ref_dat$tx==mod_tx])
}

get_DNA<-function(mod_tx,ref_dat){
  return(ref_dat$DNA_sequence[ref_dat$tx==mod_tx])
}

get_gene<-function(ref_tx,mod_dat){
  if (ref_tx %in% unique(mod_dat$tx_id)){
    return(unique(mod_dat$gene[mod_dat$tx_id == ref_tx]))
  } else {
    return("NONE")
  }
}

get_cds<-function(ref_tx,cdsBytranscript){
  if (ref_tx %in% cdsBytranscript){
    return(cdsBytranscript)
  } else {
    return("NONE")
  }
}

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
  if (length(start)>0 & length(stop)>0) {
    if (stop[1] < start[1]){
      out<-1
      return(as.character(out))
    } else if (length(stop)>1){
      out<-2
    } else if (stop[1] != length(a[[1]])){
      out<-3
    } else if (start[1] != 1){
      out<-4
    } else {
      out<-substr(as.character(trans),start[1],stop[1])
      return(as.character(out))
    }
  } else if (length(start)>0 & !length(stop)>0){
    out<-5
  } else if (!length(start)>0 & length(stop)>0){
    out<-6
  }
  return(as.character(out))
}
#-------------------------------------------------------------------------------#
# loading in the data

Bsgenome<-BSgenome.Hsapiens.GENCODE.GRCh38.p13

# since ref_dat only has tx information, must use both ref_dat and mod data
ref_dat<-readRDS("/media/theron/My_Passport/head_and_neck_DARIA/data/leafcutter_10_21_2020/ref/peptides/peptides.rds") # reading in the reference data
ref_dat_err<-ref_dat[nchar(ref_dat$peptide) == 1,]
ref_dat_complete<-ref_dat[nchar(ref_dat$peptide) > 1,]
mod_dat<-readRDS("/media/theron/My_Passport/head_and_neck_DARIA/data/leafcutter_10_21_2020/introns/peps.rds")
mod_dat<-mod_dat[mod_dat$frame == "1",]

# loading in the gtf files
gtf_file<-import("/media/theron/My_Passport/reference_genomes/GTF_GFF/GENCODE/gencode.v35_copy.gtf") #importing the gtf file for the 
gtf_frame<-as.data.frame(gtf_file)
txdb_gtf<-makeTxDbFromGFF("/media/theron/My_Passport/reference_genomes/GTF_GFF/GENCODE/gencode.v35_copy.gtf") # making the txdb from gtf
cdsBytranscript_gtf<-cdsBy(txdb_gtf, by="tx", use.names=FALSE)


gff_file<-import("/media/theron/My_Passport/reference_genomes/GTF_GFF/GENCODE/gencode.v35.annotation.gff3")
gff_frame<-as.data.frame(gff_file)
txdb_gff<-makeTxDbFromGFF("/media/theron/My_Passport/reference_genomes/GTF_GFF/GENCODE/gencode.v35_copy.gff3") # making the txdb from gtf
cdsBytranscript_gff<-cdsBy(txdb_gff, by="tx", use.names=FALSE)

# annotating the ref_dat_err with gene information
genes<-apply(ref_dat_err,1,function(x){
  get_gene(x[1],mod_dat)
})

ref_dat_err$gene<-genes
ref_dat_err_wg=ref_dat_err[ref_dat_err$gene != "NONE",]

#-------------------------------------------------------------------------------#
# PUSL1

# PUSL1: PUSL1-201
ensembl_NT<-"ATGAGTTCGGCGCCGGCCTCAGGCTCCGTGCGCGCGCGCTATCTTGTGTACTTCCAGTACGTGGGCACCG
ACTTTAACGGGGTCGCGGCCGTCAGGGGCACTCAGCGCGCCGTCGGGGTCCAGAACTACCTGGAGGAGGC
CGCCGAGCGGCTGAATTCCGTGGAGCCGGTCAGGTTCACCATCTCCAGCCGCACGGACGCCGGGGTCCAC
GCCCTGAGCAACGCGGCGCACCTGGACGTCCAGCGCCGCTCAGGCCGGCCGCCCTTCCCGCCCGAGGTCC
TGGCCGAGGCCCTCAACACACACCTGCGGCACCCGGCCATCAGGGTCCTGCGGGCCTTCCGAGTGCCCAG
CGACTTCCACGCTCGTCACGCAGCCACGTCCCGGACCTACCTGTACCGCCTGGCCACTGGCTGTCACCGG
CGTGATGAGCTGCCGGTGTTTGAACGCAACCTATGCTGGACTCTCCCGGCAGACTGCCTGGATATGGTCG
CCATGCAGGAAGCCGCCCAGCACCTCCTCGGCACACACGACTTCAGCGCCTTCCAGTCCGCTGGCAGCCC
GGTGCCGAGCCCCGTGCGAACGCTGCGCCGGGTCTCCGTTTCCCCAGGCCAAGCCAGCCCCTTGGTCACC
CCCGAGGAGAGCAGGAAGCTGCGGTTCTGGAACCTGGAGTTTGAGAGCCAGTCTTTCCTGTATAGACAGG
TACGGAGGATGACGGCTGTGCTGGTGGCCGTGGGGCTGGGGGCTTTGGCACCTGCCCAGGTGAAGACGAT
TCTGGAGAGCCAAGATCCCCTGGGCAAGCACCAGACACGTGTAGCCCCAGCCCACGGCTTATTCCTCAAG
TCAGTGCTGTACGGGAACCTCGGTGCTGCCTCCTGCACCCTGCAGGGGCCACAGTTCGGGAGCCACGGAT
GA"
ensembl_AA<-"MSSAPASGSVRARYLVYFQYVGTDFNGVAAVRGTQRAVGVQNYLEEAAERLNSVEPVRFTISSRTDAGVH
ALSNAAHLDVQRRSGRPPFPPEVLAEALNTHLRHPAIRVLRAFRVPSDFHARHAATSRTYLYRLATGCHR
RDELPVFERNLCWTLPADCLDMVAMQEAAQHLLGTHDFSAFQSAGSPVPSPVRTLRRVSVSPGQASPLVT
PEESRKLRFWNLEFESQSFLYRQVRRMTAVLVAVGLGALAPAQVKTILESQDPLGKHQTRVAPAHGLFLK
SVLYGNLGAASCTLQGPQFGSHG"

# PUSL1:PUSL1-203, ensembl
ensemble_NT <- NA
ensembl_AA <- "MVAMQEAAQHLLGTHDFSAFQSAGSPVPSPVRTLRRVSVSPGQASPLVTPEESRKLRFWNLEFESQSFLYRQVRRMTAVLVAVGLGALAPAQVKTILESQDPLGKHQTRVAPAHGLFLKSVLYGN" 
ensembl_AA_split<-strsplit(ensembl_AA,"")[[1]]

# PUSL1:PUSL1-203, splicemutr
cds<-cdsBytranscript_gtf[["174"]]
seq<-getSeq(Bsgenome, cds)
seq<-paste(as.character(seq[1:length(seq)]), collapse="")
trans<-translate(DNAStringSet(seq),no.init.codon = FALSE)
AA<-as.character(trans)
AA_match<- all((strsplit(AA,"")[[1]] == ensembl_AA_split) == TRUE)

#-------------------------------------------------------------------------------#
# ZBTB48-207

ensembl_NT<-NA
ensembl_AA<-"MSEPEAVLTRRKSNVIRKPCAAEPALSAGSLAAEPAENRKGTAVPVECPTCHKKFLSKYYLKVHNRKHTGEKPFECPKCGKCYFRKENLLEHEARNCMNRSEQVFTCSVCQETFRRRMELRVHMVSHTGEMPYKCSSCSQQFMQKKDLQSHMIKLHGAPKPHACPTCAKCFLSRTELQLHEAFKHRGEKLFVCEECGHRASSRNGLQMHIKAKHRNERPHVCEFCSHAFT"
ensembl_AA_split<-strsplit(ensembl_AA,"")[[1]]

cds<-cdsBytranscript_gtf[["567"]]
seq<-getSeq(Bsgenome, cds)
seq<-paste(as.character(seq[1:length(seq)]), collapse="")
trans<-translate(DNAStringSet(seq),no.init.codon = FALSE)
AA<-as.character(trans)
AA_match<- all((strsplit(AA,"")[[1]] == ensembl_AA_split) == TRUE)

# txdb cds match gtf cds
# txdb cds match gff cds
# cds 3' incomplete

#-------------------------------------------------------------------------------#
# FBXO44-207

ensembl_NT<-NA
ensembl_AA<-"MAVGNINELPENILLELFTHVPARQLLLNCRLVCSLWRDLIDLVTLWKRKCLREGFITEDWDQPVADWKIFYFLRSLHRNLLHNPCAEEGFEFWSLDVNGGDEWKVEDLSRDQRKEFPNDQVRSQARLRVQVPAVRSAPVVRARASGDLPARPGDHPAEERCQVE"
ensembl_AA_split<-strsplit(ensembl_AA,"")[[1]]

cds<-cdsBytranscript_gtf[["823"]]
seq<-getSeq(Bsgenome, cds)
seq<-paste(as.character(seq[1:length(seq)]), collapse="")
trans<-translate(DNAStringSet(seq),no.init.codon = FALSE)
AA<-as.character(trans)
AA_match<- all((strsplit(AA,"")[[1]] == ensembl_AA_split) == TRUE)

# txdb cds match gtf cds
# txdb cds match gff cds
# cds 3' incomplete

#-------------------------------------------------------------------------------#
# FBXO44-209

ensembl_NT<-"ATGGCTGTGGGGAACATCAACGAGCTGCCCGAGAACATCCTGCTGGAGCTGTTCACGCACGTGCCCGCCC
GCCAGCTGCTGCTGAACTGCCGCCTGGTCTGCAGCCTCTGGCGGGACCTCATCGACCTCGTGACCCTCTG
GAAACGCAAGTGCCTGCGAGAGGGCTTCATCACTGAGGACTGGGACCAGCCCGTGGCCGACTGGAAGATC
TTCTACTTCTTACGGAGCCTGCACAGGAACCTCCTGCACAACCCGTGCGCTGAAGAGGGGTTCGAGTTCT
GGAGCCTGGATGTGAATGGAGGCGATGAGTGGAAGGTGGAGGATCTCTCTCGAGACCAGAGGAAGGAATT
CCCCAATGACCAGGTTCGCAGCCAGGCCAGATTGCGGGTCCAAGTACCAGCTGTGCGTTCAGCTCCTGTC
GTCCGCGCACGCGCCTCTGGGGACCTTCCAGCCAGACCCGGCGACCATCCAGCAGAAGAGCGATGCCAAG
TGGAGGGAGGTCTCCCACACATTCTCCAACTACCCGCCCGGCGTCCGCTACATCTGGTTTCAGCACGGCG
GCGTGGACACTCATTACTGGGCCGGCTGGTACGGCCCGAGGGTCACCAACAGCAGCATCACCATCGGGCC
CCCGCTGCCCTGACACCCCCTGAGCCCCCATCTGCTGAACCCTGA"
ensembl_AA<-"XHRNLLHNPCAEEGFEFWSLDVNGGDEWKVEDLSRDQRKEFPNDQVRSQARLRVQVPAVRSAPVVRARASGDLPARPGDHPAEERCQVEGGLPHILQLPARRPLHLVSARRRGHSLLGRLVRPEGHQQQHHHRAPAALTPPEPPSAEP"
ensembl_AA_split<-strsplit(ensembl_AA,"")[[1]]

cds<-cdsBytranscript_gtf[["829"]]
seq<-getSeq(Bsgenome, cds)
seq<-paste(as.character(seq[1:length(seq)]), collapse="")
trans<-translate(DNAStringSet(seq),no.init.codon = FALSE)
AA<-as.character(trans)
AA_match<- all((strsplit(AA,"")[[1]] == ensembl_AA_split) == TRUE)

# txdb cds mismatch gtf cds
# txdb cds match gff cds
# cds 5' incomplete

#-------------------------------------------------------------------------------#
# FBXO44-208

ensembl_NT<-NA
ensembl_AA<-"QVPGGGPQGRRVLGGADGYHTAGHRGQGLVRSQARLRVQVPAVRSAPVVRARASGDLPARPGDHPAEERCQVEGGLPHILQLPARRPLHLVSARRRGHSLLGRLVRPEGHQQQHHHRAPAALTPPEPPSAEP*"
ensembl_AA_split<-strsplit(ensembl_AA,"")[[1]]

cds<-cdsBytranscript_gtf[["830"]]
seq<-getSeq(Bsgenome, cds)
seq<-paste(as.character(seq[1:length(seq)]), collapse="")
trans<-translate(DNAStringSet(seq),no.init.codon = FALSE)
AA<-as.character(trans)
AA_match<- all((strsplit(AA,"")[[1]] == ensembl_AA_split) == TRUE)

# txdb cds mismatch gtf cds
# txdb cds match gff cds
# cds 5' incomplete

#-------------------------------------------------------------------------------#
# KAZN-202

ensembl_NT<-NA
ensembl_AA<-"MKEMLAKDLEESQGGKSSEVLSATELRVQLAQKEQELARAKEALQAMKADRKRLKGEKTDLVSQMQQLYATLESREEQLRDFIRNYEQHRKESEDAVKALAKEKDLLEREKWELRRQAKEATDHATALRSQLDLKDNRMKELEAELAMAKQSLATLTKDVPKRHSL"
ensembl_AA_split<-strsplit(ensembl_AA,"")[[1]]

cds<-cdsBytranscript_gtf[["1024"]]
seq<-getSeq(Bsgenome, cds)
seq<-paste(as.character(seq[1:length(seq)]), collapse="")
trans<-translate(DNAStringSet(seq),no.init.codon = FALSE)
AA<-as.character(trans)
AA_match<- all((strsplit(AA,"")[[1]] == ensembl_AA_split) == TRUE)

# txdb cds match gtf cds
# txdb cds match gff cds
# cds 5' incomplete

#-------------------------------------------------------------------------------#
# SPEN-206

ensembl_NT<-NA
ensembl_AA<-"MKEMLAKDLEESQGGKSSEVLSATELRVQLAQKEQELARAKEALQAMKADRKRLKGEKTDLVSQMQQLYATLESREEQLRDFIRNYEQHRKESEDAVKALAKEKDLLEREKWELRRQAKEATDHATALRSQLDLKDNRMKELEAELAMAKQSLATLTKDVPKRHSL"
ensembl_AA_split<-strsplit(ensembl_AA,"")[[1]]

cds<-cdsBytranscript_gtf[["1115"]]
seq<-getSeq(Bsgenome, cds)
seq<-paste(as.character(seq[1:length(seq)]), collapse="")
trans<-translate(DNAStringSet(seq),no.init.codon = FALSE)
AA<-as.character(trans)
AA_match<- all((strsplit(AA,"")[[1]] == ensembl_AA_split) == TRUE)

# txdb cds match gtf cds
# txdb cds match gff cds
# cds 5' incomplete

#-------------------------------------------------------------------------------#
# SPEN-203

ensembl_NT<-NA
ensembl_AA<-"XQASRPTRSPSGSGSRSRSSSSDSISSSSSTSSDSDSSSSSSDDSPARSVQSAAVPAPTSQLLSSLEKDEPRKSFGIKVQNLPVRSTDTSLKDGLFHEFKKFGKVTSVQIHGTSEERYGLVFFRQQEDQEKALTASKGKLFFGMQIEVTAWIGPETESENEFRPLDERIDEFHPKATRTLFIGNLEKTTTYHDLRNIFQRFGEIVDIDIKKVNGVPQYAFLQYCDIASVCKAIKKMDGEYLGNNRLKLGFGKSMPTNCVWLDGLSSNVSDQYLTRHFCRYGPVVKVVFDRLKGMALV"
ensembl_AA_split<-strsplit(ensembl_AA,"")[[1]]

cds<-cdsBytranscript_gtf[["1119"]]
seq<-getSeq(Bsgenome, cds)
seq<-paste(as.character(seq[1:length(seq)]), collapse="")
trans<-translate(DNAStringSet(seq),no.init.codon = TRUE)
AA<-as.character(trans)
AA_match<- all((strsplit(AA,"")[[1]] == ensembl_AA_split) == TRUE)

# txdb cds match gtf cds
# cds 5' and 3' incomplete

#-------------------------------------------------------------------------------#
# EPB41-215

ensembl_NT<-NA
ensembl_AA<-"XDTPTHEDLTKNKERTSESRGLSRLFSSFLKRPKSQVSEEEGKEVESDKEKGEGGQKEIEFGTSLDEEIILKAPIAAPEPELKTDPSLDLHSLSSAETQPAQEELREDPDFEIKEGEGLEECSKIEVKEESPQSKAETELKASQKPIRKHRNMHCKVSLLDDTVYECVVEKHAKGQDLLKRVCEHLNLLEEDYFGLAIWDNATSKTWLDSAKEIKKQVRGVPWNFTFNVKFYPPDPAQLTEDITRYYLCLQLRQDIVAGRLPCSFATLALLGSYTIQSELGDYDPELHGVDYVSDFKLAPNQTKELEEKVMELHKSYRSMTPAQADLEFLENAKKLSMYGVDLHKAKLLKLLGRSSARRGGSHL"
ensembl_AA_split<-strsplit(ensembl_AA,"")[[1]]

cds<-cdsBytranscript_gtf[["1981"]]
seq<-getSeq(Bsgenome, cds)
seq<-paste(as.character(seq[1:length(seq)]), collapse="")
trans<-translate(DNAStringSet(seq),no.init.codon = TRUE)
AA<-as.character(trans)
AA_match<- all((strsplit(AA,"")[[1]] == ensembl_AA_split) == TRUE)

# txdb cds match gtf cds
# nonsense mediated decay; cds 5' incomplete

#-------------------------------------------------------------------------------#
# RNF220-206

ensembl_NT<-NA
ensembl_AA<-"XEGESPTASPHSSATDDLHHSDRYQTFLRVRANRQTRLNARIGKMKRRKQDEGQVCPLCNRPLAGSEQEMSRHVEHCLSKREGSCMAEDDAVDIEHENNNRFEEYEWCGQKRIRATTLLEG"
ensembl_AA_split<-strsplit(ensembl_AA,"")[[1]]

cds<-cdsBytranscript_gtf[["3080"]]
seq<-getSeq(Bsgenome, cds)
seq<-paste(as.character(seq[1:length(seq)]), collapse="")
trans<-translate(DNAStringSet(seq),no.init.codon = TRUE)
AA<-as.character(trans)
AA_match<- all((strsplit(AA,"")[[1]] == ensembl_AA_split) == TRUE)

# txdb cds match gtf cds
# cds 5' incomplete

#-------------------------------------------------------------------------------#
# SCP2-211

ensembl_NT<-NA
ensembl_AA<-"XALQHNLGIGGAVVVTLYKMGFPEAASSFRTHQIEAVPTSSASDGFKANLVFKEIEKKLEEIRRLTAQSQWLTQTSWL"
ensembl_AA_split<-strsplit(ensembl_AA,"")[[1]]

cds<-cdsBytranscript_gtf[["3525"]]
seq<-getSeq(Bsgenome, cds)
seq<-paste(as.character(seq[1:length(seq)]), collapse="")
trans<-translate(DNAStringSet(seq),no.init.codon = TRUE)
AA<-as.character(trans)
AA_match<- all((strsplit(AA,"")[[1]] == ensembl_AA_split) == TRUE)

# txdb cds mismatch gtf cds
# cds 5' incomplete

#-------------------------------------------------------------------------------#
# RSRC1-209

ensembl_NT<-NA
ensembl_AA<-"XKKHRRRSSSSSSSDSRTYSRKKGGRKSRSKSRSWSRDLQPRSHSYDRRRRHRSSSSSSYGSRRKRSRSRSRGRGKSYRVQRSRSKSRTRRYPPTSTSWVAGTTGTHRHTQL"
ensembl_AA_split<-strsplit(ensembl_AA,"")[[1]]

cds<-cdsBytranscript_gtf[["42375"]]
seq<-getSeq(Bsgenome, cds)
seq<-paste(as.character(seq[1:length(seq)]), collapse="")
trans<-translate(DNAStringSet(seq),no.init.codon = TRUE)
AA<-as.character(trans)
AA_match<- all((strsplit(AA,"")[[1]] == ensembl_AA_split) == TRUE)

# txdb cds mismatch gtf cds (stp codon)
# nonsense mediated decay; cds 5' incomplete

#-------------------------------------------------------------------------------#
# ALG3-203

ensembl_NT<-NA
ensembl_AA<-"XGLRKRGRSGSAAQAEGLCKQWLQRAWQERRLLLREPRYTLLVAACLCLAEVGITFWVIHRVACTQLVSCTSLWGCTMPPAEALTSAWPRTSLLCSTWLPCCLSS"
ensembl_AA_split<-strsplit(ensembl_AA,"")[[1]]

cds<-cdsBytranscript_gtf[["50358"]]
seq<-getSeq(Bsgenome, cds)
seq<-paste(as.character(seq[1:length(seq)]), collapse="")
trans<-translate(DNAStringSet(seq),no.init.codon = TRUE)
AA<-as.character(trans)
AA_match<- all((strsplit(AA,"")[[1]] == ensembl_AA_split) == TRUE)

# txdb cds mismatch gtf cds (cds2 start shifted right 1 codon)
# nonsense mediated decay; cds 5' incomplete

#-------------------------------------------------------------------------------#
# PRC1-217

ensembl_NT<-NA
ensembl_AA<-"XSIRPIFGGTVYHSPVSRLPPSGSKPVAASTCSGKKTPRTGRHGANKENLELNGSILSEGSVPL"
ensembl_AA_split<-strsplit(ensembl_AA,"")[[1]]

cds<-cdsBytranscript_gtf[["165035"]]
seq<-getSeq(Bsgenome, cds)
seq<-paste(as.character(seq[1:length(seq)]), collapse="")
trans<-translate(DNAStringSet(seq),no.init.codon = TRUE)
AA<-as.character(trans)
AA_match<- all((strsplit(AA,"")[[1]] == ensembl_AA_split) == TRUE)

# txdb cds mismatch gtf cds (cds2 start shifted right 1 codon)
# cds 5' incomplete

#-------------------------------------------------------------------------------#
# FCGR2A-205

ensembl_NT<-NA
ensembl_AA<-"AAPPKAVLKLEPPWINVLQEDSVTLTCQGARSPESDSIQWFHNGNLIPTHTQPSYRFKANNNDSGEYTCQTGQTSLSDPVHLTVLSEWLVLQTPHLEFQEGETIMLRCHSWKDKPLVKVTFFQNGKSQKFSHLDPTFSIPQANHSHSGDYHCTGNIGYTLFSSKPVTITVQAWAALHQWGSLWLWSLRLL*"
ensembl_AA_split<-strsplit(ensembl_AA,"")[[1]]

cds<-cdsBytranscript_gtf[["7376"]]
seq<-getSeq(Bsgenome, cds)
seq<-paste(as.character(seq[1:length(seq)]), collapse="")
trans<-translate(DNAStringSet(seq),no.init.codon = TRUE)
AA<-as.character(trans)
AA_match<- all((strsplit(AA,"")[[1]] == ensembl_AA_split) == TRUE)

# txdb cds mismatch gtf cds (stop codon)
# nonsense mediated decay; cds 5' incomplete

#-------------------------------------------------------------------------------#
# MPZL1-202

ensembl_NT<-NA
ensembl_AA<-"PFSSVTAGVSALEVYTPKEIFVANGTQGKLTCKFKSTSTTGGLTSVSWSFQPEGADTTVSFFHYSQGQVYLGNYPPFKDRISWAGDLDKKDASINIENMQFIHNGTYICDVKNPPDIVVQPGHIRLYVVEKENLPVFPVWVVVGIVTAVVLGLTLLISMILAVLYRRKNSKRDYTGAQSYMHS*"
ensembl_AA_split<-strsplit(ensembl_AA,"")[[1]]

cds<-cdsBytranscript_gtf[["7646"]]
seq<-getSeq(Bsgenome, cds)
seq<-paste(as.character(seq[1:length(seq)]), collapse="")
trans<-translate(DNAStringSet(seq),no.init.codon = TRUE)
AA<-as.character(trans)
AA_match<- all((strsplit(AA,"")[[1]] == ensembl_AA_split) == TRUE)

# txdb cds mismatch gtf cds (stop codon)
# cds 5' incomplete

#-------------------------------------------------------------------------------#
# RABGAP1L-215

ensembl_NT<-NA
ensembl_AA<-"SSWSMTFEERENRRLQEASMRLEQENDDLAHELVTSKIALRNDLDQAEDKADVLNKELLLTKQRLVETEEEKRKQEEETAQENHMCN*"
ensembl_AA_split<-strsplit(ensembl_AA,"")[[1]]

cds<-cdsBytranscript_gtf[["7902"]]
seq<-getSeq(Bsgenome, cds)
seq<-paste(as.character(seq[1:length(seq)]), collapse="")
trans<-translate(DNAStringSet(seq),no.init.codon = TRUE)
AA<-as.character(trans)
AA_match<- all((strsplit(AA,"")[[1]] == ensembl_AA_split) == TRUE)

# txdb cds mismatch gtf cds (stop codon)
# nonsense mediated decay; cds 5' incomplete

#-------------------------------------------------------------------------------#
# PPP1R12B-205

ensembl_NT<-NA
ensembl_AA<-"PQTIAPSTYVSTYLKRSIKM*"
ensembl_AA_split<-strsplit(ensembl_AA,"")[[1]]

cds<-cdsBytranscript_gtf[["8611"]]
seq<-getSeq(Bsgenome, cds)
seq<-paste(as.character(seq[1:length(seq)]), collapse="")
trans<-translate(DNAStringSet(seq),no.init.codon = TRUE)
AA<-as.character(trans)
AA_match<- all((strsplit(AA,"")[[1]] == ensembl_AA_split) == TRUE)

#-------------------------------------------------------------------------------#
# FBXO44-208

ensembl_NT<-NA
ensembl_AA<-"QVPGGGPQGRRVLGGADGYHTAGHRGQGLVRSQARLRVQVPAVRSAPVVRARASGDLPARPGDHPAEERCQVEGGLPHILQLPARRPLHLVSARRRGHSLLGRLVRPEGHQQQHHHRAPAALTPPEPPSAEP*"
ensembl_AA_split<-strsplit(ensembl_AA,"")[[1]]

cds<-cdsBytranscript_gtf[["830"]]
seq<-getSeq(Bsgenome, cds)
seq<-paste(as.character(seq[1:length(seq)]), collapse="")
trans<-translate(DNAStringSet(seq),no.init.codon = TRUE)
AA<-as.character(trans)
AA_match<- all((strsplit(AA,"")[[1]] == ensembl_AA_split) == TRUE)

#-------------------------------------------------------------------------------#
# SCP2-209

ensembl_NT<-NA
ensembl_AA<-"XALQHNLGIGGAVVVTLYKMGFPEAASSFRTHQIEAVPTSSASDGFKANLVFKEIEKKLEEEGEQFVKKIGGIFAFKVKDGPGGKEATWVVDVKNGKGSVLPNSDKKADCTITMADSDFLALMTGKMNPQSVSMMELISC"
ensembl_AA_split<-strsplit(ensembl_AA,"")[[1]]

cds<-cdsBytranscript_gtf[["3524"]]
seq<-getSeq(Bsgenome, cds)
seq<-paste(as.character(seq[1:length(seq)]), collapse="")
trans<-translate(DNAStringSet(seq),no.init.codon = TRUE)
AA<-as.character(trans)
AA_match<- all((strsplit(AA,"")[[1]] == ensembl_AA_split) == TRUE)

# txdb cds mismatch gtf cds (stop codon)
# cds 5' incomplete

#-------------------------------------------------------------------------------#
# AKR7A2-203

ensembl_NT<-NA
ensembl_AA<-"GTPVEETLHACQRLHQEGKFVELGLSNYASWEVAEICTLCKSNGWILPTVYQGAC*"
ensembl_AA_split<-strsplit(ensembl_AA,"")[[1]]

cds<-cdsBytranscript_gtf[["11694"]]
seq<-getSeq(Bsgenome, cds)
seq<-paste(as.character(seq[1:length(seq)]), collapse="")
trans<-translate(DNAStringSet(seq),no.init.codon = TRUE)
AA<-as.character(trans)
AA_match<- all((strsplit(AA,"")[[1]] == ensembl_AA_split) == TRUE)

