#### WARNING: this is the main file which generates the simulated data
library(polyester)
library(Biostrings)
library(refGenome)
library(org.Hs.eg.db)
library(annotate)
library(GenomicFeatures)




#loading isoform informations
load("/Users/bahman/Dropbox/SEVApaper/PaperSuppl/Data/JoeData/isoforms.rda")
#transfer to log 2.
log2isos <- log2(isos+1)
#hist(apply(log2isos[,normsamp],MARGIN=1,mean))
# calculating means

txdb <- makeTxDbFromGFF(
  file="/Users/bahman/Documents/Postdoc Research/EVA/SimulatedData/unc_hg19.gtf.txt",
  format = "gtf") 


## two normal samples
normsamp <- pheno[which(pheno[,"classes"]=="Normal"),"junctionSample"]
tumorsamp <- pheno[which(pheno[,"classes"]=="Tumor"),"junctionSample"]

readnumavenorm <- rowMeans(log2isos[,normsamp])
readnumaveall <- rowMeans(log2isos)

expressediso <- names(which(readnumaveall>3))

isos2genes <- read.delim(
  "/Users/bahman/Dropbox/SEVApaper/PaperSuppl/Scripts/Simulation/unc_knownToLocus.txt",
                         header=F,sep="\t")



isos2genesvect <- as.vector(isos2genes[,1])
names(isos2genesvect) <- isos2genes[,2]

isos2genesexpressed <- sapply(strsplit(x=isos2genesvect[expressediso],split = "[|]"),
                              function(x) x[1])
#### converting gene symbol to gene id
x <- org.Hs.egSYMBOL2EG
# Get the entrez gene identifiers that are mapped to a gene symbol
mapped_genesx <- mappedkeys(x)
# Convert to a list
xy <- unlist(as.list(x[mapped_genesx]))

xx <- org.Hs.egCHR
# Get the entrez gene identifiers that are mapped to a chromosome
mapped_genesxx <- mappedkeys(xx)
# Convert to a list
xxy <- unlist(as.list(xx[mapped_genesxx]))


xxx <- org.Hs.egSYMBOL
# Get the gene symbol that are mapped to an entrez gene identifiers
mapped_genesxxx <- mappedkeys(xxx)
# Convert to a list
xxxy <- as.list(xxx[mapped_genesxxx])

genesChr1 <- unlist(xxxy[names(which(xxy[xy]==1))])
###################
###finding isos


rownames(pheno)<-pheno[,"junctionSample"]
genesmultiso <- names(which(tapply(X=names(isos2genesexpressed),
                    INDEX=isos2genesexpressed, FUN=length ) >1 ))
                    

genesmultisoChr1 <- intersect(genesmultiso,genesChr1)
isos2genesmod <- isos2genesexpressed[which(isos2genesexpressed %in% genesmultisoChr1)]
genesisos <- tapply(X=names(isos2genesmod),INDEX=isos2genesmod,
                    simplify=F, FUN=function(x) readnumaveall[x] )
#### The isoforms must have standard deviation at least 1.
genesisossign <- genesisos[which(sapply(genesisos,sd)>1)]

#### Choose 4 batches of 150 random genes. Use permutation to 
EachSetN <- 150 #number of genes in each class (DS&DE, DS&nonDE, nonDE&DS, nonDE&nonDS)
set.seed(2)
simulatedgenes <- sample(genesisossign,4*EachSetN,replace=F) #all the genes

indexstart <- 1
indexend <- EachSetN
neutralgenes <- simulatedgenes[indexstart:indexend]

indexstart <- EachSetN+1
indexend <- 2*EachSetN
# myp <- min(c(mean(sapply(simulatedgenes[indexstart:indexend],max)<5),1))
# DEnonDS <- sapply(simulatedgenes[indexstart:indexend],
#                   FUN=function(x) x+
#                     ifelse(test=min(x)<5,yes=1,no=sample(c(1,-1),
#                     prob=c(1-myp,myp),size=1)))

DEnonDS <- sapply(simulatedgenes[indexstart:indexend],
                FUN=function(x) x+sample(c(1,-1),size=1))


#Only DE genes with 50% chance up or down regulated in tumors
indexstart <- 2*EachSetN+1
indexend <- 3*EachSetN
nonDEDS <- sapply(simulatedgenes[indexstart:indexend],
                  FUN=function(x){y<-x; 
                                  y[which.max(x)] <- min(x);
                                  y[which.min(x)] <- max(x);
                                  return(y)})#Only DS genes

indexstart <- 3*EachSetN+1
indexend <- 4*EachSetN
# myp <- min(c(mean(sapply(simulatedgenes[indexstart:indexend],min)<5),1))
# 
# DEDS <- sapply(simulatedgenes[indexstart:indexend],
#                   FUN=function(x){y<-x; 
#                                   y[which.max(x)] <- min(x);
#                                   y[which.min(x)] <- max(x);
#                                   return(y+ifelse(test=min(y)<5,yes=1,no=sample(c(1,-1),
#                                            prob=c(myp,1-myp),size=1)))})                                           
#                                            #Only DS genes


DEDS <- sapply(simulatedgenes[indexstart:indexend],
               FUN=function(x){y<-x; 
                               y[which.max(x)] <- min(x);
                               y[which.min(x)] <- max(x);
                               return(y+sample(c(1,-1),size=1))})                                           
#Only DS genes


####generating two normal samples

simulatedgenes.norm.vector <- unlist(simulatedgenes)#isoforms average expression as a vector
simulatedgenes.tum.vector <- c(unlist(neutralgenes),unlist(DEnonDS),unlist(nonDEDS),unlist(DEDS))
#isoform average expression
fold_changes <- cbind(simulatedgenes.norm.vector,simulatedgenes.tum.vector)
isoformfinalnames<-sapply(names(simulatedgenes.tum.vector),FUN=function(x) {y<-strsplit(x,"[.]"); return(paste(y[[1]][2],y[[1]][3],sep="."))})
rownames(fold_changes) <- isoformfinalnames
colnames(fold_changes) <- c("normal","tumor")


#reading transcript
myfasta_file <- '/Users/bahman/Documents/Postdoc Research/EVA/SimulatedData/mygenome/hg19_M_rCRS_ref.transcripts.fa.txt'
myfasta <- readDNAStringSet(myfasta_file)

small_myfasta <- '/Users/bahman/Documents/Postdoc Research/EVA/SimulatedData/small_transcripts.fa'
#small_myfasta <- '/Users/bahman/Dropbox/SEVApaper/PaperSuppl/Results/Simulation/small_transcripts.fa'
writeXStringSet(myfasta[rownames(fold_changes)], file=small_myfasta)


i <- 2
curfold <- paste('/Users/bahman/Documents/Postdoc Research/EVA/SimulatedData/SimulatedNormalVSDEDS_',i,'/',sep="")
dir.create(curfold)

save(list=ls(),file=paste(curfold,"groundtruth.rda",sep=""))

ptm <- proc.time()
simulate_experiment( fasta=small_myfasta,reads_per_transcript=2, #reads_per_transcript=readspertx, 
                     num_reps=c(25,25), fold_changes= (2^fold_changes-1),
                     outdir=curfold)
#outdir='/Users/bahman/Dropbox/SEVApaper/PaperSuppl/Results/Simulation/FirstNiormals_2')
print(proc.time()-ptm)




#simulate 




#for(i in 1:25){


#}
##### End simulation 

# simulate_experiment_countmat( fasta=small_myfasta,readmat= fold_changes[names(myfasta)[1:200],]/100, 
#                      outdir='/Users/bahman/Dropbox/SEVApaper/PaperSuppl/Results/Simulation/FirstNiormals')
# 
# 
# #Three classes: The way we go
# load("/Users/bahman/Dropbox/SEVApaper/PaperSuppl/Data/JoeData/pvalAlpha.rda.RData")
# fold_changes = matrix(c(4,4,rep(1,18),1,1,4,4,rep(1,16),1.5,1.5,4,4,rep(2,16)), nrow=20)
# head(fold_changes)
# 
# 
# # FASTA annotation
# fasta_file = system.file('extdata', 'chr22.fa', package='polyester')
# fasta = readDNAStringSet(fasta_file)
# 
# 
# 
# 
# 
# genetransmap <- read.csv(file='/Users/bahman/Documents/Postdoc Research/EVA/SimulatedData/mygenome/unc_knownToLocus.txt',
#                          sep="\t")
# 
# 
# # subset the FASTA file to first 20 transcripts
# small_fasta = fasta[1:20]
# writeXStringSet(small_fasta, 'chr22_small.fa')
# 
# # ~20x coverage ----> reads per transcript = transcriptlength/readlength * 20
# # here all transcripts will have ~equal FPKM
# readspertx = round(20 * width(small_fasta) / 100)
# 
# # simulation call:
# simulate_experiment('chr22_small.fa', reads_per_transcript=readspertx, 
#                     num_reps=c(10,10,3), fold_changes=fold_changes, outdir='/Users/bahman/simulated_reads')
# 
# 
# ### If we can generate the exact read-counts
# # set up transcript-by-timepoint matrix:
# num_timepoints = 12
# countmat = matrix(readspertx, nrow=length(small_fasta), ncol=num_timepoints)
# 
# # add spikes in expression at certain timepoints to certain transcripts:
# up_early = c(1,2) 
# up_late = c(3,4)
# countmat[up_early, 2] = 3*countmat[up_early, 2]
# countmat[up_early, 3] = round(1.5*countmat[up_early, 3])
# countmat[up_late, 10] = 6*countmat[up_late, 10]
# countmat[up_late, 11] = round(1.2*countmat[up_late, 11])
# 
# # simulate reads:
# simulate_experiment_countmat('chr22_small.fa', readmat=countmat, 
#                              outdir='timecourse_reads') 
# 
# #### simulation of experiment
# source("http://bioconductor.org/biocLite.R")
# biocLite("ballgown")
# 
# library(ballgown)
# data(bg)
# bg = subset(bg, "chr=='22'")
# 
# # load gtf file for annotation:
# gtfpath = system.file('extdata', 'bg.gtf.gz', package='polyester')
# gtf = subset(gffRead(gtfpath), seqname=='22')
# 
# # load chromosome sequence corresponding to gtf file (just for this example)
# system('wget https://www.dropbox.com/s/04i6msi9vu2snif/chr22seq.rda')
# load('chr22seq.rda')
# names(chr22seq) = '22'
# 
# # simulate reads based on this experiment's FPKMs
# simulate_experiment_empirical(bg, grouplabels=pData(bg)$group, gtf=gtf,
#                               seqpath=chr22seq, mean_rps=5000, outdir='empirical_reads', seed=1247)