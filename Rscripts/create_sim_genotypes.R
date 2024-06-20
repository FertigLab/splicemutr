#!/usr/bin/env Rscript

#------------------------------------------------------------------------------#
# creating the simulation genotypes file

genotypes <- list()
hlas <- c("HLA-A01-01")
names(hlas) <- "A*01:01"
genotypes["sample_01.filt"] <- hlas
genotypes["sample_02.filt"] <- hlas
genotypes["sample_03.filt"] <- hlas
genotypes["sample_04.filt"] <- hlas
genotypes["sample_05.filt"] <- hlas
genotypes["sample_06.filt"] <- hlas
genotypes["sample_07.filt"] <- hlas
genotypes["sample_08.filt"] <- hlas

#------------------------------------------------------------------------------#
# saving the genotypes file

saveRDS(genotypes,file="./genotypes.rds")
