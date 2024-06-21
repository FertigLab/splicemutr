#!/usr/bin/env Rscript

#------------------------------------------------------------------------------#
# creating the simulation genotypes file

genotypes <- list()
hlas <- c("HLA-A01-01")
names(hlas) <- "A*01:01"
genotypes["sample_01"] <- hlas
genotypes["sample_02"] <- hlas
genotypes["sample_03"] <- hlas
genotypes["sample_04"] <- hlas
genotypes["sample_05"] <- hlas
genotypes["sample_06"] <- hlas
genotypes["sample_07"] <- hlas
genotypes["sample_08"] <- hlas

#------------------------------------------------------------------------------#
# saving the genotypes file

saveRDS(genotypes,file="./genotypes.rds")
