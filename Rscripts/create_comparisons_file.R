#!/usr/bin/env Rscript

#------------------------------------------------------------------------------#
# creating the simulation genotypes file

comparisons <- list()
comparisons[["data"]]$targets <- "sample_01"
comparisons[["data"]]$comparators <- c("sample_02",
                                       "sample_03",
                                       "sample_04",
                                       "sample_05",
                                       "sample_06",
                                       "sample_07",
                                       "sample_08")

#------------------------------------------------------------------------------#
# saving the genotypes file

saveRDS(comparisons,"comparisons.rds")
