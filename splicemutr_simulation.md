---
title: "simulated_reads"
author: "Theron Palmer"
date: "2024-01-19"
output: html_document
---

Polyester simulations


```r
packages <- c("BiocManager", "BiocFileCache", "polyester", "Biostrings", "tidyverse", "dplyr")
lapply(packages, function(x){
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})
```

```
## Loading required package: BiocManager
```

```
## Loading required package: BiocFileCache
```

```
## Loading required package: dbplyr
```

```
## Loading required package: polyester
```

```
## Loading required package: Biostrings
```

```
## Loading required package: BiocGenerics
```

```
## 
## Attaching package: 'BiocGenerics'
```

```
## The following objects are masked from 'package:stats':
## 
##     IQR, mad, sd, var, xtabs
```

```
## The following objects are masked from 'package:base':
## 
##     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
##     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
##     get, grep, grepl, intersect, is.unsorted, lapply, Map, mapply,
##     match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
##     Position, rank, rbind, Reduce, rownames, sapply, setdiff, sort,
##     table, tapply, union, unique, unsplit, which.max, which.min
```

```
## Loading required package: S4Vectors
```

```
## Loading required package: stats4
```

```
## 
## Attaching package: 'S4Vectors'
```

```
## The following object is masked from 'package:utils':
## 
##     findMatches
```

```
## The following objects are masked from 'package:base':
## 
##     expand.grid, I, unname
```

```
## Loading required package: IRanges
```

```
## Loading required package: XVector
```

```
## Loading required package: GenomeInfoDb
```

```
## 
## Attaching package: 'Biostrings'
```

```
## The following object is masked from 'package:base':
## 
##     strsplit
```

```
## Loading required package: tidyverse
```

```
## Error: package or namespace load failed for 'tidyverse':
##  .onAttach failed in attachNamespace() for 'tidyverse', details:
##   call: NULL
##   error: package or namespace load failed for 'lubridate' in loadNamespace(j <- i[[1L]], c(lib.loc, .libPaths()), versionCheck = vI[[j]]):
##  there is no package called 'timechange'
```

```
## Error in contrib.url(repos, "source"): trying to use CRAN without setting a mirror
```

```r
set.seed(27)

# FASTA annotation
bfc <- BiocFileCache::BiocFileCache(ask = FALSE)
url <- "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/gencode.v39.pc_transcripts.fa.gz"
fasta_file <- BiocFileCache::bfcrpath(bfc, url)
fasta <- readDNAStringSet(fasta_file)

# generating transcript gene map
fasta_names <- names(fasta)
transcripts <- unname(vapply(fasta_names,function(fname){str_split(fname,"[|]")[[1]][1]},character(1)))
genes <- unname(vapply(fasta_names,function(fname){str_split(fname,"[|]")[[1]][2]},character(1)))
genes_trans_map <- data.frame(gene=genes,transcript=transcripts,index=seq(length(transcripts)))
gene_counts <- genes_trans_map %>% group_by(gene) %>% dplyr::count()

# counting transcripts per gene
gene_counts_filt <- gene_counts %>% dplyr::filter(n>1) 
gene_counts_filt <- gene_counts_filt[order(gene_counts_filt$n,decreasing=F),] 
gene_counts_filt <- gene_counts_filt[seq(200),]

genes_trans_map_target <- genes_trans_map %>% dplyr::filter(gene %in% gene_counts_filt$gene)

# 50 genes differentially expressed, 50 differentially spliced, 50 both, the rest neither

gene_counts_filt[seq(1,200),"type"] <- "diff_splice" # isoform usage varies across groups
rownames(gene_counts_filt) <- gene_counts_filt$gene
```

```
## Warning: Setting row names on a tibble is deprecated.
```

```r
genes_trans_map_target$type <- "diff_splice"

# ~20x coverage ----> reads per transcript = transcriptlength/readlength * 20
# here all transcripts will have ~equal FPKM
readspertx <- round(30 * width(fasta) / 100)

num_timepoints <- 2
countmat <- data.frame(matrix(0, nrow = length(fasta), ncol = num_timepoints))

a <- lapply(seq(nrow(gene_counts_filt)),function(row_val){
  gene_val <- rownames(gene_counts_filt)[row_val]
  genes_trans_map_small <- genes_trans_map_target %>% dplyr::filter(gene==gene_val)
  transcripts <- genes_trans_map_small$transcript
  indeces <- genes_trans_map_small$index
  type <- unique(genes_trans_map_small$type)
  up_side <- as.numeric(vapply(runif(1),function(rand_val){
    if (rand_val > 0.5) {
      return(2)
    } else {
      return(1)
    }
  },numeric(1)))
  down_side <- seq(2)[!(seq(2) %in% up_side)]
  up_indeces <- sample(indeces,floor(length(indeces)/2))
  down_indeces <- indeces[!(indeces %in% up_indeces)]
  if (type == "diff_splice") {
    total_reads_distr_1 <- sum(countmat[down_indeces, 1]) / length(up_indeces)
    total_reads_distr_2 <- sum(countmat[up_indeces, 1]) / length(down_indeces)

    countmat[up_indeces, up_side] <<- countmat[up_indeces, up_side]+4
    countmat[down_indeces, up_side] <<- 1
    countmat[down_indeces, down_side] <<- countmat[down_indeces, down_side] + 4
    countmat[up_indeces, down_side] <<- 1
  }
  return(c(indeces, type))
})

countmat_filt <- countmat[genes_trans_map_target$index, ]
countmat_labelled <- countmat_filt
countmat_labelled$gene <- genes_trans_map_target$gene
countmat_labelled$transcript <- genes_trans_map_target$transcript
countmat_filt <- as.matrix(countmat_filt)

readspertx_filt <- readspertx[genes_trans_map_target$index]

# subset the FASTA file to first 20 transcripts
small_fasta <- fasta[genes_trans_map_target$index]
writeXStringSet(small_fasta, "./gencode.v39.pc_transcripts_small.fa")

simulate_experiment("./gencode.v39.pc_transcripts_small.fa",
                    reads_per_transcript = readspertx_filt,
                    num_reps = c(4, 4),
                    fold_changes = countmat_filt,
                    outdir = './simulated_reads')
```


Forming groups file


```r
library(tidyverse)

sim_rep_file <- "./simulated_reads/sim_rep_info.txt"
sim_rep <- read.table(sim_rep_file, sep = "\t", header = T)

groups_file <-  sim_rep %>% select(rep_id, group)

write.table(groups_file,
            file = "./groups_file.txt",
            sep = "\t", col.names = FALSE, quote = FALSE, row.names = FALSE)
```

Forming genotypes_specific.txt file


HLA-A01-01,HLA-B02-01,HLA-C03-01  sample_1
HLA-A01-01,HLA-B02-01,HLA-C03-01  sample_2
HLA-A01-01,HLA-B02-01,HLA-C03-01  sample_3


```r
library(tidyverse)

sim_rep_file <- "./simulated_reads/sim_rep_info.txt"
sim_rep <- read.table(sim_rep_file, sep = "\t", header = TRUE)

hla <- paste(rep("HLA-A01-01", 6), collapse = ",")

sim_rep$HLA <- hla

genotypes_specific <- sim_rep %>% select(HLA, rep_id)
genotypes_specific$rep_id <- sprintf("%s.filt", genotypes_specific$rep_id)

write.table(
  genotypes_specific,
  sep = "\t",
  quote = FALSE,
  col.names = FALSE,
  row.names = FALSE,
  file = "./simulated_reads/genotypes_specific.txt"
)
```
