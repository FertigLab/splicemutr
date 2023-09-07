# job submission params

#!/bin/sh
#$ -N recount3_expression
#$ -S /bin/sh
#$ -l mem_free=30G,h_vmem=30G
#$ -o /recount3_juncs/gene_expression.o
#$ -e /TCGA/recount3_juncs/gene_expression.e
#$ -M tpalme15@jhmi.edu

echo $(date)

module load conda
source activate /miniconda3/envs/R-4.0.2

OUT=/recount3_juncs/TCGA_juncs
CANCER_TYPE="PAAD" # an example of a TCGA cancer type to input
SCRIPT_DIR=/splicemute/scripts

$SCRIPT_DIR/recount3_extract_gene_expression.R -o $OUT -c $CANCER_TYPE

echo $(date)