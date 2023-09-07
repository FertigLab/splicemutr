# job submission params

#!/bin/sh
#$ -N recount3_form
#$ -S /bin/sh
#$ -l mem_free=30G,h_vmem=30G
#$ -o /recount3_juncs/juncs_form.o
#$ -e /TCGA/recount3_juncs/juncs_form.e
#$ -M tpalme15@jhmi.edu

echo $(date)

module load conda
source activate /miniconda3/envs/R-4.0.2

OUT=/recount3_juncs/TCGA_juncs
SCRIPT_DIR=/splicemute/scripts

$SCRIPT_DIR/recount3_tcga_juncs.R -o $OUT

echo $(date)
