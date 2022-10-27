# job submission params
#!/bin/sh
#$ -N recount3_init
#$ -S /bin/sh
#$ -l mem_free=30G,h_vmem=30G
#$ -o /users/tpalmer/TCGA/recount3_juncs/juncs_init.o
#$ -e /users/tpalmer/TCGA/recount3_juncs/juncs_init.e
#$ -M tpalme15@jhmi.edu

echo $(date)

module load conda
source activate /users/tpalmer/miniconda3/envs/R-4.0.2

OUT=/users/tpalmer/TCGA/recount3_juncs/TCGA_juncs
SCRIPT_DIR=/users/tpalmer/splicemute/scripts

$SCRIPT_DIR/recount3_tcga_juncs_init.R -o $OUT

echo $(date)
