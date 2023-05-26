# job submission params
#!/bin/sh
#$ -N create_genotypes
#$ -S /bin/sh
#$ -l mem_free=5G,h_vmem=5G
#$ -o /users/tpalmer/TCGA/create_genotypes/create_genotypes.o
#$ -e /users/tpalmer/TCGA/create_genotypes/create_genotypes.e
#$ -M tpalme15@jhmi.edu
#$ -m ea
#$ -t 1-15 -tc 15

echo $(date)

module load conda
source activate /users/tpalmer/miniconda3/envs/R-4.0.2

DAT_FILES=/dcs04/fertig/data/theron/splicemutr_TCGA/cancer_dirs_filt.txt
DAT_FILE=/dcs04/fertig/data/theron/splicemutr_TCGA/$(sed -n ${SGE_TASK_ID}p $DAT_FILES)
#DAT_FILE=/dcs04/fertig/data/theron/splicemutr_TCGA/$(sed -n 1p $DAT_FILES)
CANCER=$(basename $DAT_FILE)
GENOTYPES_FILE=/users/tpalmer/TCGA/create_genotypes/OptiTypeCallsHLA_20171207.tsv

SCRIPT_DIR=/users/tpalmer/splicemute/scripts

$SCRIPT_DIR/TCGA_create_genotypes.R -m $DAT_FILE/junc_metadata.rds -g $GENOTYPES_FILE
