# job submission params
#!/bin/sh
#$ -N create_junc_expr
#$ -S /bin/sh
#$ -l mem_free=10G,h_vmem=10G
#$ -o /save_introns/introns_split.o
#$ -e /save_introns/introns_split.e
#$ -M tpalme15@jhmi.edu
#$ -t 1-15 -tc 15

echo $(date)

module load conda
source activate /users/tpalmer/miniconda3/envs/R-4.0.2

TCGA_ROOT_DIR=/splicemutr_TCGA
TCGA_CANCER_FILE=/splicemutr_TCGA/cancer_dirs.txt # a file containing the one TCGA cancer subtype per line
TCGA_CANCER=$(sed -n ${SGE_TASK_ID}p $TCGA_CANCER_FILE)
INTRON_FILE=$TCGA_ROOT_DIR/$TCGA_CANCER/data.Rdata
OUT_DIR=$TCGA_ROOT_DIR/$TCGA_CANCER
SPLIT_NUM=5000

SCRIPT_DIR=/users/tpalmer/splicemute/scripts

$SCRIPT_DIR/split_introns.R -i $INTRON_FILE -o $OUT_DIR -s $SPLIT_NUM

echo $(date)