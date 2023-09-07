# job submission params

#!/bin/sh
#$ -N mhcnuggets
#$ -S /bin/sh
#$ -l mem_free=15G,h_vmem=15G
#$ -o /running_mhcnuggets/mhc_run.o
#$ -e /running_mhcnuggets/mhc_run.e
#$ -M tpalme15@jhmi.edu
#$ -m ea
#$ -t 1-16 -tc 16

echo $(date)

module load conda
source activate /users/tpalmer/miniconda3/envs/R-4.0.2

TCGA_ROOT_DIR=/splicemutr_TCGA
TCGA_CANCER_FILE=/splicemutr_TCGA/cancer_dirs.txt
TCGA_CANCER=$(sed -n ${SGE_TASK_ID}p $TCGA_CANCER_FILE)

TYPE="I"
INPUT_KMERS=$TCGA_ROOT_DIR/$TCGA_CANCER/process_peptides_out/peps_9.txt # an ouput from the previous step
MHC_ALLELE_FILE=$TCGA_ROOT_DIR/$TCGA_CANCER/${TCGA_CANCER}_class1_alleles.txt # this is a file containing all unique class 1 HLA alleles extracted from The Immune Landscape of Cancer Optitype calls 
OUT_DIR=$TCGA_ROOT_DIR/$TCGA_CANCER/mhcnuggets_out
mkdir $OUT_DIR
SCRIPT_DIR=/users/tpalmer/splicemute/inst

$SCRIPT_DIR/runMHCnuggets.py -t $TYPE -k $INPUT_KMERS -m $MHC_ALLELE_FILE -o $OUT_DIR

echo $(date)
