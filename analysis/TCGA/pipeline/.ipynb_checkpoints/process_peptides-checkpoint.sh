# job submission params
#!/bin/sh
#$ -N process_peptides
#$ -S /bin/sh
#$ -l mem_free=5G,h_vmem=5G
#$ -o /process_peptides/peps_proc.o
#$ -e /process_peptides/peps_proc.e
#$ -M tpalme15@jhmi.edu
#$ -m ea
#$ -t 1-15 -tc 15

module load conda
source activate /users/tpalmer/miniconda3/envs/R-4.0.2

echo $(date)

SCRIPT_DIR=/splicemute/inst
TCGA_ROOT_DIR=/splicemutr_TCGA
TCGA_CANCER_FILE=/splicemutr_TCGA/cancer_dirs.txt
TCGA_CANCER=$(sed -n ${SGE_TASK_ID}p $TCGA_CANCER_FILE)
PEPTIDES=$TCGA_ROOT_DIR/$TCGA_CANCER/combine_splicemutr_out/proteins.txt
OUT_DIR=$TCGA_ROOT_DIR/$TCGA_CANCER/process_peptides_out
mkdir $OUT_DIR
KMER_LENGTH=9

$SCRIPT_DIR/process_peptides.py -p $PEPTIDES -o $OUT_DIR -k $KMER_LENGTH 

echo $(date)
