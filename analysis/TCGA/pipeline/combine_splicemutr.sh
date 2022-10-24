# job submission params
#!/bin/sh
#$ -N combine_splicemutr
#$ -S /bin/sh
#$ -l mem_free=5G,h_vmem=10G
#$ -o /users/tpalmer/TCGA/combine_splicemutr/splice_comb_cp.o
#$ -e /users/tpalmer/TCGA/combine_splicemutr/splice_comb_cp.e
#$ -M tpalme15@jhmi.edu
#$ -t 1-16 -tc 16

echo $(date)

module load conda
source activate /users/tpalmer/miniconda3/envs/R-4.0.2

TCGA_ROOT_DIR=/dcs04/fertig/data/theron/splicemutr_TCGA
TCGA_CANCER_FILE=/dcs04/fertig/data/theron/splicemutr_TCGA/cancer_dirs.txt
TCGA_CANCER=$(sed -n ${SGE_TASK_ID}p $TCGA_CANCER_FILE)

SPLICE_FILES=$TCGA_ROOT_DIR/$TCGA_CANCER/formed_transcripts/filenames_cp.txt
OUT=$TCGA_ROOT_DIR/$TCGA_CANCER/combine_splicemutr_out_cp
mkdir $OUT
SCRIPT_DIR=/users/tpalmer/splicemute/scripts

$SCRIPT_DIR/combine_splicemutr.R -o $OUT -s $SPLICE_FILES

echo $(date)
