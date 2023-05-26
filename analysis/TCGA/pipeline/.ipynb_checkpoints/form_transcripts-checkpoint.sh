# job submission params
#!/bin/sh
#$ -N form_trans
#$ -S /bin/sh
#$ -l mem_free=5G,h_vmem=10G
#$ -o /form_transcripts/trans_form.o
#$ -e /form_transcripts/trans_form.e
#$ -M tpalme15@jhmi.edu
#$ -t 1-222 -tc 222

echo $(date)

module load conda
source activate /users/tpalmer/miniconda3/envs/R-4.0.2

INTRONS=/splicemutr_TCGA/intron_files.txt # this file is a conglomeration of all intron files for all TCGA cancer subtype
INTRON_FILE=$(sed -n ${SGE_TASK_ID}p $INTRONS)
TCGA_CANCER_DIR=$(dirname $(dirname $(sed -n ${SGE_TASK_ID}p $INTRONS)))
OUT=$TCGA_CANCER_DIR/formed_transcripts
TXDB=/reference/recount3/G026_txdb.sqlite
OUT_PREFIX=$OUT/$(echo $(basename $INTRON_FILE) | sed s/'.rds'/''/g)
SCRIPT_DIR=/splicemute/scripts

$SCRIPT_DIR/form_transcripts.R -o $OUT_PREFIX -t $TXDB -j $INTRON_FILE -b BSgenome.Hsapiens.GENCODE.GRCh38.p10

echo $(date)
