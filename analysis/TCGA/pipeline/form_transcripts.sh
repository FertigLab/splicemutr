# job submission params
#!/bin/sh
#$ -N form_trans
#$ -S /bin/sh
#$ -l mem_free=5G,h_vmem=10G
#$ -o /users/tpalmer/TCGA/form_transcripts/trans_form.o
#$ -e /users/tpalmer/TCGA/form_transcripts/trans_form.e
#$ -M tpalme15@jhmi.edu
#$ -t 1-222 -tc 222

echo $(date)

module load conda
source activate /users/tpalmer/miniconda3/envs/R-4.0.2

#for VAR in {1..16}
#do
#	TCGA_ROOT_DIR=/dcs04/fertig/data/theron/splicemutr_TCGA
#	TCGA_CANCERS_FILE=/dcs04/fertig/data/theron/splicemutr_TCGA/cancer_dirs.txt
#	TCGA_CANCER=$(sed -n ${VAR}p $TCGA_CANCERS_FILE)
#	mkdir $TCGA_ROOT_DIR/$TCGA_CANCER/formed_transcripts
#done

INTRONS=/dcs04/fertig/data/theron/splicemutr_TCGA/intron_files.txt
INTRON_FILE=$(sed -n ${SGE_TASK_ID}p $INTRONS)
#INTRON_FILE=$(sed -n 1p $INTRONS)
TCGA_CANCER_DIR=$(dirname $(dirname $(sed -n ${SGE_TASK_ID}p $INTRONS)))
OUT=$TCGA_CANCER_DIR/formed_transcripts
TXDB=/users/tpalmer/TCGA/reference/recount3/G026_txdb.sqlite
OUT_PREFIX=$OUT/$(echo $(basename $INTRON_FILE) | sed s/'.rds'/''/g)
SCRIPT_DIR=/users/tpalmer/splicemute/scripts

$SCRIPT_DIR/form_transcripts.R -o $OUT_PREFIX -t $TXDB -j $INTRON_FILE -b BSgenome.Hsapiens.GENCODE.GRCh38.p10

echo $(date)
