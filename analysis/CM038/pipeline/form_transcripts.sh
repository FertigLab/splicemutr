# job submission params
#!/bin/sh
#$ -N form_trans
#$ -S /bin/sh
#$ -l mem_free=5G,h_vmem=10G
#$ -o /users/tpalmer/valsamo/analysis/form_transcripts/trans_form.o
#$ -e /users/tpalmer/valsamo/analysis/form_transcripts/trans_form.e
#$ -M tpalme15@jhmi.edu
#$ -t 1-14 -tc 14

echo $(date)

module load conda
source activate /users/tpalmer/miniconda3/envs/R-4.0.2

OUT=/users/tpalmer/valsamo/analysis/form_transcripts/formed_transcripts
TXDB=/users/tpalmer/valsamo/reference/GRCh38_Ensembl99_sparseD3_sjdbOverhang99/GRCh38_p13_txdb.sqlite
INTRONS=/users/tpalmer/valsamo/analysis/form_transcripts/introns/intron_files.txt
INTRON_FILE=$(sed -n ${SGE_TASK_ID}p $INTRONS)
OUT_PREFIX=$OUT/$(echo $(basename $INTRON_FILE) | sed s/'.rds'/''/g)
SCRIPT_DIR=/users/tpalmer/splicemute/scripts

$SCRIPT_DIR/form_transcripts.R -o $OUT_PREFIX -t $TXDB -j $INTRON_FILE -b BSgenome.Hsapiens.GENCODE.GRCh38.p13

echo $(date)
