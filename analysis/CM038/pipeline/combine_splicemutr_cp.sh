# job submission params
#!/bin/sh
#$ -N combine_splicemutr
#$ -S /bin/sh
#$ -l mem_free=5G,h_vmem=10G
#$ -o /users/tpalmer/valsamo/analysis/combine_splicemutr/splice_comb.o
#$ -e /users/tpalmer/valsamo/analysis/combine_splicemutr/splice_comb.e
#$ -M tpalme15@jhmi.edu

echo $(date)

module load conda
source activate /users/tpalmer/miniconda3/envs/R-4.0.2

OUT=/users/tpalmer/valsamo/analysis/combine_splicemutr/combine_splicemutr_out_cp
SPLICE_FILES=/users/tpalmer/valsamo/analysis/form_transcripts/formed_transcripts/filenames_cp.txt
SCRIPT_DIR=/users/tpalmer/splicemute/scripts

$SCRIPT_DIR/combine_splicemutr.R -o $OUT -s $SPLICE_FILES

echo $(date)
