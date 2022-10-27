# job submission params
#!/bin/sh
#$ -N combine_splicemutr
#$ -S /bin/sh
#$ -l mem_free=2G,h_vmem=2G
#$ -o /users/tpalmer/valsamo/analysis/calc_coding_potential/cod_pot_calc.o
#$ -e /users/tpalmer/valsamo/analysis/calc_coding_potential/cod_pot_calc.e
#$ -M tpalme15@jhmi.edu
#$ -m ea
#$ -t 1-14 -tc 14

echo $(date)

module load conda
source activate /users/tpalmer/miniconda3/envs/R-4.0.2

OUT=/users/tpalmer/valsamo/analysis/form_transcripts/formed_transcripts
SPLICE_FILES=/users/tpalmer/valsamo/analysis/form_transcripts/formed_transcripts/filenames.txt
SPLICE_FILE=$(sed -n ${SGE_TASK_ID}p $SPLICE_FILES)
TRANSCRIPT_FILE=$(echo $SPLICE_FILE | sed s/'_data_splicemutr.rds'/'_sequences.fa'/g)
FUNCS=/users/tpalmer/splicemute/R/functions.R

SCRIPT_DIR=/users/tpalmer/splicemute/scripts

$SCRIPT_DIR/calc_coding_potential.R -o $OUT -s $SPLICE_FILE -t $TRANSCRIPT_FILE -f $FUNCS

echo $(date)
