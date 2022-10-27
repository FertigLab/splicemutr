# job submission params
#!/bin/sh
#$ -N analyze_splicemutr
#$ -S /bin/sh
#$ -l mem_free=3G,h_vmem=10G
#$ -o /users/tpalmer/valsamo/analysis/analyze_splicemutr/splice_analyze.o
#$ -e /users/tpalmer/valsamo/analysis/analyze_splicemutr/splice_analyze.e
#$ -M tpalme15@jhmi.edu
#$ -m ea
#$ -t 1-114 -tc 114

echo $(date)

module load conda
source activate /users/tpalmer/miniconda3/envs/R-4.0.2

GENOTYPES=/users/tpalmer/valsamo/analysis/analyze_splicemutr/genotypes.rds
SUMMARY_DIR=/users/tpalmer/valsamo/analysis/process_bindaff/process_bindaff_out
SPLICE_DAT_FILE=/users/tpalmer/valsamo/analysis/combine_splicemutr/combine_splicemutr_out/data_splicemutr_all_pep.rds
COUNTS_FILES=/dcs04/fertig/data/theron/share/juncs/inner_juncs/filenames.txt
COUNTS_FILE=$(sed -n ${SGE_TASK_ID}p $COUNTS_FILES)
#COUNTS_FILE=$(sed -n 1p $COUNTS_FILES)
OUT_DIR=/users/tpalmer/valsamo/analysis/analyze_splicemutr/analyze_splicemutr_out
SCRIPT_DIR=/users/tpalmer/splicemute/scripts

$SCRIPT_DIR/valsamo_analyze_splicemutr.R -g $GENOTYPES -s $SUMMARY_DIR -d $SPLICE_DAT_FILE -c $COUNTS_FILE -o $OUT_DIR

echo $(date)
