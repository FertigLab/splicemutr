# job submission params
#!/bin/sh
#$ -N create_comparisons
#$ -S /bin/sh
#$ -l mem_free=3G,h_vmem=10G
#$ -o /users/tpalmer/valsamo/analysis/create_comparisons/comp_create.o
#$ -e /users/tpalmer/valsamo/analysis/create_comparisons/comp_create.e
#$ -M tpalme15@jhmi.edu
#$ -m ea
#$ -t 1-24 -tc 24

echo $(date)

module load conda
source activate /users/tpalmer/miniconda3/envs/R-4.0.2

COMP_NUM=${SGE_TASK_ID}
#COMP_NUM=1
COMPS_JUNCS_FILE=/users/tpalmer/valsamo/analysis/create_comparisons/comparison_juncs.rds
SPLICE_DAT_FILE=/users/tpalmer/valsamo/analysis/combine_splicemutr/combine_splicemutr_out_cp/data_splicemutr_all_pep.rds
COMPARISONS_FILE=/users/tpalmer/valsamo/analysis/create_comparisons/comparisons.rds
JUNC_DIR=/users/tpalmer/valsamo/analysis/analyze_splicemutr/analyze_splicemutr_out
OUT_DIR=/users/tpalmer/valsamo/analysis/create_comparisons/create_comparisons_out_cp
LEAF_DIR=/dcs04/fertig/data/theron/share/juncs/run_12072021/PERIND_COUNTS
SCRIPT_DIR=/users/tpalmer/splicemute/scripts

$SCRIPT_DIR/valsamo_create_comparisons.R -c $COMPS_JUNCS_FILE -d $SPLICE_DAT_FILE -e $COMPARISONS_FILE -n $COMP_NUM -j $JUNC_DIR -o $OUT_DIR -l $LEAF_DIR

echo $(date)
