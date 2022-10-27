# job submission params
#!/bin/sh
#$ -N create_junc_expression
#$ -S /bin/sh
#$ -l mem_free=20G,h_vmem=25G
#$ -o /users/tpalmer/valsamo/analysis/create_junc_expression/expr_junc.o
#$ -e /users/tpalmer/valsamo/analysis/create_junc_expression/expr_junc.e
#$ -M tpalme15@jhmi.edu
#$ -m ea 

echo $(date)

module load conda
source activate /users/tpalmer/miniconda3/envs/R-4.0.2

SCRIPT_DIR=/users/tpalmer/splicemute/scripts
JUNC_DIR=/dcs04/fertig/data/theron/share/juncs/inner_juncs
JUNC_FILES=/dcs04/fertig/data/theron/share/juncs/inner_juncs/filenames.txt
OUT_DIR=/users/tpalmer/valsamo/analysis/create_junc_expression/create_junc_expression_out

$SCRIPT_DIR/create_junc_expr.R -j $JUNC_DIR -f $JUNC_FILES -o $OUT_DIR

echo $(date)
