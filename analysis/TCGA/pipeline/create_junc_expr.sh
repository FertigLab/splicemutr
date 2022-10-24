# job submission params
#!/bin/sh
#$ -N create_junc_expr
#$ -S /bin/sh
#$ -l mem_free=30G,h_vmem=30G
#$ -o /users/tpalmer/TCGA/create_junc_expr/junc_expr.o
#$ -e /users/tpalmer/TCGA/create_junc_expr/junc_expr.e
#$ -M tpalme15@jhmi.edu
#$ -t 2-15 -tc 15

echo $(date)

module load conda
source activate /users/tpalmer/miniconda3/envs/R-4.0.2

TCGA_ROOT_DIR=/dcs04/fertig/data/theron/splicemutr_TCGA
TCGA_CANCER_FILE=/dcs04/fertig/data/theron/splicemutr_TCGA/cancer_dirs.txt
TCGA_CANCER=$(sed -n ${SGE_TASK_ID}p $TCGA_CANCER_FILE)
#TCGA_CANCER=$(sed -n 1p $TCGA_CANCER_FILE)
OUT_DIR=$TCGA_ROOT_DIR/$TCGA_CANCER/junction_counts

if [[ ! -d $OUT_DIR ]]
then
        mkdir $OUT_DIR
fi

JUNC_RSE_FILE=$TCGA_ROOT_DIR/$TCGA_CANCER/junc_rse.rds
SPLICE_DAT_FILE=$TCGA_ROOT_DIR/$TCGA_CANCER/combine_splicemutr_out/data_splicemutr_all_pep.rds

SCRIPT_DIR=/users/tpalmer/splicemute/scripts

$SCRIPT_DIR/create_junc_expr_TCGA.R -j $JUNC_RSE_FILE -s $SPLICE_DAT_FILE -o $OUT_DIR

echo $(date)
