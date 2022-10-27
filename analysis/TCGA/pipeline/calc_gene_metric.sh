# job submission params
#!/bin/sh
#$ -N calc_gene_metric
#$ -S /bin/sh
#$ -l mem_free=5G,h_vmem=10G
#$ -o /users/tpalmer/TCGA/calc_gene_metric/gene_calc.o
#$ -e /users/tpalmer/TCGA/calc_gene_metric/gene_calc.e
#$ -M tpalme15@jhmi.edu
#$ -t 1-81 -tc 81

echo $(date)

module load conda
source activate /users/tpalmer/miniconda3/envs/R-4.0.2

JUNC_EXPR_FILES=/dcs04/fertig/data/theron/splicemutr_TCGA/junc_expr_vst.txt
JUNC_EXPR_FILE=$(sed -n ${SGE_TASK_ID}p $JUNC_EXPR_FILES)
#JUNC_EXPR_FILE=$(sed -n 1p $JUNC_EXPR_FILES)
FILE=$(echo $(basename $JUNC_EXPR_FILE) | sed "s/.rds//g")
CANCER_DIR_PRE=$(dirname $JUNC_EXPR_FILE)
CANCER_DIR=$(dirname $CANCER_DIR_PRE)

GENE_EXPR_FILE=$CANCER_DIR/gene_expression_vst.rds
SPLICE_DAT_FILE=$CANCER_DIR/combine_splicemutr_out_cp/data_splicemutr_all_pep.rds
KMER_COUNTS=$CANCER_DIR/kmer_counts/all_kmers.txt
OUT_DIR=$CANCER_DIR/GENE_METRIC_CP

if [[ ! -d $OUT_DIR ]]
then
        mkdir $OUT_DIR
fi

OUT_PREFIX=$OUT_DIR/$FILE
SCRIPT_DIR=/users/tpalmer/splicemute/scripts

$SCRIPT_DIR/calc_gene_metric_len_norm.R -g $GENE_EXPR_FILE -s $SPLICE_DAT_FILE -k $KMER_COUNTS -j $JUNC_EXPR_FILE -o $OUT_PREFIX

echo $(date)

