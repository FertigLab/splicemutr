# job submission params
#!/bin/sh
#$ -N calc_gene_metric
#$ -S /bin/sh
#$ -l mem_free=5G,h_vmem=10G
#$ -o /users/tpalmer/valsamo/analysis/calc_gene_metric/updated_gene_calc_len_norm.o
#$ -e /users/tpalmer/valsamo/analysis/calc_gene_metric/updated_gene_calc_len_norm.e
#$ -M tpalme15@jhmi.edu
#$ -m ea
#$ -t 1-22 -tc 22

echo $(date)

module load conda
source activate /users/tpalmer/miniconda3/envs/R-4.0.2

SCRIPT_DIR=/users/tpalmer/splicemute/scripts
GENE_EXPRESSION=/users/tpalmer/valsamo/analysis/create_gene_expression/featurecounts_all_vst.rds
SPLICE_DAT_FILES=/users/tpalmer/valsamo/analysis/create_comparisons/create_comparisons_out_cp/filenames.txt
SPLICE_DAT_FILE=$(sed -n ${SGE_TASK_ID}p $SPLICE_DAT_FILES)
#SPLICE_DAT_FILE=$(sed -n 22p $SPLICE_DAT_FILES)
COMPARISON=$(echo $(echo $(basename $SPLICE_DAT_FILE) | sed 's/splice_dat_//g') | sed 's/.rds//g')
KMER_COUNTS=/users/tpalmer/valsamo/analysis/create_comparisons/create_comparisons_out_cp/kmers_specific_${COMPARISON}.rds
JUNC_EXPR_FILE=/users/tpalmer/valsamo/analysis/create_junc_expression/create_junc_expression_out/junc_expr_combined_vst.rds
OUT_PREFIX=/dcs04/fertig/data/theron/share/calc_gene_metric_out_cp/${COMPARISON}

$SCRIPT_DIR/calc_gene_metric_len_norm.R -g $GENE_EXPRESSION -s $SPLICE_DAT_FILE -k $KMER_COUNTS -j $JUNC_EXPR_FILE -o $OUT_PREFIX

echo $(date)
