# job submission params
#!/bin/sh
#$ -N tcga_gene_expr
#$ -S /bin/sh
#$ -q zappa
#$ -l mem_free=90,h_vmem=95G
#$ -o /home/tpalme15/splicemutr_project/scripts/tcga/gene_expression/expr_gene.o
#$ -e /home/tpalme15/splicemutr_project/scripts/tcga/gene_expression/expr_gene.e
#$ -M tpalme15@jhmi.edu
#$ -t 1-15 -tc 15

echo $(date)

module load sharedapps conda
source activate /home/tpalme15/miniconda3/envs/R-4.0.2

SCRIPT_DIR=/home/tpalme15/splicemutr_project/splicemute/scripts
DAT_FILES=/home/tpalme15/splicemutr_project/TCGA_junctions/filenames.txt
DAT_FILE=$(sed -n ${SGE_TASK_ID}p $DAT_FILES)
$SCRIPT_DIR/TCGA_gene_expr_per_sample.R -o $DAT_FILE

