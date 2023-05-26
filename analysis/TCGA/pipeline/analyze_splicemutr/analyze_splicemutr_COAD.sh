# job submission params
#!/bin/sh
#$ -N COAD_analyze_splicemutr
#$ -S /bin/sh
#$ -l mem_free=5G,h_vmem=10G
#$ -o /users/tpalmer/TCGA/analyze_splicemutr/splice_analyze_COAD.o
#$ -e /users/tpalmer/TCGA/analyze_splicemutr/splice_analyze_COAD.e
#$ -M tpalme15@jhmi.edu
#$ -t 1-546 -tc 100

echo $(date)

module load conda
source activate /users/tpalmer/miniconda3/envs/R-4.0.2

CANCER=COAD
echo $CANCER $SGE_TASK_ID

ROOT_DIR=/dcs04/fertig/data/theron/splicemutr_TCGA
GENOTYPES=$ROOT_DIR/$CANCER/${CANCER}_genotypes_specific.txt
SUMMARY_DIR=$ROOT_DIR/$CANCER/process_bindaff_out
SPLICE_DAT_FILE=$ROOT_DIR/$CANCER/combine_splicemutr_out/data_splicemutr_all_pep_nov_corr.txt

OUT_DIR=$ROOT_DIR/$CANCER/analyze_splicemutr_out
SUMMARY_TYPE='IC50'
SCRIPT_DIR=/users/tpalmer/splicemute/inst

$SCRIPT_DIR/analyze_splicemutr.py -g $GENOTYPES -s $SUMMARY_DIR -d $SPLICE_DAT_FILE -o $OUT_DIR -t $SUMMARY_TYPE -n $SGE_TASK_ID

echo $(date)
