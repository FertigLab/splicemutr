# job submission params
#!/bin/sh
#$ -N process_bindaff
#$ -S /bin/sh
#$ -l mem_free=5G,h_vmem=10G
#$ -o /users/tpalmer/valsamo/analysis/process_bindaff/bindaff_proc.o
#$ -e /users/tpalmer/valsamo/analysis/process_bindaff/bindaff_proc.e
#$ -M tpalme15@jhmi.edu
#$ -m ea
#$ -t 1-82 -tc 82

echo $(date)

module load conda
source activate /users/tpalmer/miniconda3/envs/R-4.0.2

ALLELE_FILES=/users/tpalmer/valsamo/analysis/running_mhcnuggets/mhcnuggets_out/allele_files.txt
OUT_DIR=/users/tpalmer/valsamo/analysis/process_bindaff/process_bindaff_out
ALLELE=$(sed -n ${SGE_TASK_ID}p $ALLELE_FILES)
#ALLELE=$(sed -n 1p $ALLELE_FILES)
BINDERS=$OUT_DIR/$(echo $(basename $ALLELE) | sed 's/.txt/_filt.txt/g')
awk -F "," '{ if ($2 <= 500) { print } }' $ALLELE > $BINDERS
PICKLE_DIR=/users/tpalmer/valsamo/analysis/process_peptides/process_peptides_out
KMER_LENGTH=9

SCRIPT_DIR=/users/tpalmer/splicemute/inst

$SCRIPT_DIR/process_bindaff.py -b $BINDERS -p $PICKLE_DIR -o $OUT_DIR -k $KMER_LENGTH

echo $(date)
