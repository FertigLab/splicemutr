# job submission params
#!/bin/sh
#$ -N extract_data
#$ -S /bin/sh
#$ -l mem_free=20G,h_vmem=20G
#$ -o /users/tpalmer/TCGA/extract_data/data_extract_large.o
#$ -e /users/tpalmer/TCGA/extract_data/data_extract_large.e
#$ -M tpalme15@jhmi.edu
#$ -m ea
#$ -t 1-17 -tc 17

echo $(date)

module load conda
source activate /users/tpalmer/miniconda3/envs/R-4.0.2

ALLELE_FILES=/users/tpalmer/TCGA/extract_data/faulty_alleles.txt
ALLELE=$(sed -n ${SGE_TASK_ID}p $ALLELE_FILES)
CANCER_DIR=$(dirname $(dirname $ALLELE))
OUT_DIR=$CANCER_DIR/process_bindaff_out
if [[ ! -d $OUT_DIR ]]
then
        mkdir $OUT_DIR
fi
BINDERS=$OUT_DIR/$(echo $(basename $ALLELE) | sed 's/.txt/_filt.txt/g')
awk -F "," '{ if ($2 <= 500) { print } }' $ALLELE > $BINDERS
PICKLE_DIR=$CANCER_DIR/process_peptides_out
KMER_LENGTH=9

SCRIPT_DIR=/users/tpalmer/splicemute/inst

$SCRIPT_DIR/process_bindaff.py -b $BINDERS -p $PICKLE_DIR -o $OUT_DIR -k $KMER_LENGTH

ALLELE_FILES=/users/tpalmer/TCGA/extract_data/faulty_alleles.txt
ALLELE=$(sed -n ${SGE_TASK_ID}p $ALLELE_FILES)
CANCER_DIR=$(dirname $(dirname $ALLELE))
PICKLE_DIR=$CANCER_DIR/process_bindaff_out
ALLELE_VAL=$(basename $ALLELE)
ALLELE_VAL=$(echo $ALLELE_VAL | sed 's/_peps_9.txt//g')

SCRIPT_DIR=/users/tpalmer/splicemute/inst

$SCRIPT_DIR/extract_data.py -a $ALLELE_VAL -p $PICKLE_DIR -b 9 -e 10

echo $(date)
