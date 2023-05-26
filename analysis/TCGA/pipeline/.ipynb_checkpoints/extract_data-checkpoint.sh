# job submission params

#!/bin/sh
#$ -N extract_data
#$ -S /bin/sh
#$ -l mem_free=5G,h_vmem=10G
#$ -o /extract_data/data_extract.o
#$ -e /extract_data/data_extract.e
#$ -M tpalme15@jhmi.edu
#$ -m ea
#$ -t 1-15 -tc 15

echo $(date)

module load conda
source activate /users/tpalmer/miniconda3/envs/R-4.0.2

TCGA_ROOT_DIR=/splicemutr_TCGA
TCGA_CANCER_FILE=/splicemutr_TCGA/cancer_dirs.txt
MHC_ALLELE_FILE=$TCGA_ROOT_DIR/$TCGA_CANCER/${TCGA_CANCER}_class1_alleles.txt # this is a file containing all unique class 1 HLA alleles extracted from The Immune Landscape of Cancer Optitype calls

LINES_IN_FILE=$(wc -l $MHC_ALLELE_FILE | awk '{print $1}')
START=1

for (( VAL=$START; VAL<=$LINES_IN_FILE); VAL++ ))
do
    ALLELE=$(sed -n ${VAL}p $MHC_ALLELE_FILE)
    PICKLE_DIR=/process_bindaff/process_bindaff_out

    SCRIPT_DIR=/users/tpalmer/splicemute/inst

    $SCRIPT_DIR/extract_data.py -a $ALLELE -p $PICKLE_DIR -b 9 -e 10
done

echo $(date)