# job submission params
#!/bin/sh
#$ -N extract_data
#$ -S /bin/sh
#$ -l mem_free=5G,h_vmem=10G
#$ -o /users/tpalmer/valsamo/analysis/extract_data/data_extract.o
#$ -e /users/tpalmer/valsamo/analysis/extract_data/data_extract.e
#$ -M tpalme15@jhmi.edu
#$ -m ea
#$ -t 1-82 -tc 82

echo $(date)

module load conda
source activate /users/tpalmer/miniconda3/envs/R-4.0.2

ALLELE_FILES=/users/tpalmer/valsamo/analysis/running_mhcnuggets/class_1_alleles.txt
ALLELE=$(sed -n ${SGE_TASK_ID}p $ALLELE_FILES)
#ALLELE=$(sed -n 1p $ALLELE_FILES)
PICKLE_DIR=/users/tpalmer/valsamo/analysis/process_bindaff/process_bindaff_out

SCRIPT_DIR=/users/tpalmer/splicemute/inst

$SCRIPT_DIR/extract_data.py -a $ALLELE -p $PICKLE_DIR -b 9 -e 10

echo $(date)
