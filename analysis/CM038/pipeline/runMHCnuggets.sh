# job submission params
#!/bin/sh
#$ -N mhcnuggets
#$ -S /bin/sh
#$ -l mem_free=10G,h_vmem=15G
#$ -o /users/tpalmer/valsamo/analysis/running_mhcnuggets/mhc_run.o
#$ -e /users/tpalmer/valsamo/analysis/running_mhcnuggets/mhc_run.e
#$ -M tpalme15@jhmi.edu
#$ -m ea
#$ -t 1-82 -tc 82

echo $(date)

module load conda
source activate /users/tpalmer/miniconda3/envs/R-4.0.2

TYPE="I"
INPUT_KMERS=/users/tpalmer/valsamo/analysis/process_peptides/process_peptides_out/peps_9.txt
MHC_ALLELE_FILE=/users/tpalmer/valsamo/analysis/running_mhcnuggets/class_1_alleles.txt
ALLELE=$(sed -n ${SGE_TASK_ID}p $MHC_ALLELE_FILE)
#ALLELE=$(sed -n 1p $MHC_ALLELE_FILE)
OUT_DIR=/users/tpalmer/valsamo/analysis/running_mhcnuggets/mhcnuggets_out
SCRIPT_DIR=/users/tpalmer/splicemute/inst

$SCRIPT_DIR/runMHCnuggets_ind.py -t $TYPE -k $INPUT_KMERS -m $ALLELE -o $OUT_DIR

echo $(date)
