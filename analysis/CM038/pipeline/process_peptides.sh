# job submission params
#!/bin/sh
#$ -N process_peptides
#$ -S /bin/sh
#$ -l mem_free=5G,h_vmem=5G
#$ -o /users/tpalmer/valsamo/analysis/process_peptides/peps_proc.o
#$ -e /users/tpalmer/valsamo/analysis/process_peptides/peps_proc.e
#$ -M tpalme15@jhmi.edu
#$ -m ea

module load conda
source activate /users/tpalmer/miniconda3/envs/R-4.0.2

SCRIPT_DIR=/users/tpalmer/splicemute/inst
PEPTIDES=/users/tpalmer/valsamo/analysis/combine_splicemutr/combine_splicemutr_out/proteins.txt
OUT_DIR=/users/tpalmer/valsamo/analysis/process_peptides/process_peptides_out
KMER_LENGTH=9

$SCRIPT_DIR/process_peptides.py -p $PEPTIDES -o $OUT_DIR -k $KMER_LENGTH 
