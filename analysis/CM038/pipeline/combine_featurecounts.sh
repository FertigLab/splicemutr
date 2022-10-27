# script to align fastq files and then sort and index the resulting alignment files

# job submission params
#!/bin/bash
#$ -N featurecounts
#$ -S /bin/sh
#$ -l mem_free=10G,h_vmem=15G
#$ -o /users/tpalmer/valsamo/featurecounts/feature_comb.o
#$ -e /users/tpalmer/valsamo/featurecounts/feature_comb.e
#$ -M tpalme15@jhmi.edu
#$ -m ea

module load conda
source activate /users/tpalmer/miniconda3/envs/R-4.0.2
 
echo $(date)

# parse variables
FEATURECOUNTS_FILE=/users/tpalmer/valsamo/featurecounts_out/filenames.txt
SCRIPT_DIR=/users/tpalmer/splicemute/scripts

$SCRIPT_DIR/combine_featurecounts.R -f $FEATURECOUNTS_FILE

