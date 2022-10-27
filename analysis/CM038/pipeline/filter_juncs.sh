# job submission params
#!/bin/sh
#$ -N prepare_leafcutter
#$ -S /bin/sh
#$ -l mem_free=10G,h_vmem=12G
#$ -o /users/tpalmer/valsamo/leafcutter_prep/junc_filt.o 
#$ -e /users/tpalmer/valsamo/leafcutter_prep/junc_filt.e
#$ -M tpalme15@jhmi.edu

module load sharedapps
module load conda

source activate /users/tpalmer/miniconda3/envs/R-4.0.2

SJ_FILES=/dcs04/fertig/data/theron/share/bams/SJ_files.txt
OUT_DIR=/dcs04/fertig/data/theron/share/bams
SCRIPT_DIR=/users/tpalmer/splicemute/scripts

for VAR in {1..114}
do

STAR_JUNCFILE=$(sed -n ${VAR}p $SJ_FILES)

$SCRIPT_DIR/filter_juncs.R -o $OUT_DIR -s $STAR_JUNCFILE

done
