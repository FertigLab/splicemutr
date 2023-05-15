# job submission params
#!/bin/bash
#$ -N leafcutter
#$ -S /bin/sh
#$ -N leafcutter
#$ -l mem_free=20G,h_vmem=25G
#$ -o /users/tpalmer/valsamo/analysis/runLeafcutter/leaf_run_12072021.o
#$ -e /users/tpalmer/valsamo/analysis/runLeafcutter/leaf_run_12072021.e
#$ -M tpalme15@jhmi.edu
#$ -t 1-117 -tc 15

echo $(date)

module load conda
source activate /users/tpalmer/miniconda3/envs/R-4.0.2

JUNC_DIR=/dcs04/fertig/data/theron/share/juncs
LEAF_SCRIPTS=/users/tpalmer/leafcutter/scripts
OUTLIER_FILES=/users/tpalmer/valsamo/analysis/leafcutter_prep/outlier_files/run_12072021/outlier_files.txt
JUNC_FILE=$(sed -n ${SGE_TASK_ID}p $OUTLIER_FILES)
JUNC_FILE=$(echo $JUNC_FILE)
SAMPLE=$(basename $JUNC_FILE | sed s/'.txt'/''/g)

cd $JUNC_DIR

echo "leafcutter_cluster_regtools"
python2 $LEAF_SCRIPTS/splicemutr_leafcutter_cluster_regtools.py -j $JUNC_FILE -o $SAMPLE -l 500000

echo "leafcutter_ds"
$LEAF_SCRIPTS/leafcutterMD.R --num_threads $NSLOTS -o $SAMPLE ${SAMPLE}_perind_numers.counts.gz

echo $(date)
