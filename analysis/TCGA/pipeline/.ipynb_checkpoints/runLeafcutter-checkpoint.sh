# job submission params

#!/bin/bash
#$ -N leafcutter
#$ -S /bin/sh
#$ -l mem_free=25G,h_vmem=25G
#$ -o /users/tpalmer/HNSCC_TCGA/runLeafCutter/leaf_run.o
#$ -e /users/tpalmer/HNSCC_TCGA/runLeafCutter/leaf_run.e
#$ -M tpalme15@jhmi.edu

echo $(date)

module load conda
source activate /users/tpalmer/miniconda3/envs/R-4.0.2

JUNC_DIR=./splice_junctions
LEAF_SCRIPTS=./leafcutter/scripts
REF_DIR=./leafcutter_annotations
LEAFVIZ_DIR=./leafcutter/leafviz
GROUPS_FILE=./groups_file.txt

echo "leafcutter_cluster_regtools"
python2 $LEAF_SCRIPTS/splicemutr_leafcutter_cluster_regtools.py -j $JUNC_DIR/junc_file.txt -r $JUNC_DIR -o data -l 500000

echo "leafcutter_ds"
LEAF_SCRIPTS/leafcutter_ds.R --num_threads 1 --exon_file=$REF_DIR/G026.exons.txt.gz -o $JUNC_DIR/leafcutter_ds $JUNC_DIR/data_perind_numers.counts.gz $GROUPS_FILE

echo "prepare_results"
$LEAFVIZ_DIR/prepare_results.R -o $JUNC_DIR/data.Rdata -m $GROUPS_FILE $JUNC_DIR/data_perind_numers.counts.gz $JUNC_DIR/leafcutter_ds_cluster_significance.txt $JUNC_DIR/leafcutter_ds_effect_sizes.txt $REF_DIR/G026

echo $(date)
