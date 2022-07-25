CPTAC_FILE=/mnt/f/TCGA_junctions/ext_dat/CPTAC/cptac_sample_peptides.txt
BRCA_DIR=/mnt/f/TCGA_junctions/ext_dat/CPTAC
KMERS_FILES=/mnt/f/TCGA_junctions/ext_dat/CPTAC/filenames.txt
OUT_DIR=$BRCA_DIR/CPTAC_analysis
if [ ! -d $OUT_DIR ]
then
        mkdir $OUT_DIR
fi
SCRIPT_DIR=/mnt/f/splicemute/inst

$SCRIPT_DIR/process_CPTAC_splicemutr_overlap.py -c $CPTAC_FILE -s $KMERS_FILES -o $OUT_DIR -r 0 #SGE_TSK_ID
