# running create_ref_kmer.sh locally

PEP_FA=/media/theron/My_Passport/reference_genomes/SEQUENCES/GENCODE/gencode.v19.pc_translations.fa
SCRIPT_DIR=/media/theron/My_Passport/scripts/Daria_HNSCC
OUTPUT_FILE=/media/theron/My_Passport/head_and_neck_DARIA/data/splicemutr_05_26_2021/ref_peps.txt

$SCRIPT_DIR/create_ref_kmers.py -k 200000 $PEP_FA $OUTPUT_FILE
