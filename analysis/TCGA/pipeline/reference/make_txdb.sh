# job submission params
#!/bin/sh
#$ -N make_txdb
#$ -S /bin/sh
#$ -q zappa
#$ -l mem_free=30G,h_vmem=36G
#$ -o /home/tpalme15/splicemutr_project/scripts/tcga/form_transcripts/txdb_make.o
#$ -e /home/tpalme15/splicemutr_project/scripts/tcga/form_transcripts/txdb_make.e
#$ -M tpalme15@jhmi.edu

echo $(date)

module load sharedapps conda
source activate /home/tpalme15/miniconda3/envs/R-4.0.2

OUT=/home/tpalme15/splicemutr_project/reference_genomes/SEQUENCES/GENCODE/recount3/G026_txdb.sqlite
GTF=/home/tpalme15/splicemutr_project/reference_genomes/GTF_GFF/GENCODE/gencode.v26.chr_patch_hapl_scaff.annotation.gtf.gz
SCRIPT_DIR=/home/tpalme15/splicemutr_project/splicemute/scripts

$SCRIPT_DIR/make_txdb.R -o $OUT -g $GTF
