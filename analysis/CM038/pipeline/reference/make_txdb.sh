# making the txdb object from the associated gtf

# job submission params
#!/bin/sh
#$ -N make_txdb
#$ -S /bin/sh
#$ -q zappa
#$ -l mem_free=10G,h_vmem=16G
#$ -o /home/tpalme15/splicemutr_project/scripts/running_splicemutr/generate_references/txdb_make.o
#$ -e /home/tpalme15/splicemutr_project/scripts/running_splicemutr/generate_references/txdb_make.e
#$ -M tpalme15@jhmi.edu

echo $(date)

module load sharedapps conda
source activate /home/tpalme15/miniconda3/envs/R-4.0.2

OUT=/home/tpalme15/splicemutr_project/reference_genomes/SEQUENCES/GENCODE/GRCh38_Ensembl99_sparseD3_sjdbOverhang99
GTF=$OUT/Homo_sapiens.GRCh38.99.gtf
SCRIPT_DIR=/home/tpalme15/splicemutr_project/splicemute/scripts


$SCRIPT_DIR/make_txdb.R -o $OUT/GRCh38_p13_txdb.sqlite -g $GTF
