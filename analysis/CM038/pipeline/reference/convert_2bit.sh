# job submission params
#!/bin/sh
#$ -N convert
#$ -S /bin/sh
#$ -q zappa
#$ -l mem_free=30G,h_vmem=36G
#$ -o /home/tpalme15/splicemutr_project/reference_genomes/SEQUENCES/GENCODE/recount3/2bit_convert.o
#$ -e /home/tpalme15/splicemutr_project/reference_genomes/SEQUENCES/GENCODE/recount3/2bit_convert.e
#$ -M tpalme15@jhmi.edu

echo $(date)

module load sharedapps conda
source activate /home/tpalme15/miniconda3/envs/R-4.0.2

GZIPPED_FILE=/home/tpalme15/splicemutr_project/reference_genomes/SEQUENCES/GENCODE/recount3/GRCh38.primary_assembly.genome.fa.gz

gunzip $GZIPPED_FILE

FASTA_FILE=/home/tpalme15/splicemutr_project/reference_genomes/SEQUENCES/GENCODE/recount3/GRCh38.primary_assembly.genome.fa
TWOBIT_FILE=/home/tpalme15/splicemutr_project/reference_genomes/SEQUENCES/GENCODE/recount3/GRCh38.primary_assembly.genome.2bit

faToTwoBit $FASTA_FILE $TWOBIT_FILE
