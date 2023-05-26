# The purpose of this file is to get the scripts for turning the STAR genome fasta file into the two bit format

# job submission params
#!/bin/sh
#$ -N build_2bit
#$ -S /bin/sh
#$ -q zappa
#$ -l mem_free=20G,h_vmem=26G
#$ -o /home/tpalme15/splicemutr_project/scripts/running_splicemutr/generate_references/bit_make.o
#$ -e /home/tpalme15/splicemutr_project/scripts/running_splicemutr/generate_references/bit_make.e
#$ -M tpalme15@jhmi.edu

echo $(date)

module load sharedapps conda
source activate /home/tpalme15/miniconda3/envs/R-4.0.2

OUT=/home/tpalme15/splicemutr_project/reference_genomes/SEQUENCES/GENCODE/GRCh38_Ensembl99_sparseD3_sjdbOverhang99
FASTA=$OUT/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

faToTwoBit $FASTA $OUT/Homo_sapiens.GRCh38.dna.primary_assembly.2bit
