#!/bin/bash
#SBATCH --mail-user=tpalme15@jhmi.edu
#SBATCH --mail-type=end,fail
#SBATCH --job-name=index_STAR
#SBATCH --mem=50G
#SBATCH --output index_star.o
#SBATCH --error index_star.e

echo $(date)

# load modules
#module load sharedapps
module load star

REF_DIR=/GRCh38.gencode.v39
OUT_DIR=$REF_DIR/STAR_GRCh38.gencode.v39

mkdir -p $OUT_DIR

cd $REF_DIR

FASTA=$REF_DIR/GRCh38.primary_assembly.genome.fa
GFF3=$REF_DIR/gencode.v39.primary_assembly.annotation.gtf
OVERHANG=99

STAR --runMode genomeGenerate --sjdbGTFtagExonParentTranscript Parent --genomeDir $OUT_DIR --genomeFastaFiles $FASTA --sjdbGTFfile $GFF3 --sjdbOverhang $OVERHANG

echo $(date)
