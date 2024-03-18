#!/bin/bash
#SBATCH --mail-type=end,fail
#SBATCH --job-name=index_STAR
#SBATCH --mem=50G
#SBATCH --output index_star.o
#SBATCH --error index_star.e

echo $(date)

# load modules
#module load sharedapps
module load star

REF_DIR=$(pwd)/GRCh38.gencode.v39
OUT_DIR=$REF_DIR/STAR_GRCh38.gencode.v39

mkdir -p $OUT_DIR

cd $REF_DIR

FASTA_URL=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/GRCh38.primary_assembly.genome.fa.gz
GTF_URL=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/gencode.v39.primary_assembly.annotation.gtf.gz

wget $FASTA_URL
wget $GTF_URL
gunzip *

FASTA=$REF_DIR/GRCh38.primary_assembly.genome.fa
GTF=$REF_DIR/gencode.v39.primary_assembly.annotation.gtf
OVERHANG=99

STAR --runMode genomeGenerate --sjdbGTFtagExonParentTranscript Parent --genomeDir $OUT_DIR --genomeFastaFiles $FASTA --sjdbGTFfile $GTF --sjdbOverhang $OVERHANG

echo $(date)
