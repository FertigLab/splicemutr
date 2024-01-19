# build star index using grid node

# job submission params
#!/bin/sh
#$ -N buildSTARindex
#$ -S /bin/sh
#$ -l mem_free=45G,h_vmem=50G,h_fsize=50G
#$ -o /index_star_human.o
#$ -e /index_star_human.e
#$ -M tpalme15@jhmi.edu
#$ -m ea

echo $(date)

# load modules
#module load sharedapps
module load star

REF_DIR=/GRCh38.gencode.v39
OUT_DIR=$REF_DIR/STAR_GRCh38.gencode.v39

mkdir -p $OUT_DIR

cd $OUT_DIR

FASTA=$REF_DIR/GRCh38.primary_assembly.genome.fa
GFF3=$REF_DIR/gencode.v39.primary_assembly.annotation.gtf
OVERHANG=99

STAR --runMode genomeGenerate --sjdbGTFtagExonParentTranscript Parent --genomeDir $OUT_DIR --genomeFastaFiles $FASTA --sjdbGTFfile $GFF3 --sjdbOverhang $OVERHANG

echo $(date)
