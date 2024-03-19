#!/bin/bash
#SBATCH --mail-type=end,fail
#SBATCH --job-name=run_STAR
#SBATCH --mem=50G
#SBATCH --output star_run.o
#SBATCH --error star_run.e
#SBATCH --cpus-per-task 10
#SBATCH --array=1-8

module load star

echo $(date)

FASTAFILES=$(pwd)/simulated_reads/fasta_files.txt
FASTAFILE_1=$(sed -n ${SLURM_ARRAY_TASK_ID}p $FASTAFILES)
FILENAME=$(basename ${FASTAFILE_1//'_1.fasta'/})
FASTAFILE_2=${FASTAFILE_1//'1.fasta'/'2.fasta'}
OUT_DIR=./simulated_reads/bams
mkdir -p $OUT_DIR
OUT_DIR=$OUT_DIR/$FILENAME

THREADS=10

REFERENCE=STAR_GRCh38.gencode.v39

STAR --genomeDir $GENOME_DIR --readFilesIn ${fastq1} ${fastq2} --twopassMode Basic --outSAMstrandField intronMotif --outFileNamePrefix $FILENAME --runThreadN 6 --readFilesCommand zcat --outSAMtype BAM Unsorted --quantMode GeneCounts
