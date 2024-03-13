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

FASTAFILES=/simulated_reads/fasta_files.txt
FASTAFILE_1=$(sed -n ${SLURM_ARRAY_TASK_ID}p $FASTAFILES)
FILENAME=$(basename ${FASTAFILE_1//'_1.fasta'/})
FASTAFILE_2=${FASTAFILE_1//'1.fasta'/'2.fasta'}
OUT_DIR=/simulated_reads/bams
mkdir -p $OUT_DIR
OUT_DIR=$OUT_DIR/$FILENAME

THREADS=10

REFERENCE=STAR_GRCh38.gencode.v39

STAR --readFilesIn $FASTAFILE_1 $FASTAFILE_2 --genomeDir $REFERENCE --runThreadN $THREADS --twopassMode Basic --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.1 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outFilterType BySJout --outFilterScoreMinOverLread 0.66 --outFilterMatchNminOverLread 0.66 --limitSjdbInsertNsj 1200000 --outFileNamePrefix $OUT_DIR --outSAMstrandField intronMotif --outFilterIntronMotifs None --alignSoftClipAtReferenceEnds Yes --quantMode GeneCounts --outSAMtype BAM Unsorted --outSAMunmapped Within --genomeLoad NoSharedMemory --outSJtype Standard
