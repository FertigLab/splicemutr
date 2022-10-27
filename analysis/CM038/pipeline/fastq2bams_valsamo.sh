# script to align fastq files and then sort and index the resulting alignment files

# job submission params
#!/bin/bash
#$ -N fastq2bams
#$ -S /bin/sh
#$ -l mem_free=25G,h_vmem=30G
#$ -o /users/tpalmer/valsamo/running_STAR/fastq2bams/bam_fastq.o
#$ -e /users/tpalmer/valsamo/running_STAR/fastq2bams/bam_fastq.e
#$ -M tpalme15@jhmi.edu
#$ -t 1-114 -tc 10

module load conda
source activate /users/tpalmer/miniconda3/envs/R-4.0.2

echo $(date)

# parse variables
NUM=$SGE_TASK_ID # number of the line from input file that has filename you want to use as input to STAR
#NUM=1
fastq1=`sed "${NUM}q;d" /dcs04/fertig/data/theron/share/filenames.txt`
fastq2="${fastq1/1.clipped/2.clipped}"
samplename=${fastq1%"_1.clipped.fastq.gz"}
samplename=${samplename##*/}
outprefix="/dcs04/fertig/data/theron/share/bams/${samplename}"

GENOME_DIR=/users/tpalmer/valsamo/GRCh38_Ensembl99_sparseD3_sjdbOverhang99

# align fastqs using STAR
STAR --genomeDir $GENOME_DIR --readFilesIn ${fastq1} ${fastq2} --twopassMode Basic --outSAMstrandField intronMotif --outFileNamePrefix $outprefix --runThreadN 6 --readFilesCommand zcat --outSAMtype BAM Unsorted 

# convenience variable for sam file created by STAR aligner
samfile="${outprefix}Aligned.out.sam"

# delete any reads where CIGAR and sequence length are inconsistent 
fixedsam=`echo $samfile | sed s/sam/fixed.sam/`
cat $samfile | awk '{if (length($10)==length($11)) print $0}' > $fixedsam

# convert sam file to bam file
samtools view -S -b $fixedsam > ${outprefix}.bam

# sort bam file
#samtools sort ${outprefix}.bam ${outprefix}.sorted

# index bam file
#samtools index ${outprefix}.sorted.bam

# rm sam/fixedsam/bam/ file now that sorted bam has been created from it
rm $samfile
rm $fixedsam
#rm ${outprefix}.bam
