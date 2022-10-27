# script to align fastq files and then sort and index the resulting alignment files

# job submission params
#!/bin/bash
#$ -N featurecounts
#$ -S /bin/sh
#$ -l mem_free=1G,h_vmem=1G
#$ -o /users/tpalmer/valsamo/featurecounts/feature_run.o
#$ -e /users/tpalmer/valsamo/featurecounts/feature_run.e
#$ -M tpalme15@jhmi.edu
#$ -m ea
#$ -pe local 20
#$ -R y
#$ -t 1-117 -tc 1

module load featurecounts/2.0.0

echo $(date)

# parse variables
NUM=$SGE_TASK_ID # number of the line from input file that has filename you want to use as input to STAR
#NUM=1
BAM_FILE=/dcs04/fertig/data/theron/share/bams/bamfiles.txt
cd /dcs04/fertig/data/theron/share/bams
BAM=$(sed -n ${SGE_TASK_ID}p $BAM_FILE)

GTF_FILE=/users/tpalmer/valsamo/GRCh38_Ensembl99_sparseD3_sjdbOverhang99/Homo_sapiens.GRCh38.99.gtf
OUT=/users/tpalmer/valsamo/featurecounts_out/$(basename $BAM)_feature_counts.txt

featureCounts -F GTF -a $GTF_FILE -O -s 0 -M -T 20 --largestOverlap --minOverlap 8 -p -C --donotsort -o $OUT $BAM
