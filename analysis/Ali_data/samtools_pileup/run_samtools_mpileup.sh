# running samtools mpileup locally using input BED locations

BED_FILE=/mnt/f/Ali_data/WES_Data_for_4T1MIS/MB01JHU503/MB01JHU503_000_analysis/MAF/mutation_coords.txt
BAM_1=/mnt/f/Ali_data/HTMBmice_RNAseq_Data/MB01JHU504/MB01JHU504_000_analysis/RSEM/alignments/C1_rsem_mm10.STAR.genome.sorted.bam
BAM_2=/mnt/f/Ali_data/HTMBmice_RNAseq_Data/MB01JHU504/MB01JHU504_000_analysis/RSEM/alignments/C2_rsem_mm10.STAR.genome.sorted.bam
BAM_3=/mnt/f/Ali_data/HTMBmice_RNAseq_Data/MB01JHU504/MB01JHU504_000_analysis/RSEM/alignments/C3_rsem_mm10.STAR.genome.sorted.bam
BAM_4=/mnt/f/Ali_data/HTMBmice_RNAseq_Data/MB01JHU504/MB01JHU504_000_analysis/RSEM/alignments/C5_rsem_mm10.STAR.genome.sorted.bam
BAM_5=/mnt/f/Ali_data/HTMBmice_RNAseq_Data/MB01JHU504/MB01JHU504_000_analysis/RSEM/alignments/C6_rsem_mm10.STAR.genome.sorted.bam
BAM_6=/mnt/f/Ali_data/HTMBmice_RNAseq_Data/MB01JHU504/MB01JHU504_000_analysis/RSEM/alignments/C7_rsem_mm10.STAR.genome.sorted.bam
BAM_7=/mnt/f/Ali_data/HTMBmice_RNAseq_Data/MB01JHU504/MB01JHU504_000_analysis/RSEM/alignments/P1_rsem_mm10.STAR.genome.sorted.bam
BAM_8=/mnt/f/Ali_data/HTMBmice_RNAseq_Data/MB01JHU504/MB01JHU504_000_analysis/RSEM/alignments/P2_rsem_mm10.STAR.genome.sorted.bam
BAM_9=/mnt/f/Ali_data/HTMBmice_RNAseq_Data/MB01JHU504/MB01JHU504_000_analysis/RSEM/alignments/P3_rsem_mm10.STAR.genome.sorted.bam
BAM_10=/mnt/f/Ali_data/HTMBmice_RNAseq_Data/MB01JHU504/MB01JHU504_000_analysis/RSEM/alignments/P5_rsem_mm10.STAR.genome.sorted.bam
BAM_11=/mnt/f/Ali_data/HTMBmice_RNAseq_Data/MB01JHU504/MB01JHU504_000_analysis/RSEM/alignments/P6_rsem_mm10.STAR.genome.sorted.bam
BAM_12=/mnt/f/Ali_data/HTMBmice_RNAseq_Data/MB01JHU504/MB01JHU504_000_analysis/RSEM/alignments/P7_rsem_mm10.STAR.genome.sorted.bam
BAM_13=/mnt/f/Ali_data/HTMBmice_RNAseq_Data/MB01JHU504/MB01JHU504_000_analysis/RSEM/alignments/U1_rsem_mm10.STAR.genome.sorted.bam
BAM_14=/mnt/f/Ali_data/HTMBmice_RNAseq_Data/MB01JHU504/MB01JHU504_000_analysis/RSEM/alignments/U2_rsem_mm10.STAR.genome.sorted.bam
BAM_15=/mnt/f/Ali_data/HTMBmice_RNAseq_Data/MB01JHU504/MB01JHU504_000_analysis/RSEM/alignments/U3_rsem_mm10.STAR.genome.sorted.bam
OUT_FILE=/mnt/f/Ali_data/HTMBmice_RNAseq_Data/MB01JHU504/MB01JHU504_000_analysis/RSEM/alignments/all_samples.pileup

samtools mpileup -l $BED_FILE -o $OUT_FILE -A -B $BAM_1 $BAM_2 $BAM_3 $BAM_4 $BAM_5 $BAM_6 $BAM_7 $BAM_8 $BAM_9 $BAM_10 $BAM_11 $BAM_12 $BAM_13 $BAM_14 $BAM_15
