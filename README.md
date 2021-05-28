splicemutr is a tool for evaluating how individual splice junctions affect coding transcripts from RNAseq data.

# Software USED

STAR: 2.7.3a
R: 4.0.2

# STAR Alignment of RNAseq data


STAR --genomeDir $GENOME_DIR --readFilesIn ${fastq1} ${fastq2} --twopassMode Basic --outSAMstrandField intronMotif --outFileNamePrefix $outprefix --runThreadN 6 --readFilesCommand zcat --outSAMtype BAM Unsorted

