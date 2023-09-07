
The reference data used for the TCGA analysis is as follows. All SpliceMutr analysis requires a genomic fasta file as well as a GTF or GFF annotation. The references, preferably the orignals, used for alignment of the RNA-seq data MUST be the references used for SpliceMutr analysis. It is difficult to adapt unmatched references appropriately. 

# Reference file GENCODE links

GENCODE fasta: https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_26/GRCh38.primary_assembly.genome.fa.gz\
GENCODE GTF: https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_26/gencode.v26.chr_patch_hapl_scaff.annotation.gtf.gz\

# Scripts used to build references for LeafCutter and SpliceMutr

make_txdb.sh - creating the .txdb file necessary to access the genome annotation via SpliceMutr
prep_references.sh - creating the reference files necessary to run LeafCutter and LeafCutterMD

## Building the BSgenome reference used by SpliceMutr

convert_2bit.sh - converting the fasta reference to a .2bit\
Building the BSgenome package is done using the following directions: https://www.bioconductor.org/packages/devel/bioc/vignettes/BSgenome/inst/doc/BSgenomeForge.pdf

The "DESCRIPTION" file functions as the "seed" file described in the BSgenome documentation above. When you look at this file you will see conflicting reference names. This is only a naming issue. The correct GENCODE reference version (V26) was used to create the BSgenome and associated .txdb reference files for the TCGA analysis. Gencode V26 matches the splice-junction recount3 metadata information.
