# running process_cptac_peptides.py locally

CPTAC_PEPTIDES=/mnt/f/TCGA_junctions/ext_dat/CPTAC/TCGA_Breast_BI_Proteome.peptides.tsv
GENOTYPES=/mnt/f/TCGA_junctions/TCGA_cancers/BRCA/JHPCE/BRCA_genotypes_specific.txt
OUT_DIR=/mnt/f/TCGA_junctions/ext_dat/CPTAC
SCRIPT_DIR=/mnt/f/splicemute/inst

$SCRIPT_DIR/process_CPTAC_peptides.py -c $CPTAC_PEPTIDES -g $GENOTYPES -o $OUT_DIR