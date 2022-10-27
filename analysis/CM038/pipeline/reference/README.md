Prebuilt STAR indeces were obtained from the following location: https://labshare.cshl.edu/shares/gingeraslab/www-data/dobin/STAR/STARgenomes/Human/GRCh38_Ensembl99_sparseD3_sjdbOverhang99/

# Building the references necessary for running LeafCutterMD and SpliceMutr

prep_references.sh - Preparing the LeafCutterMD reference files.

make_txdb.sh - Making the .txdb object necessary for using the GTF annotation file within SpliceMutr transcript formation.


## Building the BSgenome object

make_2bit.sh - An example of converting a fasta file to a .2bit file for generation of the BSgenome.

BSgenome_seed.dcf - The seed for the BSgenome object.

Building the BSgenome package is done using the following directions: https://www.bioconductor.org/packages/devel/bioc/vignettes/BSgenome/inst/doc/BSgenomeForge.pdf

