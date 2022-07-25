# making the txdb reference for the mouse data

SCRIPT_DIR=/media/theron/My_Passport/splicemute/scripts
OUT_FILE=/media/theron/My_Passport/reference_genomes/GTF_GFF/ENSEMBL/Mus_musculus.GRCm38.102.chr.sqlite
GTF=/media/theron/My_Passport/reference_genomes/GTF_GFF/ENSEMBL/Mus_musculus.GRCm38.102.chr.gtf

$SCRIPT_DIR/make_txdb.R -o $OUT_FILE -g $GTF
