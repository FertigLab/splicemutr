SJ_files=/media/theron/One_Touch/KenOlive_CUMC_LCM_RNASeq_SJout/SJout/filenames.txt
TXDB=/media/theron/One_Touch/reference_genomes/SEQUENCES/GENCODE/GRCh38_Ensembl99_sparseD3_sjdbOverhang99/GRCh38_p13_txdb.sqlite
./find_CTLA4.R -o $PWD -t $TXDB -j $SJ_files
