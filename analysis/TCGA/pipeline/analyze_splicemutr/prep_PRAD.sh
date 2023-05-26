# job submission params
#!/bin/sh
#$ -N PREP_splicemutr
#$ -S /bin/sh
#$ -l mem_free=5G,h_vmem=10G
#$ -o /users/tpalmer/TCGA/analyze_splicemutr/splice_analyze_PRAD.o
#$ -e /users/tpalmer/TCGA/analyze_splicemutr/splice_analyze_PRAD.e
#$ -M tpalme15@jhmi.edu

echo $(date)


CANCER=PRAD

ROOT_DIR=/dcs04/fertig/data/theron/splicemutr_TCGA
SPLICE_DAT=$ROOT_DIR/$CANCER/combine_splicemutr_out/data_splicemutr_all_pep.txt
SPLICE_DAT_FILE=$(sed 's/.txt/_nov_corr.txt/g' <<< "$SPLICE_DAT")
sed 's/novel annotated pair/novel_annotated_pair/g' $SPLICE_DAT > $SPLICE_DAT_FILE

OUT_DIR=$ROOT_DIR/$CANCER/analyze_splicemutr_out
if [[ ! -d $OUT_DIR ]]
then
        mkdir $OUT_DIR
fi
echo $(date)

