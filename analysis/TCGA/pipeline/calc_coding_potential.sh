# job submission params
#!/bin/sh
#$ -N calc_coding_potential
#$ -S /bin/sh
#$ -l mem_free=2G,h_vmem=2G
#$ -o /users/tpalmer/TCGA/calc_coding_potential/cod_pot_calc.o
#$ -e /users/tpalmer/TCGA/calc_coding_potential/cod_pot_calc.e
#$ -M tpalme15@jhmi.edu
#$ -m ea
# -t 1-14 -tc 14
#$ -t 15-16 -tc 2

echo $(date)

module load conda
source activate /users/tpalmer/miniconda3/envs/R-4.0.2

TCGA_ROOT_DIR=/dcs04/fertig/data/theron/splicemutr_TCGA
TCGA_CANCER_FILE=/dcs04/fertig/data/theron/splicemutr_TCGA/cancer_dirs.txt
TCGA_CANCER=$(sed -n ${SGE_TASK_ID}p $TCGA_CANCER_FILE)

OUT=$TCGA_ROOT_DIR/$TCGA_CANCER/formed_transcripts
SPLICE_FILES=$OUT/filenames.txt

START=1
END=$(wc -l $SPLICE_FILES | awk '{print $1}')
for (( VAL=$START; VAL<=$END; VAL++ ))
do

SPLICE_FILE=$(sed -n ${VAL}p $SPLICE_FILES)
TRANSCRIPT_FILE=$(echo $SPLICE_FILE | sed s/'_data_splicemutr.rds'/'_sequences.fa'/g)
FUNCS=/users/tpalmer/splicemute/R/functions.R

SCRIPT_DIR=/users/tpalmer/splicemute/scripts

$SCRIPT_DIR/calc_coding_potential.R -o $OUT -s $SPLICE_FILE -t $TRANSCRIPT_FILE -f $FUNCS

done

echo $(date)
