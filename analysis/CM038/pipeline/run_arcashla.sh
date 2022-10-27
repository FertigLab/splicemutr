#this script is used to run splicemutr on subsets of the intron data and save the corresponding peps output

# job submission params
#!/bin/sh
#$ -N arcashla_geno
#$ -S /bin/sh
#$ -l mem_free=10G,h_vmem=15G
#$ -o /users/tpalmer/valsamo/arcashla/arcas_geno.o
#$ -e /users/tpalmer/valsamo/arcashla/arcas_geno.e
#$ -M tpalme15@jhmi.edu
#$ -m e
# -t 1-117 -tc 30

echo $(date)

module load conda
source activate /users/tpalmer/miniconda3/envs/arcashla

GENOTYPES_DIR=/dcs04/fertig/data/theron/share/genotypes
FILENAMES_FILE=/dcs04/fertig/data/theron/share/bams/bamfiles.txt
#FILE=$(sed -n ${SGE_TASK_ID}p $FILENAMES_FILE)
FILE=$(sed -n 1p $FILENAMES_FILE)
FILE_BASE=$(basename $FILE)
FILE_DIR=$GENOTYPES_DIR/${FILE_BASE}_dir
mkdir $FILE_DIR

ARCAS_SCRIPTS=/users/tpalmer/arcasHLA/scripts

# sort bam file
#samtools sort -o ${FILE}.sorted $FILE

python $ARCAS_SCRIPTS/extract.py ${FILE} -o $FILE_DIR -v

cd $FILE_DIR
#python $ARCAS_SCRIPTS/reference.py --update
FASTQ1=$(ls *.extracted.1*)
FASTQ2=$(ls *.extracted.2*)
python $ARCAS_SCRIPTS/genotype.py $FASTQ1 $FASTQ2 -g A,B,C,DPA1,DPB1,DQA1,DQB1,DRA,DRB1 -o $FILE_DIR -v
