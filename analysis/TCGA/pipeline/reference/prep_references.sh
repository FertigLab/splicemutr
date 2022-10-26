# job submission params
#!/bin/sh
#$ -N prepare_ref
#$ -S /bin/sh
#$ -q zappa
#$ -l mem_free=40G,h_vmem=45G
#$ -o /home/tpalme15/splicemutr_project/scripts/tcga/ref_prep.o
#$ -e /home/tpalme15/splicemutr_project/scripts/tcga/ref_prep.e
#$ -M tpalme15@jhmi.edu

echo $(date)

module load sharedapps conda
source activate /home/tpalme15/miniconda3/envs/R-4.0.2

SCRIPT_DIR=/home/tpalme15/splicemutr_project/leafcutter/scripts
LEAFVIZ_DIR=/home/tpalme15/splicemutr_project/leafcutter/leafviz
GTF=/home/tpalme15/splicemutr_project/reference_genomes/GTF_GFF/GENCODE/gencode.v26.chr_patch_hapl_scaff.annotation.gtf.gz
ANN_DIR=/home/tpalme15/splicemutr_project/TCGA_junctions/annotations

# creating the exons file necessary for differential analysis
$SCRIPT_DIR/gtf_to_exons.R $GTF $ANN_DIR/G026.exons.txt.gz

cd $ANN_DIR

# creating the necessary annotations for leafcutter visualizations
$LEAFVIZ_DIR/gtf2leafcutter.pl -o G026 $GTF
