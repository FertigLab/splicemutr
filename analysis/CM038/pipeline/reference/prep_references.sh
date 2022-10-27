# job submission params
#!/bin/sh
#$ -N prepare_leafcutter
#$ -S /bin/sh
#$ -q zappa
#$ -l mem_free=24G,h_vmem=30G
#$ -o /home/tpalme15/splicemutr_project/scripts/runningLeafcutter/run_10_20_2020/prep_references.o
#$ -e /home/tpalme15/splicemutr_project/scripts/runningLeafcutter/run_10_20_2020/prep_references.e
#$ -M tpalme15@jhmi.edu

module load sharedapps conda
source activate /home/tpalme15/miniconda3/envs/R-4.0.2

# creating the exons file necessary for differential analysis
/home/tpalme15/splicemutr_project/leafcutter/scripts/gtf_to_exons.R /home/tpalme15/splicemutr_project/reference_genomes/GTF_GFF/GENCODE/gencode.v35.annotation.gtf.gz /home/tpalme15/splicemutr_project/reference_genomes/leafcutter/gencode.v35.exons.txt.gz

cd /home/tpalme15/splicemutr_project/reference_genomes/leafcutter

# creating the necessary annotations for leafcutter visualizations
/home/tpalme15/splicemutr_project/leafcutter/leafviz/gtf2leafcutter.pl -o gencode.v35 /home/tpalme15/splicemutr_project/reference_genomes/GTF_GFF/GENCODE/gencode.v35.annotation.gtf.gz
