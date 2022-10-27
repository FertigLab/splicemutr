# Running leafcutter on junction files. The junction files are contained in a single location. The file "junc_files.txt" contains a list of the directories for the junction files. The function leafcutter_cluster_regtools.py runs on the junc files from "junc_files.txt." 


# job submission params
#!/bin/sh
#$ -N prepare_leafcutter
#$ -S /bin/sh
#$ -q zappa
#$ -l mem_free=30G,h_vmem=36G
#$ -o /home/tpalme15/splicemutr_project/scripts/runningLeafcutter/run_10_21_2020/runleaf.3.o
#$ -e /home/tpalme15/splicemutr_project/scripts/runningLeafcutter/run_10_21_2020/runleaf.3.e
#$ -M tpalme15@jhmi.edu

module load sharedapps conda
source activate /home/tpalme15/miniconda3/envs/R-4.0.2

cd /home/tpalme15/splicemutr_project/sorted_bams/leafcutter_10_21_2020

#echo "leafcutter_cluster_regtools"
# The script leafcutter_cluster_regtools.py runs on the junc files from "junc_files.txt."
#python2 /home/tpalme15/splicemutr_project/leafcutter/scripts/leafcutter_cluster_regtools.py -j /home/tpalme15/splicemutr_project/sorted_bams/junc_files.txt -o data -l 500000

#cd /home/tpalme15/splicemutr_project/sorted_bams/leafcutter_10_21_2020

#echo "leafcutter_ds"
# The script leafcutter_ds.R takes the output of leafcutter_cluster_regtools.py,reference exons, and a groups file and determines differentially used intron clusters.
#/home/tpalme15/splicemutr_project/leafcutter/scripts/leafcutter_ds.R --num_threads 4 --exon_file=/home/tpalme15/splicemutr_project/reference_genomes/leafcutter/gencode.v35.exons.txt.gz /home/tpalme15/splicemutr_project/sorted_bams/leafcutter_10_20_2020/data_perind_numers.counts.gz /home/tpalme15/splicemutr_project/sorted_bams/groups_file_noid.txt

echo "prepare_results"
# The script prepare_results.R prepares the results of the scripts leafcutter_cluster_regtools.py and prepare_resutls.R for visualization. The output gives a list of differentially used intron clusters and individual clusters as well as metadata for the types of introns analyzed. This output, specifically introns, is what is being used as input into transcript formation
/home/tpalme15/splicemutr_project/leafcutter/leafviz/prepare_results.R -o /home/tpalme15/splicemutr_project/sorted_bams/leafcutter_10_21_2020/data.Rdata -m /home/tpalme15/splicemutr_project/sorted_bams/groups_file_noid.txt /home/tpalme15/splicemutr_project/sorted_bams/leafcutter_10_21_2020/data_perind_numers.counts.gz /home/tpalme15/splicemutr_project/sorted_bams/leafcutter_10_21_2020/leafcutter_ds_cluster_significance.txt /home/tpalme15/splicemutr_project/sorted_bams/leafcutter_10_21_2020/leafcutter_ds_effect_sizes.txt gencode.v35

