# Running the splicemutr simulation

## prepare simulation data
splicemutr_simulation.Rmd runs the polyester simulations generating the fasta files that will be processed. To prepare on the cluster, run the `prepare-simulation.slurm` script on cluster or `splicemutr_simulation.Rmd` locally.

## run STAR
Within the running_STAR folder there are two files, run_STAR.sh and buildSTARindex_human.sh. buildSTARindex_human.sh must be run before run_STAR.sh. You will need to download the appropriate annotation and genome fasta files, but the buildSTARindex_human.sh SLURM script, once customized for your specific annotation and genome fasta file, will generate the STAR index necessary to run STAR. The run_STAR.sh SLURM script, once modified to fit your exact file structure, will run STAR on your fasta files. You will need to install an anaconda version of STAR to run this script properly. That or ensure STAR is installed on the machine you will be running this pipeline on. 

Within the prep_references folder, you will find config, DESCRIPTION, and .smk files. The .smk file is a snakemake file that will prepare the referneces that will be used by SpliceMutr and LeafCutter for analysis. The DESCRIPTION file is a BSgenome file that will need to be modified to fit your BSgenome reference file. The specifications within the config file will need to be modified to fit your reference files that are to be downloaded and used by SpliceMutr. 

You will need to activate the leafcutter_pacakge.yml conda environment within the envs folder in order to run this portion of the SpliceMutr pipeline. Within the running_leafcutter folder, you will find a config, a groups_file, and .smk file. Once the config file is modified to fit your personal specification, the snakemake file can be run which will modify and then run LeafCutter on your SJ.out.tab files output from the STAR alignment step. 

You will need to activate the splicemutr_packages.yml conda environment within the envs folder in order to run this portion of the SpliceMutr pipeline. Within the running_splicemutr folder, you will find a config file and .smk file. Once the config file is modified to fit your personal file structure, the snakemake file can be run. This snakemake will run the SpliceMutr pipeline. 

The genotyping_samples directory is not necessary to run this SpliceMutr simulation. It is necessary if you are choosing to run your own samples though and specifically want to use arcasHLA to genotype your samples. The config and genotype_samples.smk file have been included within this directory for your reference.  
