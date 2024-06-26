# Running the splicemutr simulation

Create the splicemutr environment first:
```
module load conda
conda env create -f ./splicemutr/envs/splicemutr_packages.yml #takes time
```

N.B. There is a known limitation that results in errors during the prep_ref.smk run unless addressed with a workaround. The workaround is mannuall installation of `optparse` and `BSgenome`.

```
conda activate splicemutr
R
install.packages("optparse")
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("BSgenome")
```

Create the leafcutter environment next:
```
module load conda
conda env create -f ./splicemutr/envs/leafcutter_package.yml #takes time
```

## prepare simulation data
splicemutr_simulation.Rmd runs the polyester simulations generating the fasta files that will be processed. 

run:
```
conda activate splicemutr
sbatch splicemutr/simulation/prepare-simulation.slurm
conda deactivate
```

## run STAR
Within the running_STAR folder there are two files, run_STAR.sh and buildSTARindex_human.sh. buildSTARindex_human.sh must be run before run_STAR.sh. You will need to download the appropriate annotation and genome fasta files, but the buildSTARindex_human.sh SLURM script, once customized for your specific annotation and genome fasta file, will generate the STAR index necessary to run STAR. The run_STAR.sh SLURM script, once modified to fit your exact file structure, will run STAR on your fasta files. You will need to install an anaconda version of STAR to run this script properly. That or ensure STAR is installed on the machine you will be running this pipeline on. 

run:
```
sbatch splicemutr/simulation/running_STAR/buildSTARindex_human.sh
sbatch splicemutr/simulation/running_STAR/run_STAR.sh
```

## prep references 

After the environment has been built, prep_references can be ran, for example:
```
conda activate splicemutr
snakemake --snakefile ./splicemutr/simulation/prep_references/prep_ref.smk \
          --configfile ./splicemutr/simulation/prep_references/config.yaml \
          --cores 4
```


Within the prep_references folder, you will find config, DESCRIPTION, and .smk files. The .smk file is a snakemake file that will prepare the referneces that will be used by SpliceMutr and LeafCutter for analysis. The DESCRIPTION file is a BSgenome file that will need to be modified to fit your BSgenome reference file. The specifications within the config file will need to be modified to fit your reference files that are to be downloaded and used by SpliceMutr. 

## LeafCutter

You will need to activate the leafcutter_package conda environment in order to run this portion of the SpliceMutr pipeline.

```
module load conda
conda activate leafcutter_package
```

After activating the LeafCutter package conda environment you can run the leafcutter snakemake using:

```
conda activate leafcutter_package
snakemake --snakefile ./splicemutr/simulation/running_leafcutter/running_leafcutter.smk \
          --configfile ./splicemutr/simulation/running_leafcutter/config.yaml \
          --cores 4
```

Within the running_leafcutter folder, you will find a config, a groups_file, and .smk file. The snakemake file can be run which will then run LeafCutter on your SJ.out.tab files output from the STAR alignment step. 

## SpliceMutr

You will need to activate the splicmutr conda environment in order to run this portion of the SpliceMutr pipeline.

```
module load conda
conda activate splicemutr
```

After activating the splicemutr package conda environment you can run the splicemutr snakemake using:

```
conda activate splicemutr
snakemake --snakefile ./splicemutr/simulation/running_splicemutr/run_splicemutr.smk \
          --configfile ./splicemutr/simulation/running_splicemutr/config.yaml \
          --cores 4
```

Within the running_splicemutr folder, you will find a config file and .smk file. The snakemake file will run the remainder of the SpliceMutr pipeline.

## Genotyping samples

The genotyping_samples directory is not necessary to run this SpliceMutr simulation. It is necessary if you are choosing to run your own samples though and specifically want to use arcasHLA to genotype your samples. The config and genotype_samples.smk file have been included within this directory for your reference.  
