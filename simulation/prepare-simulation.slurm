#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --job-name="splicemutr-prepare-simulation"
#SBATCH --mem=8GB
#SBATCH --output=prepare-simulation.out

module load R

Rscript -e 'knitr::knit("splicemutr_simulation.Rmd")'
