#!/bin/bash
#SBATCH -o slurm/slurm_files/experiment.out
#SBATCH -e slurm/slurm_files/experiment.err
#SBATCH -p common,statdept
#SBATCH --mem=15G # GB RAM
#SBATCH --mail-type=ALL
#SBATCH --mail-user=raphael.morsomme@duke.edu
hostname # print hostname
module load R/4.1.1-rhel8
Export R_LIBS_USER = ~/R/x86_64-pc-linux-gnu-library
R CMD BATCH --no-save R/mcmc_simulation.R