#!/bin/bash
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH -a 1-5
#SBATCH -o slurm/slurm_files/experiment.out
#SBATCH -e slurm/slurm_files/experiment.err
#SBATCH -p common,statdept
#SBATCH --mem=15G # GB RAM
#SBATCH --mail-type=ALL
#SBATCH --mail-user=raphael.morsomme@duke.edu
hostname # print hostname
module load R/4.1.1-rhel8
Export R_LIBS_USER = ~/R/x86_64-pc-linux-gnu-library/4.1
R CMD BATCH --no-save R/mcmc_real_data.R