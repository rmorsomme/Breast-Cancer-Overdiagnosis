#!/bin/bash
#SBATCH -o slurm_files/experiment.out
#SBATCH -e slurm_files/experiment.err
#SBATCH -p common,statdept
#SBATCH --mem=2G # 4 GB RAM
#SBATCH --mail-type=ALL
#SBATCH --mail-user=raphael.morsomme@duke.edu
hostname # print hostname
module load R
Export R_LIBS_USER = ~/R/x86_64-pc-linux-gnu-library/4.1
R CMD BATCH --no-save mcmc_run.R