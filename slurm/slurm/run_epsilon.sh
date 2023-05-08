#!/bin/bash
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=5
#SBATCH -a 1-12
#SBATCH -o slurm_files/experiment.out
#SBATCH -e slurm_files/experiment.err
#SBATCH -p common,statdept
#SBATCH --mail-type=ALL
#SBATCH --mail-user=raphael.morsomme@duke.edu
hostname # print hostname
module load R
Export R_LIBS_USER = ~/R/x86_64-pc-linux-gnu-library/4.1
R CMD BATCH --no-save mcmc_run_epsilon.R