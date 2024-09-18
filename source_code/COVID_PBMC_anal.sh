#!/bin/bash
#SBATCH --job-name=covid
#SBATCH -o WORKDIRCOVID_PBMC_anal_DATE.out
#SBATCH -e WORKDIRCOVID_PBMC_anal_DATE.err
#SBATCH -p msismall
#SBATCH -c NCORES
#SBATCH --mem-per-cpu=MEM
#SBATCH -t 48:00:00

source /common/software/install/migrated/anaconda/miniconda3_4.8.3-jupyter/etc/profile.d/conda.sh
conda activate r_env
Rscript WORKDIRCOVID_PBMC_anal_DATE.R

