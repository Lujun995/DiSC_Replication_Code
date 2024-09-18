#!/bin/bash
#SBATCH --job-name=covidsub
#SBATCH -o /home/guanwh/zhan7474/mayo/real_data_anal/COVID_PBMC_subset.out
#SBATCH -e /home/guanwh/zhan7474/mayo/real_data_anal/COVID_PBMC_subset.err
#SBATCH -p msismall
#SBATCH -c 1
#SBATCH --mem-per-cpu=200000
#SBATCH -t 24:00:00

source /common/software/install/migrated/anaconda/miniconda3_4.8.3-jupyter/etc/profile.d/conda.sh
conda activate r_env
Rscript /home/guanwh/zhan7474/mayo/real_data_anal/COVID_PBMC_subset.R
