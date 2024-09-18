#!/bin/bash
#SBATCH --job-name=covid
#SBATCH -o /home/guanwh/zhan7474/mayo/real_data_anal/COVID_PBMC_anal_0612.out
#SBATCH -e /home/guanwh/zhan7474/mayo/real_data_anal/COVID_PBMC_anal_0612.err
#SBATCH -p msismall
#SBATCH -c 40
#SBATCH --mem-per-cpu=4000
#SBATCH -t 48:00:00

source /common/software/install/migrated/anaconda/miniconda3_4.8.3-jupyter/etc/profile.d/conda.sh
conda activate r_env
Rscript /home/guanwh/zhan7474/mayo/real_data_anal/COVID_PBMC_anal_0612.R

