#!/bin/bash
#SBATCH --job-name=resimtype1
#SBATCH -o /home/guanwh/zhan7474/mayo/real_data_sim/type1_deseq_0814.out
#SBATCH -e /home/guanwh/zhan7474/mayo/real_data_sim/type1_deseq_0814.err
#SBATCH -p msismall
#SBATCH -c 80
#SBATCH --mem-per-cpu=2500
#SBATCH -t 60:00:00

source /common/software/install/migrated/anaconda/miniconda3_4.8.3-jupyter/etc/profile.d/conda.sh
conda activate r_env
Rscript /home/guanwh/zhan7474/mayo/real_data_sim/type1_deseq_0814.R