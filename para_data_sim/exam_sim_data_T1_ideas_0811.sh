#!/bin/bash
#SBATCH --job-name=pasimT1ideas
#SBATCH -o /home/guanwh/zhan7474/mayo/para_data_sim/exam_sim_data_T1_ideas_0811.out
#SBATCH -e /home/guanwh/zhan7474/mayo/para_data_sim/exam_sim_data_T1_ideas_0811.err
#SBATCH -p msismall
#SBATCH -c 80
#SBATCH --mem-per-cpu=2500
#SBATCH -t 96:00:00

source /common/software/install/migrated/anaconda/miniconda3_4.8.3-jupyter/etc/profile.d/conda.sh
conda activate r_env
Rscript /home/guanwh/zhan7474/mayo/para_data_sim/exam_sim_data_T1_ideas_0811.R