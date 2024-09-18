#!/bin/bash
#SBATCH --job-name=pasim3deseq
#SBATCH -o /home/guanwh/zhan7474/mayo/para_data_sim/exam_sim_data_3_deseq_0811.out
#SBATCH -e /home/guanwh/zhan7474/mayo/para_data_sim/exam_sim_data_3_deseq_0811.err
#SBATCH -p msismall
#SBATCH -c 30
#SBATCH --mem-per-cpu=7000
#SBATCH -t 96:00:00

source /common/software/install/migrated/anaconda/miniconda3_4.8.3-jupyter/etc/profile.d/conda.sh
conda activate r_env
Rscript /home/guanwh/zhan7474/mayo/para_data_sim/exam_sim_data_3_deseq_0811.R