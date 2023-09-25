#!/bin/bash
#SBATCH --job-name=pasim4
#SBATCH -o /home/guanwh/zhan7474/mayo/para_data_sim/exam_sim_data_4_0811.out
#SBATCH -e /home/guanwh/zhan7474/mayo/para_data_sim/exam_sim_data_4_0811.err
#SBATCH -p msismall
#SBATCH -c 20
#SBATCH --mem-per-cpu=10000
#SBATCH -t 96:00:00

module load R/4.3.0-openblas
Rscript /home/guanwh/zhan7474/mayo/para_data_sim/exam_sim_data_4_0811.R