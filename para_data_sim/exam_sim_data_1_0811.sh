#!/bin/bash
#SBATCH --job-name=pasim1
#SBATCH -o /home/guanwh/zhan7474/mayo/para_data_sim/exam_sim_data_1_0811.out
#SBATCH -e /home/guanwh/zhan7474/mayo/para_data_sim/exam_sim_data_1_0811.err
#SBATCH -p msismall
#SBATCH -c 80
#SBATCH --mem-per-cpu=2500
#SBATCH -t 96:00:00

module load R/4.3.0-openblas
Rscript /home/guanwh/zhan7474/mayo/para_data_sim/exam_sim_data_1_0811.R