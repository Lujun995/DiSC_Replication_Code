#!/bin/bash
#SBATCH --job-name=pasim3
#SBATCH -o /home/guanwh/zhan7474/mayo/para_data_sim/exam_sim_data_3_0811.out
#SBATCH -e /home/guanwh/zhan7474/mayo/para_data_sim/exam_sim_data_3_0811.err
#SBATCH -p msismall
#SBATCH -c 30
#SBATCH --mem-per-cpu=7000
#SBATCH -t 96:00:00

module load R/4.3.0-openblas
Rscript /home/guanwh/zhan7474/mayo/para_data_sim/exam_sim_data_3_0811.R