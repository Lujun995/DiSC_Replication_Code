#!/bin/bash
#SBATCH --job-name=resimType1
#SBATCH -o /home/guanwh/zhan7474/mayo/real_data_sim/Type1_0715.out
#SBATCH -e /home/guanwh/zhan7474/mayo/real_data_sim/Type1_0715.err
#SBATCH -p msismall
#SBATCH -c 80
#SBATCH --mem-per-cpu=2500
#SBATCH -t 60:00:00

module load R/4.3.0-openblas
Rscript /home/guanwh/zhan7474/mayo/real_data_sim/Type1_0715.R