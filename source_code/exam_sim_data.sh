#!/bin/bash
#SBATCH --job-name=pasimCAT
#SBATCH -o WORKINGDIRECTORYexam_sim_data_CAT_DATE.out
#SBATCH -e WORKINGDIRECTORYexam_sim_data_CAT_DATE.err
#SBATCH -p msismall
#SBATCH -c NCORE
#SBATCH --mem-per-cpu=MEM
#SBATCH -t 96:00:00

module load R/4.3.0-openblas
Rscript WORKINGDIRECTORYexam_sim_data_CAT_DATE.R