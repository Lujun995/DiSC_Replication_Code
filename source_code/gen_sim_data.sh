#!/bin/bash
#SBATCH --job-name=gen_d
#SBATCH -o WORKINGDIRECTORYgen_sim_data_DATE.out
#SBATCH -e WORKINGDIRECTORYgen_sim_data_DATE.err
#SBATCH -p msismall
#SBATCH -c NCORE
#SBATCH --mem-per-cpu=MEM
#SBATCH -t 96:00:00

module load R/4.3.0-openblas
Rscript WORKINGDIRECTORYgen_sim_data_DATE.R
