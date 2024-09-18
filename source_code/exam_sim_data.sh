#!/bin/bash
#SBATCH --job-name=pasimCATMETHOD
#SBATCH -o WORKINGDIRECTORYexam_sim_data_CAT_METHOD_DATE.out
#SBATCH -e WORKINGDIRECTORYexam_sim_data_CAT_METHOD_DATE.err
#SBATCH -p msismall
#SBATCH -c NCORE
#SBATCH --mem-per-cpu=MEM
#SBATCH -t 96:00:00

source /common/software/install/migrated/anaconda/miniconda3_4.8.3-jupyter/etc/profile.d/conda.sh
conda activate r_env
Rscript WORKINGDIRECTORYexam_sim_data_CAT_METHOD_DATE.R