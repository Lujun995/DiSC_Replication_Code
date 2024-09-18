#!/bin/bash
#SBATCH --job-name=resimCLASS
#SBATCH -o WORKINGDIRECTORYCLASS_METHOD_DATE.out
#SBATCH -e WORKINGDIRECTORYCLASS_METHOD_DATE.err
#SBATCH -p msismall
#SBATCH -c NCORES
#SBATCH --mem-per-cpu=MEM
#SBATCH -t 60:00:00

source /common/software/install/migrated/anaconda/miniconda3_4.8.3-jupyter/etc/profile.d/conda.sh
conda activate r_env
Rscript WORKINGDIRECTORYCLASS_METHOD_DATE.R