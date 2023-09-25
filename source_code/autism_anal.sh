#!/bin/bash
#SBATCH --job-name=autism
#SBATCH -o WORKINGDIRECTORYautism_anal_DATE.out
#SBATCH -e WORKINGDIRECTORYautism_anal_DATE.err
#SBATCH -p msismall
#SBATCH -c NCORE
#SBATCH --mem-per-cpu=MEM
#SBATCH -t 24:00:00

module load R/4.3.0-openblas
# #################
# # the following steps are for dca
# #################
# conda deactivate
# cd /home/guanwh/zhan7474/mayo/ideas_pipeline-main/Autism/

# Rscript /home/guanwh/zhan7474/mayo/ideas_pipeline-main/Autism/step1a_dca_prepare_data.R

# # assume using singularity
# # to install tensorflow environment and dca
# module load singularity
# singularity pull docker://tensorflow/tensorflow:2.4.3
# singularity run tensorflow_2.4.3.sif
# singularity exec -i tensorflow_2.4.3.sif  dpkg -l
# singularity exec tensorflow_2.4.3.sif pip install dca
# # need to download micro (a text editor) and edit a python file.

# # assume dca has been installed
# module load singularity
# mkdir /home/guanwh/zhan7474/mayo/ideas_pipeline-main/Autism/dca_PFC_all
# singularity exec --bind /home/guanwh/zhan7474/mayo/ideas_pipeline-main/Autism/data:/container/input \
                 # --bind /home/guanwh/zhan7474/mayo/ideas_pipeline-main/Autism/dca_PFC_all:/container/output \
                 # /home/guanwh/zhan7474/tensorflow_2.4.3.sif \
                 # /home/guanwh/zhan7474/.local/bin/dca /container/input/PFC_all.csv /container/output/dca_PFC_all --threads NCORE --type zinb-conddisp

# Rscript /home/guanwh/zhan7474/mayo/ideas_pipeline-main/Autism/step1a_dca_recover_mean_norm.R
# Rscript /home/guanwh/zhan7474/mayo/ideas_pipeline-main/Autism/step1a_split_dca_outputs.R
# #################
# # dca completed
# #################

Rscript WORKINGDIRECTORYautism_anal_DATE.R

