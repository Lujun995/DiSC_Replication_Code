##########################
#parameters to be specified:
#general parameters: WORKDIR, NCORES, DATE, DATADIR
##########################

library(tidyverse)
library(Matrix)
library(parallel)
library(doParallel)
library(DESeq2)
library(ideas)
library(qvalue)
# library(data.table)
# rm(list=ls())

setwd("WORKDIR")
source("DiSC.R")

set.seed(seed = 123456)
n.cores = NCORES
DDDD = "DATE"
data_dir = "DATADIR"
# CTYPES = str_c(str_glue("'{str_remove(list.files(DATADIR), ".rds")}'"), collapse = ",")
files = list.files(data_dir)
ctypes = str_remove(files[str_detect(files, ".RData$")], ".RData")
pval_array_list <- vector(mode = "list", length = length(ctypes))
names(pval_array_list) <- ctypes

for(ctp in ctypes){
  load(str_glue("{data_dir}/{ctp}.RData")) 
  # this dataset should have:
  # count_matrix: genes in rows, cells in cols, rowname: gene_name, colname: cell_id
  # meta_cell: with variables called "individual" (string), "cell_id" (string) and "read_depth" (numeric)
  # meta_ind: with variables called "individual", the outcome and all the covariates (numeric or factor)
  # var2test, var2adjust, var2test_type: strings
  cat("Using cell type:", ctp, "\n")
  cat("Dimension of the count_matrix:", dim(count_matrix), "\n")
  cat("Using outcome:", var2test, "which is", class(meta_ind[[var2test]]), "\n")
  cat("Using covariates:", var2adjust, "\n")
  for(var in var2adjust) 
    cat(str_glue("class of {var}: {class(meta_ind[[var]])}"), "\n")
  cat("Using test type:", var2test_type, "\n")
  stopifnot(sum(colnames(count_matrix) != meta_cell$cell_id) == 0)
  stopifnot(setequal(unique(meta_cell$individual), meta_ind$individual))
  stopifnot(class(meta_cell$individual) == "character" &
              class(meta_cell$cell_id) == "character" &
              class(meta_ind$individual) == "character")
  tryCatch(
    {
      # DiSC
      mth = "DiSC"
      if(!file.exists(str_glue("res/d@COVIDPBMC-m@{mth}-v@{DDDD}-c@{ctp}.rds"))){
        # _ may have been used in ctp
        #d: dataset, m: method, v: version, c: cell type
        source(str_glue("{mth}_anal.R")) 
        saveRDS(object = pval_mat, 
                file = str_glue("res/d@COVIDPBMC-m@{mth}-v@{DDDD}-c@{ctp}.rds"))
      }
      
      # IDEAS
      mth = "IDEAS"
      if(!file.exists(str_glue("res/d@COVIDPBMC-m@{mth}-v@{DDDD}-c@{ctp}.rds"))){
        registerDoParallel(cores=n.cores)
        options(mc.cores=n.cores)
        RNGkind("L'Ecuyer-CMRG")
        source(str_glue("{mth}_anal.R")) 
        saveRDS(object = pval_mat, 
                file = str_glue("res/d@COVIDPBMC-m@{mth}-v@{DDDD}-c@{ctp}.rds"))
        stopImplicitCluster()
      }
      
      # DESeq2
      mth = "DESeq2"
      if(!file.exists(str_glue("res/d@COVIDPBMC-m@{mth}-v@{DDDD}-c@{ctp}.rds"))){
        source(str_glue("{mth}_anal.R")) 
        saveRDS(object = pval_mat, 
                file = str_glue("res/d@COVIDPBMC-m@{mth}-v@{DDDD}-c@{ctp}.rds"))
      }
      
      # More methods... (okay to rerun this script after adding new methods)
      
    }, error = function(err) {
      print(err)
      save.image(str_glue("unfinished_COVID_PBMC_{ctp}_{DDDD}.RData"))
    }) #try-catch part, when error occurs, save the current working space and go to the next ctp
  
  rm(list = c("count_matrix", "meta_cell", "meta_ind", 
              "var2test", "var2adjust", "var2test_type"))
}# go to the next ctp
