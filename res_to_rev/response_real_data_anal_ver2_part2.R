# AD_anal.R
library(tidyverse)
library(Matrix)
library(parallel)
library(doParallel)
library(DESeq2)
library(ideas)
library(qvalue)
# library(data.table)
# rm(list=ls())

setwd("/home/guanwh/zhan7474/mayo/res_to_rev/")
source("DiSC.R")

set.seed(seed = 123456)
n.cores = 80
DDDD = "0119"
data_dir = "/home/guanwh/zhan7474/mayo/res_to_rev/AD2"
dtset = "AD2"
files = list.files(data_dir)
ctypes = str_remove(files[str_detect(files, ".RData$")], ".RData")

for(ctp in ctypes){
  load(str_glue("{data_dir}/{ctp}.RData")) 
  # this dataset should have:
  # count_matrix: genes in rows, cells in cols, rowname: gene_name, colname: cell_id
  # meta_cell: with variables called "individual" (string), "cell_id" (string) and "read_depth" (numeric)
  # meta_ind: with variables called "individual", the outcome and all the covariates (numeric or factor)
  # var2test, var2adjust, var2test_type: strings
  if(inherits(count_matrix, "dgCMatrix"))
    count_matrix <- as.matrix(count_matrix)
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
      cat("Method start:", mth, "\n")
      print(Sys.time())
      if(!file.exists(str_glue("AD2/d@{dtset}-m@{mth}-v@{DDDD}-c@{ctp}.rds"))){
        # _ may have been used in ctp
        #d: dataset, m: method, v: version, c: cell type
        source(str_glue("{mth}_anal.R")) 
        saveRDS(object = pval_mat, 
                file = str_glue("AD2/d@{dtset}-m@{mth}-v@{DDDD}-c@{ctp}.rds"))
      }
      cat("Method complete:", mth, "\n")
      print(Sys.time())
      
      # IDEAS
      mth = "IDEAS"
      cat("Method start:", mth, "\n")
      print(Sys.time())
      if(!file.exists(str_glue("AD2/d@{dtset}-m@{mth}-v@{DDDD}-c@{ctp}.rds"))){
        registerDoParallel(cores=n.cores)
        options(mc.cores=n.cores)
        RNGkind("L'Ecuyer-CMRG")
        # source(str_glue("{mth}_anal.R"))
        source(str_glue("IDEAS_anal_faster.R"))	
        saveRDS(object = pval_mat, 
                file = str_glue("AD2/d@{dtset}-m@{mth}-v@{DDDD}-c@{ctp}.rds"))
        stopImplicitCluster()
      }
      cat("Method complete:", mth, "\n")
      print(Sys.time())
      
      # DESeq2
      mth = "DESeq2"
      cat("Method start:", mth, "\n")
      print(Sys.time())
      if(!file.exists(str_glue("AD2/d@{dtset}-m@{mth}-v@{DDDD}-c@{ctp}.rds"))){
        source(str_glue("{mth}_anal.R")) 
        saveRDS(object = pval_mat, 
                file = str_glue("AD2/d@{dtset}-m@{mth}-v@{DDDD}-c@{ctp}.rds"))
      }
      cat("Method complete:", mth, "\n")
      print(Sys.time())

      # BSDE
      mth = "BSDE"
      cat("Method start:", mth, "\n")
      print(Sys.time())
      if(!file.exists(str_glue("AD2/d@{dtset}-m@{mth}-v@{DDDD}-c@{ctp}.rds"))){
          cl <- makeCluster(n.cores) # BSDE_anal.R requires additionally cl in the memory
          source(str_glue("{mth}_anal.R"))
          stopCluster(cl)
          rm(cl)
          saveRDS(object = pval_mat, 
                  file = str_glue("AD2/d@{dtset}-m@{mth}-v@{DDDD}-c@{ctp}.rds"))
      }
      cat("Method complete:", mth, "\n")
      print(Sys.time())

      # iDESC
      mth = "iDESC"
      cat("Method start:", mth, "\n")
      print(Sys.time())
      if(!file.exists(str_glue("AD2/d@{dtset}-m@{mth}-v@{DDDD}-c@{ctp}.rds"))){
          cl <- makeCluster(n.cores) # iDESC_anal.R requires additionally cl in the memory
          source(str_glue("{mth}_anal.R"))
          stopCluster(cl)
          rm(cl)
          saveRDS(object = pval_mat, 
                  file = str_glue("AD2/d@{dtset}-m@{mth}-v@{DDDD}-c@{ctp}.rds"))
      }
      cat("Method complete:", mth, "\n")
      print(Sys.time())

      
      # More methods... (okay to rerun this script after adding new methods)
      
    }, error = function(err) {
      if(exists("cl")){
        stopCluster(cl)
        rm(cl)
      }
      print(err)
      save.image(str_glue("unfinished_{dtset}_{ctp}_{DDDD}.RData"))
    }) #try-catch part, when error occurs, save the current working space and go to the next ctp
  
  rm(list = c("count_matrix", "meta_cell", "meta_ind", 
              "var2test", "var2adjust", "var2test_type"))
}# go to the next ctp
