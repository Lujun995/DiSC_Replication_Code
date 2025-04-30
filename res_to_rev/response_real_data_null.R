# The authors should consider conducting additional permutations in PBMC analysis to provide 
# a more accurate comparison.

library(tidyverse)
library(Matrix)
library(parallel)
library(doParallel)
library(DESeq2)
library(ideas)
library(qvalue)

setwd("/home/guanwh/zhan7474/mayo/res_to_rev/")
source("DiSC.R")

set.seed(seed = 123456)
n.cores = 100
DDDD = "1231"
# preprocessed data coming from real_data_anal/COVID_PBMC
# this dataset should have:
# count_matrix: genes in rows, cells in cols, rowname: gene_name, colname: cell_id
# meta_cell: with variables called "individual" (string), "cell_id" (string) and "read_depth" (numeric)
# meta_ind: with variables called "individual", the outcome and all the covariates (numeric or factor)
# var2test, var2adjust, var2test_type: strings
load("/home/guanwh/zhan7474/mayo/real_data_anal/COVID_PBMC/CD4.RData") 
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
meta_ind_original <- meta_ind

boot = 1000
nindi = nrow(meta_ind_original)
perm.index = sapply(1:boot, simplify = TRUE,
                    FUN = function(xxxx) 
                      return(sample(1:nindi)) )

cl <- makeCluster(round(n.cores/1.5)) #just to be OOM safe...
clusterEvalQ(cl, {
  library(tidyverse)
  library(Matrix)
  library(DESeq2)

  setwd("/home/guanwh/zhan7474/mayo/res_to_rev/")
  source("DiSC.R")

  load("/home/guanwh/zhan7474/mayo/real_data_anal/COVID_PBMC/CD4.RData") 
  meta_ind_original <- meta_ind
  NULL # return value
})
clusterExport(cl, c("DDDD", "perm.index"))
parSapplyLB(cl = cl, X = 1:boot, FUN = function(bbbb){
# parSapply(cl = cl, X = 1:100, FUN = function(bbbb){
  meta_ind <- meta_ind_original
  meta_ind[[var2test]] <- meta_ind_original[[var2test]][perm.index[,bbbb]]
  # DiSC
  mth = "DiSC"
  if(!file.exists(str_glue("null_res/d@COVIDPBMC-m@{mth}-v@{DDDD}-b@{bbbb}.rds"))){
        #d: dataset, m: method, v: version, b: iteration
        source(str_glue("{mth}_anal.R"), local = TRUE) 
        saveRDS(object = pval_mat, 
                file = str_glue("null_res/d@COVIDPBMC-m@{mth}-v@{DDDD}-b@{bbbb}.rds"))
  }
  # DESeq2
  mth = "DESeq2"
  if(!file.exists(str_glue("null_res/d@COVIDPBMC-m@{mth}-v@{DDDD}-b@{bbbb}.rds"))){
        source(str_glue("{mth}_anal.R"), local = TRUE) 
        saveRDS(object = pval_mat, 
                file = str_glue("null_res/d@COVIDPBMC-m@{mth}-v@{DDDD}-b@{bbbb}.rds"))
  }
  return(NULL)
})
stopCluster(cl)
rm(list = "cl")

for(bbbb in 1:boot){
  cat("Current iteration:", bbbb, "\n")
  meta_ind <- meta_ind_original
  meta_ind[[var2test]] <- meta_ind_original[[var2test]][perm.index[,bbbb]]

  tryCatch(
    {
      # DiSC
      mth = "DiSC"
      if(!file.exists(str_glue("null_res/d@COVIDPBMC-m@{mth}-v@{DDDD}-b@{bbbb}.rds"))){
        #d: dataset, m: method, v: version, b: iteration
        source(str_glue("{mth}_anal.R")) 
        saveRDS(object = pval_mat, 
                file = str_glue("null_res/d@COVIDPBMC-m@{mth}-v@{DDDD}-b@{bbbb}.rds"))
      }
      
      # IDEAS
      mth = "IDEAS"
      if(!file.exists(str_glue("null_res/d@COVIDPBMC-m@{mth}-v@{DDDD}-b@{bbbb}.rds"))){
        registerDoParallel(cores=n.cores)
        options(mc.cores=n.cores)
        RNGkind("L'Ecuyer-CMRG")
        source(str_glue("{mth}_anal.R")) 
        saveRDS(object = pval_mat, 
                file = str_glue("null_res/d@COVIDPBMC-m@{mth}-v@{DDDD}-b@{bbbb}.rds"))
        stopImplicitCluster()
      }
      
      # DESeq2
      mth = "DESeq2"
      if(!file.exists(str_glue("null_res/d@COVIDPBMC-m@{mth}-v@{DDDD}-b@{bbbb}.rds"))){
        source(str_glue("{mth}_anal.R")) 
        saveRDS(object = pval_mat, 
                file = str_glue("null_res/d@COVIDPBMC-m@{mth}-v@{DDDD}-b@{bbbb}.rds"))
      }
      
      # More methods... (okay to rerun this script after adding new methods)
      
    }, error = function(err) {
      print(err)
      save.image(str_glue("unfinished_COVIDPBMC_b{bbbb}_{DDDD}.RData"))
    }) #try-catch part, when error occurs, save the current working space and go to the next ctp

}# go to the next ctp


