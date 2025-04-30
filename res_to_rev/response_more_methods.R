# response to reviewers: a comparison of power with methods of iDESC and BSDE

##########################
#parameters to be specified:
#general parameters: WORKINGDIRECTORY, 100
#PARA_OUT, SIMU_D, CAT
#parameters in simulations: DIFF, NCELL, NPG, VNCELL
##########################

setwd("/home/guanwh/zhan7474/mayo/res_to_rev/")

library(tidyverse)
library(Matrix)
library(parallel)
library(doParallel)
library(DESeq2)
library(ideas)
library(qvalue)
library(reticulate)
use_python("/home/guanwh/zhan7474/.conda/envs/r_env/bin/python", required = TRUE)
use_condaenv("/home/guanwh/zhan7474/.conda/envs/r_env", required = TRUE)
library(BSDE)
library(iDESC)
source("DiSC.R")


set.seed(seed = 12345181)
n.cores = 100
DATE = "0811"# the date when simulation datasets were created
parameters <- "VNCELL=FALSE;DIFF=1.3;NCELL=375;NPG=12;BOOT=20" # we use the first 20 of the 50 datasets

tryCatch({
  for(paras in parameters){
    cat("Using parameters:", paras)
    eval(parse(text=paras)) # VNCELL=FALSE;DIFF=1.1;NCELL=375;NPG=12;BOOT=50
    boot = BOOT
    diff = exp(abs(log(DIFF))) #promise diff >= 1
    n_cas = NPG
    n_ctr = NPG
    nall = n_cas + n_ctr
    neq_n_cell = VNCELL
    n_cell= NCELL
    
    for(bbbb in 1:boot){
        conf <- str_glue("n{n_cas}_{n_ctr}_d{diff}_c{n_cell}_v{neq_n_cell}_{DATE}_b{bbbb}")
        load(str_glue("/home/guanwh/zhan7474/mayo/para_data_sim/{conf}.RData"))
        count_matrix = dat_list$count_matrix
        meta_cell = dat_list$meta_cell
        meta_cell$phenotype = factor(ifelse(meta_cell$phenotype == 1, "case", "control"))
        names(meta_cell)[which(names(meta_cell) == "cell_rd")] <- "read_depth"
        meta_ind = dat_list$meta_ind
        meta_ind$phenotype = factor(ifelse(meta_ind$phenotype == 1, "case", "control"))
        var2test = "phenotype"
        var2adjust = NULL
        var2test_type = "binary"
        if(bbbb == 1)
            saveRDS(object = dat_list$gene_index, file = str_glue("geneidx_null_{conf}.rds"))
        #source("DiSC.R")
        # using count_matrix, meta_cell, meta_ind, var2test, var2adjust in the memory
        # determine the test type based on the type of meta_ind[[var2test]]
        # always create a "pval_mat" in the memory

        # DiSC
        mth = "DiSC"
        if(!file.exists(str_glue("results_{mth}_{conf}.rds"))){
            # _ may have been used in ctp
            #d: dataset, m: method, v: version, c: cell type
            source(str_glue("{mth}_anal.R")) 
            saveRDS(object = pval_mat, 
                    file = str_glue("results_{mth}_{conf}.rds"))
        }
      
        # IDEAS
        mth = "IDEAS"
        if(!file.exists(str_glue("results_{mth}_{conf}.rds"))){
            registerDoParallel(cores=n.cores)
            options(mc.cores=n.cores)
            RNGkind("L'Ecuyer-CMRG")
            source(str_glue("{mth}_anal.R")) 
            saveRDS(object = pval_mat, 
                    file = str_glue("results_{mth}_{conf}.rds"))
            stopImplicitCluster()
        }
      
        # DESeq2
        mth = "DESeq2"
        if(!file.exists(str_glue("results_{mth}_{conf}.rds"))){
            source(str_glue("{mth}_anal.R")) 
            saveRDS(object = pval_mat, 
                    file = str_glue("results_{mth}_{conf}.rds"))
        }

        # BSDE
        mth = "BSDE"
        if(!file.exists(str_glue("results_{mth}_{conf}.rds"))){
            cl <- makeCluster(n.cores) # BSDE_anal.R requires additionally cl in the memory
            source(str_glue("{mth}_anal.R"))
            stopCluster(cl)
            rm(cl)
            saveRDS(object = pval_mat, 
                    file = str_glue("results_{mth}_{conf}.rds"))
        }

        # iDESC
        mth = "iDESC"
        if(!file.exists(str_glue("results_{mth}_{conf}.rds"))){
            n.cores = n.cores 
            cl <- makeCluster(n.cores) # iDESC_anal.R requires additionally cl in the memory
            source(str_glue("{mth}_anal.R"))
            stopCluster(cl)
            rm(cl)
            saveRDS(object = pval_mat, 
                    file = str_glue("results_{mth}_{conf}.rds"))
        }

        cat(str_glue("finished {conf}\n\n"))
    }
}}, error=function(err){
  # exception handling. See the part above
  # if error occurs, save current work image and restart outside
  print(err)
  print("Need to be restarted...")
  save.image(file = str_glue("unfinished_{conf}.RData"))
  q(save = "no")
  return(0)
})# try-catach part

q(save = "no")
