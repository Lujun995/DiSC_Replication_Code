##########################
#parameters to be specified:
#WORKINGDIRECTORY, BOOT, NCORES, DATE
##########################

#set up the variables for the methods and for the parallel computing
setwd("WORKINGDIRECTORY")
boot = BOOT
n.cores = NCORES
ddate = "DATE"
cclass = "CLASS"
set.seed(seed = 42) # need to further set the random seed for parallel computing

library(stringr)
library(parallel)
load("setup0814.RData")

source("DiSC.R")
method = "disc" # need to specify method when creating this file.
require(matrixStats)
# Assuming a variable called "{method}_pval" will be saved

# library(Matrix)
# library(doParallel)
# library(ideas)
# library(qvalue)
# require(DESeq2)

# DiSC
cl <- makeCluster(n.cores)
clusterExport(cl, c("meta_cell", "meta_ind", 
                    "count_matrix", "perm.index",
                    "na.pad", "perm.fdr.adj", "perm.fwer.adj", 
                    "DiSC", "basic.func"))
clusterEvalQ(cl, {
  library(matrixStats)
})
disc_pval <- parSapply(cl, 1:boot, simplify = FALSE, FUN = function(xxxx){
  meta_ind_p <- meta_ind
  meta_ind_p$diagnosis <- meta_ind$diagnosis[perm.index[,xxxx]]
  obj <- DiSC(data.mat = count_matrix, 
              cell.ind = meta_cell, 
              metadata = meta_ind_p,
              outcome = "diagnosis", covariates = "RIN", 
              cell.id = "cell_id",
              individual.id = "individual", 
              verbose = FALSE, perm.no = 999, sequencing.data = TRUE)
  disc_pval_adj <- obj$p.adj.fdr
  disc_pval_raw <- obj$p.raw
  disc_pval <- matrix(c(disc_pval_raw, disc_pval_adj),
                      ncol = 2, byrow = FALSE)
  colnames(disc_pval) <- c("raw", "adjusted")
  return(disc_pval)
})
stopCluster(cl)

save(list = str_glue("{method}_pval"),
     file = str_glue("{cclass}_{method}_{ddate}.RData"))

q(save = "no")