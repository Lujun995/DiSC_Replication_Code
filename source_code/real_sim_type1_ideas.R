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

method = "ideas" # need to specify method when creating this file.
# Assuming a variable called "{method}_pval" will be saved

library(Matrix)
library(doParallel)
library(ideas)
library(qvalue)

#ideas
registerDoParallel(cores=n.cores)
options(mc.cores=n.cores)
ideas_pval <- sapply(1:boot, simplify = FALSE, FUN = function(xxxx){
  # if(xxxx %% 50 == 0) cat(xxxx)
  meta_ind_p <- meta_ind
  meta_ind_p$diagnosis <- meta_ind$diagnosis[perm.index[,xxxx]]
  # rate-limiting step
  dist1 = ideas_dist(count_input = count_matrix, 
                     meta_cell = meta_cell, meta_ind = meta_ind_p, 
                     var_per_cell = "read_depth", var2test = "diagnosis", 
                     var2test_type = "binary",
                     d_metric = "Was", fit_method = "nb")
  ideas_pval_raw <- permanova(dist_array = dist1, meta_ind = meta_ind_p, 
                              var2test = "diagnosis", var2adjust = "RIN", 
                              var2test_type = "binary", n_perm=999)
  recal_index <- which(ideas_pval_raw < 0.01)# NA values will be excluded
  # p.true = 0.001, prob(p.est >= 0.008 using n.perm = 999) = 8e-5 
  if(length(recal_index) > 0){
    if (length(recal_index) == 1){
      slice_dist = dist1[recal_index, , ]  
      recal_dist = array(dim = c(1, dim(slice_dist)))
      recal_dist[1, , ] = slice_dist
    }else{
      recal_dist = dist1[recal_index, , ]
    }
    ideas_pval_raw_r <- permanova(dist_array = recal_dist, meta_ind = meta_ind_p, 
                                  var2test = "diagnosis", var2adjust = "RIN", 
                                  var2test_type = "binary", n_perm=99999)
    ideas_pval_raw[recal_index] <- ideas_pval_raw_r
  }
  
  ideas_pval_adj <- qvalue::qvalue(ideas_pval_raw)$qvalues
  ideas_pval <- matrix(c(ideas_pval_raw, ideas_pval_adj),
                       ncol = 2, byrow = FALSE)
  colnames(ideas_pval) <- c("raw", "adjusted")
  return(ideas_pval)
})
stopImplicitCluster()

save(list = str_glue("{method}_pval"),
     file = str_glue("{cclass}_{method}_{ddate}.RData"))

q(save = "no")
