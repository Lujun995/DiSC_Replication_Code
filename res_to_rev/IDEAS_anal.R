## IDEAS
library(ideas)
library(qvalue)
# using count_matrix, meta_cell, meta_ind, var2test, var2adjust, var2test_type in the memory
# always create a "pval_mat" in the memory
pval_mat <- matrix(data = 1, nrow = 2, ncol = nrow(count_matrix), 
                   dimnames = list(c("adjusted", "raw"),
                                   rownames(count_matrix)))
dist1 = ideas_dist(count_input = count_matrix, 
                   meta_cell = meta_cell, meta_ind = meta_ind, 
                   var_per_cell = "read_depth", var2test = var2test, 
                   var2test_type = var2test_type,
                   d_metric = "Was", fit_method = "nb")
ideas_pval_raw <- permanova(dist_array = dist1, meta_ind = meta_ind, 
                            var2test = var2test, var2adjust = var2adjust, 
                            var2test_type = var2test_type, n_perm=999)
recal_index <- which(ideas_pval_raw < 0.01)# NA values will be excluded
if(length(recal_index) > 0){
  if (length(recal_index) == 1){
    slice_dist = dist1[recal_index, , ]  
    recal_dist = array(dim = c(1, dim(slice_dist)))
    recal_dist[1, , ] = slice_dist
  }else{
    recal_dist = dist1[recal_index, , ]
  }
  ideas_pval_raw_r <- permanova(dist_array = recal_dist, meta_ind = meta_ind, 
                                var2test = var2test, var2adjust = var2adjust, 
                                var2test_type = var2test_type, n_perm=99999)
  ideas_pval_raw[recal_index] <- ideas_pval_raw_r
}
pval_temp <- ideas_pval_raw
pi0 = 2*sum(pval_temp > 0.5, na.rm = TRUE)/sum(!is.na(pval_temp))
# really need this pi0, otherwise in some cases, error will occur
pi0 = ifelse(pi0>1, 1, ifelse(pi0<0.5, 0.5, pi0)) # 0.5<=pi0<=1
ideas_pval_adj <- qvalue(pval_temp, pi0=pi0)$qvalues

pval_mat["adjusted", ] <- ideas_pval_adj
pval_mat["raw", ] <- ideas_pval_raw