##########################
#parameters to be specified:
#/home/guanwh/zhan7474/mayo/real_data_sim/, 1000, 80, 0715
##########################

#set up the variables for the methods and for the parallel computing
setwd("/home/guanwh/zhan7474/mayo/real_data_sim/")
library(Matrix)
library(tidyverse)
library(parallel)
library(doParallel)
load("setup0702.RData")
set.seed(seed = 123456) # need to further set the random seed for parallel computing
boot = 1000
n.cores = 80
meta_cell = as.data.frame(meta_cell) 
meta_ind = as.data.frame(meta_ind)
count_matrix = as.matrix(count_matrix)
perm.index = sapply(1:boot, simplify = TRUE,
                    FUN = function(xxxx) 
                      return(sample(1:nindi)) )
rarefy_mat <- count_matrix %>% t() %>% 
  GUniFrac::Rarefy(.) %>%
  magrittr::extract2(., "otu.tab.rff") %>%
  t()
count_matrix_bulk = matrix(NA, nrow = ngene, ncol = nindi)
rownames(count_matrix_bulk) = rownames(count_matrix)
colnames(count_matrix_bulk) = meta_ind$individual
if(sum(meta_cell$cell_id != colnames(count_matrix))==0)
  for (i_ind in 1:length(meta_ind$individual)) {
    cur_ind   = meta_ind$individual[i_ind]
    cur_ind_m = count_matrix[, meta_cell$individual == cur_ind] 
    count_matrix_bulk[, i_ind] = rowSums(cur_ind_m, na.rm = TRUE)
  }else
    error("meta_cell$cell_id and colnames(count_matrix) dont match")

# DiSC
gc()
cl <- makeCluster(round(n.cores/2.5))
clusterExport(cl, c("meta_cell", "meta_ind", 
                    "rarefy_mat", "perm.index",
                    "na.pad", "perm.fdr.adj", "perm.fwer.adj", 
                    "DiSC", "basic.func"))
clusterEvalQ(cl, {
  library(matrixStats)
})
disc_pval <- parSapply(cl, 1:boot, simplify = FALSE, FUN = function(xxxx){
  meta_ind_p <- meta_ind
  meta_ind_p$diagnosis <- meta_ind$diagnosis[perm.index[,xxxx]]
  obj <- DiSC(ct.mat = rarefy_mat, 
              cell.ind = meta_cell, 
              metadata = meta_ind_p,
              outcome = "diagnosis", covariates = "RIN", 
              cell.id = "cell_id",
              individual.id = "individual")
  disc_pval_adj <- obj$p.adj.fdr
  disc_pval_raw <- obj$p.raw
  disc_pval <- matrix(c(disc_pval_raw, disc_pval_adj),
                      ncol = 2, byrow = FALSE)
  colnames(disc_pval) <- c("raw", "adjusted")
  return(disc_pval)
})
stopCluster(cl)

#ideas
library(ideas)
# library(qvalue)
registerDoParallel(cores=n.cores)
options(mc.cores=n.cores)
ideas_pval <- sapply(1:boot, simplify = FALSE, FUN = function(xxxx){
  if(xxxx %% 50 == 0) cat(xxxx)
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

#pseudo_bulk sample method, DEseq2
cl <- makeCluster(n.cores)
clusterExport(cl, c("meta_ind", "perm.index", "count_matrix_bulk"))
clusterEvalQ(cl, {
  library(DESeq2)
})
deseq_pval <- parSapply(cl, 1:boot, simplify = FALSE, FUN = function(xxxx){
  meta_ind_p <- meta_ind
  meta_ind_p$diagnosis <- meta_ind$diagnosis[perm.index[,xxxx]]
  
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = count_matrix_bulk,
                                        colData = meta_ind_p,
                                        design = ~ RIN + diagnosis)
  dds <- DESeq2::DESeq(dds, quiet = TRUE)
  deseq_pval_adj <- DESeq2::results(dds, name = "diagnosis_ASD_vs_Control", 
                                    independentFiltering=FALSE,
                                    cooksCutoff=FALSE)$padj
  deseq_pval_raw <- DESeq2::results(dds, name = "diagnosis_ASD_vs_Control",
                                    independentFiltering=FALSE,
                                    cooksCutoff=FALSE)$pvalue
  deseq_pval <- matrix(c(deseq_pval_raw, deseq_pval_adj),
                       ncol = 2, byrow = FALSE)
  colnames(deseq_pval) <- c("raw", "adjusted")
  
  return(deseq_pval)
})

save.image("type1_0715.RData")