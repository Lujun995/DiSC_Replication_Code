## DiSC
#source("DiSC.R")
# using count_matrix, meta_cell, meta_ind, var2test, var2adjust in the memory
# determine the test type based on the type of meta_ind[[var2test]]
# always create a "pval_mat" in the memory
pval_mat <- matrix(data = 1, nrow = 2, ncol = nrow(count_matrix), 
                   dimnames = list(c("adjusted", "raw"),
                                   rownames(count_matrix)))
obj <- DiSC(data.mat = count_matrix, 
            cell.ind = meta_cell,
            metadata = meta_ind,
            outcome = var2test, 
            covariates = var2adjust,
            cell.id = "cell_id", individual.id = "individual", 
            features = c('prev', 'nzm', 'nzsd'), 
            verbose = FALSE, perm.no = 999, sequencing.data = TRUE)
pval_mat["adjusted", ] <- obj$p.adj.fdr
pval_mat["raw", ] <- obj$p.raw