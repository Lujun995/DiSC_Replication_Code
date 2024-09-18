## DESeq2, pseudo-bulk sample method
library(DESeq2)
library(stringr)
# using count_matrix, meta_cell, meta_ind, var2test, var2adjust in the memory
# determine the test type based on the type of meta_ind[[var2test]]
# always create a "pval_mat" in the memory
pval_mat <- matrix(data = 1, nrow = 2, ncol = nrow(count_matrix), 
                   dimnames = list(c("adjusted", "raw"),
                                   rownames(count_matrix)))
ngene = nrow(count_matrix)
nindi = nrow(meta_ind)
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

design = as.formula(paste("~", paste(c(var2adjust, var2test), collapse = "+")))
dds <- DESeqDataSetFromMatrix(countData = count_matrix_bulk,
                              colData = meta_ind,
                              design = design)
dds <- DESeq(dds, quiet = TRUE)
res_name = resultsNames(dds)[str_detect(resultsNames(dds), var2test)] 
# assume var2test is a binary factor or a continuous variable
pval_mat["adjusted", ] <- results(dds, name = res_name)$padj
pval_mat["raw", ] <- results(dds, name = res_name)$pvalue

