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

method = "deseq" # need to specify method when creating this file.
# Assuming a variable called "{method}_pval" will be saved
require(DESeq2)

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
  deseq_pval_adj <- DESeq2::results(dds, name = "diagnosis_ASD_vs_Control")$padj
  deseq_pval_adj[is.na(deseq_pval_adj)] <- 1 
  # when allowing independentFiltering and cooksCutoff, p-values for some
  # genes may be returned as NA. Set them to an insignificant p-value.
  deseq_pval_raw <- DESeq2::results(dds, name = "diagnosis_ASD_vs_Control")$pvalue
  deseq_pval_raw[is.na(deseq_pval_raw)] <- 1
  deseq_pval <- matrix(c(deseq_pval_raw, deseq_pval_adj),
                       ncol = 2, byrow = FALSE)
  colnames(deseq_pval) <- c("raw", "adjusted")
  
  return(deseq_pval)
})

save(list = str_glue("{method}_pval"),
     file = str_glue("{cclass}_{method}_{ddate}.RData"))

q(save = "no")