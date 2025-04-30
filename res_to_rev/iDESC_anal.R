# iDESC
# using count_matrix, meta_cell, meta_ind, var2test, var2adjust, var2test_type in the memory
# always create a "pval_mat" in the memory
source("zp_prediction.r")
source("iDESC.r")
library(iDESC)

stopifnot(is.null(var2adjust)) # iDESC can not handle the covariates
# stopifnot(var2test_type == "binary") # iDESC can only handle binary outcome (?) 
# iDESC should test the outcome variable according to its data type
pval_mat <- matrix(data = 1, nrow = 2, ncol = nrow(count_matrix), 
                   dimnames = list(c("adjusted", "raw"),
                                   rownames(count_matrix)))
meta_cell_iDESC <- meta_cell[, c("cell_id", "individual")]
meta_cell_iDESC <- dplyr::left_join(x = meta_cell_iDESC, y = meta_ind, by = "individual")
rownames(meta_cell_iDESC) <- meta_cell_iDESC$cell_id

if(exists("cl")){
    if(inherits(cl, "cluster")){
        # crr_wd <- getwd()
        # l_cl = length(cl)
        # step = nrow(count_matrix) %/% l_cl
        # remainder = nrow(count_matrix) %% l_cl
        # cal_vec = rep(step, l_cl)
        # if(remainder > 0)
        #     cal_vec[1:remainder] <- cal_vec[1:remainder] + 1
        # rm(list = c("l_cl", "step", "remainder"))
        clusterEvalQ(cl, {
            library(iDESC)
            library(stringr)
            # library(stringr)
            # setwd(crr_wd)
        })
        # clusterExport(cl, c("count_matrix", "meta_cell_iDESC", "var2test", 
        #                     "cal_vec", "zp_prediction", "iDESC"))
        
        result <- iDESC(count_matrix, meta_cell_iDESC, 
                        subject_var="individual", group_var=var2test,
                        norm_opt="SeqDepth", span = 0.7, # according to "function usage"
                        sub_cell_filtering = 0, gene_sub_filtering = 0, 
                        gene_cell_filtering = 0, ncell_filtering = 0, 
                        cl = cl)
        pval_mat["raw", ] <- result[[stringr::str_subset(names(result), "Pval_Beta")]]
    }
} else {

    result <- iDESC(count_matrix, meta_cell_iDESC, 
                    subject_var="individual", group_var=var2test,
                    norm_opt="SeqDepth", span = 0.7, # according to "function usage"
                    sub_cell_filtering = 0, gene_sub_filtering = 0, 
                    gene_cell_filtering = 0, ncell_filtering = 0)
    pval_mat["raw", ] <- result[[stringr::str_subset(names(result), "Pval_Beta")]]
}
pval_mat["adjusted", ] <- p.adjust(pval_mat["raw", ], method = "BH")