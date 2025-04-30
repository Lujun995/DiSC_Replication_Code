# using count_matrix, meta_cell, meta_ind, var2test, var2adjust in the memory
# always create a "pval_mat" in the memory

# library(reticulate)
# use_python("/home/guanwh/zhan7474/.conda/envs/r_env/bin/python", required = TRUE) # not working. BSDE forced to use a conda env called r-reticulate
# use_condaenv("/home/guanwh/zhan7474/.conda/envs/r_env", required = TRUE)
library(BSDE)
library(qvalue)
source("cal_w2_pval.r")

stopifnot(is.null(var2adjust)) # BSDE can not handle the covariates
stopifnot(var2test_type == "binary") # BSDE can only handle binary outcome (?) 

pval_mat <- matrix(data = 1, nrow = 2, ncol = nrow(count_matrix), 
                   dimnames = list(c("adjusted", "raw"),
                                   rownames(count_matrix)))

meta_cell_BSDE <- meta_cell[, c("cell_id", "individual")]
# meta_cell_BSDE <- merge(x = meta_cell_BSDE, y = meta_ind, by = "individual", 
                        # all.x = TRUE, sort = FALSE) # even if sort=F, merge still sort the results
meta_cell_BSDE <- dplyr::left_join(x = meta_cell_BSDE, y = meta_ind, by = "individual")
if(inherits(meta_cell_BSDE[[var2test]], "factor"))
  meta_cell_BSDE[[var2test]] <- as.numeric(meta_cell_BSDE[[var2test]]) - 1 # convert to 0/1
stopifnot(inherits(meta_cell_BSDE[[var2test]], "numeric"))
stopifnot(unique(meta_cell_BSDE[[var2test]]) %in% c(0,1))
stopifnot(colnames(count_matrix) == meta_cell_BSDE$cell_id)

sim_matrix_log = log2(1 + count_matrix)

if(exists("cl")){
  if(inherits(cl, "cluster")){
    clusterEvalQ(cl, {
      library(reticulate)
      use_python("/home/guanwh/zhan7474/.conda/envs/r_env/bin/python", required = TRUE)
      use_condaenv("/home/guanwh/zhan7474/.conda/envs/r_env", required = TRUE)
      library(BSDE)
      NULL # return value
    })
    clusterExport(cl, c("sim_matrix_log", "meta_cell_BSDE", "var2test", 
                        "cal_w2_pval"))
    
    pval_mat["raw", ] <- parSapplyLB(cl = cl, X = 1:nrow(count_matrix), simplify = TRUE, FUN = function(i_g){
      cur_sim = sim_matrix_log[i_g, ]
      cur_ind = meta_cell_BSDE[["individual"]]
      cur_pheno = meta_cell_BSDE[[var2test]]
      pval <- tryCatch({
        cal_w2_pval(
          count_per_gene = cur_sim,
          meta_individual = cur_ind,
          meta_phenotype = cur_pheno,
          perm_num = 999,
          unif_round_unit = 0.5
        )[[1]] # pval
      }, error = function(e) {
        NA
      })
      return(pval)
    })
    
    recal_index <- which(pval_mat["raw", ] < 0.01) # NA values will be excluded
    if(length(recal_index) >0){
      cat("raw p-val calculated. increase the number of permutations in", length(recal_index),"genes\n")
      pval_mat["raw", recal_index] <- parSapplyLB(cl = cl, X = recal_index, simplify = TRUE, FUN = function(i_g){
        cur_sim = sim_matrix_log[i_g, ]
        cur_ind = meta_cell_BSDE[["individual"]]
        cur_pheno = meta_cell_BSDE[[var2test]]
        pval <- tryCatch({
          cal_w2_pval(
            count_per_gene = cur_sim,
            meta_individual = cur_ind,
            meta_phenotype = cur_pheno,
            perm_num = 9999,
            unif_round_unit = 0.5
          )[[1]] # pval
        }, error = function(e) {
          NA
        })
        return(pval)
      })
    }
    
  }
} else {
  for (i_g in 1:nrow(count_matrix)) {
    if(i_g %% 1000 == 0)  cat(i_g, "\n")
    cur_sim = sim_matrix_log[i_g, ]
    cur_ind = meta_cell_BSDE[["individual"]]
    cur_pheno = meta_cell_BSDE[[var2test]]
    
    pval_mat["raw", i_g] <- tryCatch({
      cal_w2_pval(
        count_per_gene = cur_sim,
        meta_individual = cur_ind,
        meta_phenotype = cur_pheno,
        perm_num = 999,
        unif_round_unit = 0.5
      )[[1]] # pval
    }, error = function(e) {
      NA
    })
  }

  recal_index <- which(pval_mat["raw", ] < 0.01) # NA values will be excluded
  if(length(recal_index) >0){
    cat("raw p-val calculated. increase the number of permutations in", length(recal_index),"genes\n")
    for (i_g in recal_index) {
      cur_sim = sim_matrix_log[i_g, ]
      cur_ind = meta_cell_BSDE[["individual"]]
      cur_pheno = meta_cell_BSDE[[var2test]]
    
      pval_mat["raw", i_g] <- tryCatch({
        cal_w2_pval( # cal_w2_pval has been parallelized through foreach and dorng on the perm_num
          count_per_gene = cur_sim,
          meta_individual = cur_ind,
          meta_phenotype = cur_pheno,
          perm_num = 9999,
          unif_round_unit = 0.5
        )[[1]] # pval
      }, error = function(e) {
        NA
      })
    }
  }
  
}



# recal_index <- which(pval_mat["raw", ] < 0.001) # NA values will be excluded
# cat("raw p-val calculated. increase the number of permutations in", length(recal_index),"genes/n")
# for (i_g in recal_index) {
#   cur_sim = sim_matrix_log[i_g, ]
#   cur_ind = meta_cell_BSDE[["individual"]]
#   cur_pheno = meta_cell_BSDE[[var2test]]

#   pval_mat["raw", i_g] <- tryCatch({
#     cal_w2_pval( # cal_w2_pval has been parallelized through foreach and dorng on the perm_num
#       count_per_gene = cur_sim,
#       meta_individual = cur_ind,
#       meta_phenotype = cur_pheno,
#       perm_num = 99999,
#       unif_round_unit = 0.5
#     )[[1]] # pval
#   }, error = function(e) {
#     NA
#   })
# }

# pval_mat["adjusted", ] <- p.adjust(pval_mat["raw", ], method = "BH")

pval_temp <- pval_mat["raw", ]
pi0 = 2*sum(pval_temp > 0.5, na.rm = TRUE)/sum(!is.na(pval_temp))
# really need this pi0, otherwise in some cases, error will occur
pi0 = ifelse(pi0>1, 1, ifelse(pi0<0.5, 0.5, pi0)) # 0.5<=pi0<=1
pval_mat["adjusted", ] <- qvalue(pval_temp, pi0=pi0)$qvalues

# # from their simulation
# cur_count = emdbook::rzinbinom(1, sample_mean_k[ig], disp_i[ig], drop_i[ig])
# sim_matrix[ig, idx_i[k]] = cur_count # count
# sim_matrix_log = log2(1 + sim_matrix) # log2 normalized

# op_pval = matrix(ncol = 1, nrow = nrow(sim_matrix_log))

# for (i_g in 1:nrow(sim_matrix_log)) {
#   cur_sim = sim_matrix_log[i_g, ]
#   cur_ind = meta$individual
#   cur_pheno = phenotype

#   op_pval[i_g] = tryCatch({
#     cal_w2_pval( # cal_w2_pval has been parallelized through foreach and dorng on the perm_num
#       count_per_gene = cur_sim,
#       meta_individual = cur_ind,
#       meta_phenotype = cur_pheno,
#       perm_num = perm_num,
#       unif_round_unit = 0.5
#     )[[1]] # pval
#   }, error = function(e) {
#     NA
#   })
# }
# print(date())
# print(gc())
# rownames(op_pval) = gene_id

# # from the vigenette
# count_per_gene=c(rpois(60,6),c(rpois(60,4)))
# meta_individual=paste0("ind",rep(1:12,each=10))
# meta_phenotype=factor(c(rep(1,60),rep(0,60)))
# # show dataset
# table(meta_individual)
# #> meta_individual
# #> ind1 ind10 ind11 ind12 ind2 ind3 ind4 ind5 ind6 ind7 ind8 ind9
# #> 10 10 10 10 10 10 10 10 10 10 10 10
# table(meta_phenotype)
# #> meta_phenotype
# #> 0 1
# #> 60 60
# df=data.frame(count=count_per_gene,ind=meta_individual, diagnosis=meta_phenotype)
# results <- cal_w2_pval(count_per_gene = count_per_gene, meta_individual = meta_individual, 
#                        meta_phenotype = meta_phenotype, perm_num = 200)
# print(results[[1]])
