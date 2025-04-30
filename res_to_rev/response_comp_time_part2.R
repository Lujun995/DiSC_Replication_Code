
setwd("/home/guanwh/zhan7474/mayo/res_to_rev/")

library(tidyverse)
library(Matrix)
library(parallel)
library(doParallel)
library(DESeq2)
library(ideas)
library(qvalue)
library(reticulate)
library(microbenchmark)
use_python("/home/guanwh/zhan7474/.conda/envs/r_env/bin/python", required = TRUE)
use_condaenv("/home/guanwh/zhan7474/.conda/envs/r_env", required = TRUE)
library(BSDE)
library(iDESC)
source("DiSC.R") # DiSC
source("zp_prediction.r") # iDESC
source("iDESC.r") # iDESC
source("cal_w2_pval.r") # BSDE

DiSC_anal <- function(count_matrix, meta_cell, meta_ind, var2test, 
                      var2adjust, perm.no = 999){
  suppressMessages(suppressWarnings({
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
                verbose = FALSE, perm.no = perm.no, sequencing.data = TRUE)
    pval_mat["adjusted", ] <- obj$p.adj.fdr
    pval_mat["raw", ] <- obj$p.raw
    return(pval_mat)
  }))
  
}

IDEAS_anal <- function(count_matrix, meta_cell, meta_ind, var2test, 
                       var2adjust, var2test_type, recal = TRUE){
  suppressMessages(suppressWarnings({
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
    if(recal){
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
    }
    
    pval_temp <- ideas_pval_raw
    pi0 = 2*sum(pval_temp > 0.5, na.rm = TRUE)/sum(!is.na(pval_temp))
    # really need this pi0, otherwise in some cases, error will occur
    pi0 = ifelse(pi0>1, 1, ifelse(pi0<0.5, 0.5, pi0)) # 0.5<=pi0<=1
    ideas_pval_adj <- qvalue(pval_temp, pi0=pi0)$qvalues
    
    pval_mat["adjusted", ] <- ideas_pval_adj
    pval_mat["raw", ] <- ideas_pval_raw
    return(pval_mat)
  }))
  
}

DESeq2_anal <- function(count_matrix, meta_cell, meta_ind, var2test, 
                        var2adjust){
  suppressMessages(suppressWarnings({
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
    return(pval_mat)
  }))
  
}

iDESC_anal <- function(count_matrix, meta_cell, meta_ind, var2test,
                       cl = NULL){
  suppressMessages(suppressWarnings({
    pval_mat <- matrix(data = 1, nrow = 2, ncol = nrow(count_matrix), 
                       dimnames = list(c("adjusted", "raw"),
                                       rownames(count_matrix)))
    meta_cell_iDESC <- meta_cell[, c("cell_id", "individual")]
    meta_cell_iDESC <- merge(x = meta_cell_iDESC, y = meta_ind, by = "individual", 
                             all.x = TRUE, sort = FALSE)
    stopifnot(colnames(count_matrix) == meta_cell_iDESC$cell_id)
    rownames(meta_cell_iDESC) <- meta_cell_iDESC$cell_id
    
    # if(exists("cl")){
    #   if(inherits(cl, "cluster"))
    if(!is.null(cl) & inherits(cl, "cluster")){
      clusterEvalQ(cl, {
        library(iDESC)
        library(stringr)
      })
      result <- iDESC(count_matrix, meta_cell_iDESC, 
                      subject_var="individual", group_var=var2test,
                      norm_opt="SeqDepth", span = 0.7, # according to "function usage"
                      sub_cell_filtering = 0, gene_sub_filtering = 0, 
                      gene_cell_filtering = 0, ncell_filtering = 0, 
                      cl = cl)
      pval_mat["raw", ] <- result[[stringr::str_subset(names(result), "Pval_Beta")]]
      #}
    } else {
      result <- iDESC(count_matrix, meta_cell_iDESC, 
                      subject_var="individual", group_var=var2test,
                      norm_opt="SeqDepth", span = 0.7, # according to "function usage"
                      sub_cell_filtering = 0, gene_sub_filtering = 0, 
                      gene_cell_filtering = 0, ncell_filtering = 0)
      pval_mat["raw", ] <- result[[stringr::str_subset(names(result), "Pval_Beta")]]
    }
    pval_temp <- pval_mat["raw", ]
    pi0 = 2*sum(pval_temp > 0.5, na.rm = TRUE)/sum(!is.na(pval_temp))
    # really need this pi0, otherwise in some cases, error will occur
    pi0 = ifelse(pi0>1, 1, ifelse(pi0<0.5, 0.5, pi0)) # 0.5<=pi0<=1
    pval_mat["adjusted", ] <- qvalue::qvalue(pval_temp, pi0=pi0)$qvalues
    # pval_mat["adjusted", ] <- p.adjust(pval_mat["raw", ], method = "BH")
    return(pval_mat)
  }))
}

BSDE_anal <- function(count_matrix, meta_cell, meta_ind, var2test,
                      cl = NULL, recal = TRUE){
  suppressMessages(suppressWarnings({
    pval_mat <- matrix(data = 1, nrow = 2, ncol = nrow(count_matrix), 
                       dimnames = list(c("adjusted", "raw"),
                                       rownames(count_matrix)))
    
    meta_cell_BSDE <- meta_cell[, c("cell_id", "individual")]
    meta_cell_BSDE <- merge(x = meta_cell_BSDE, y = meta_ind, by = "individual", 
                            all.x = TRUE, sort = FALSE)
    if(inherits(meta_cell_BSDE[[var2test]], "factor"))
      meta_cell_BSDE[[var2test]] <- as.numeric(meta_cell_BSDE[[var2test]]) - 1 # convert to 0/1
    stopifnot(inherits(meta_cell_BSDE[[var2test]], "numeric"))
    stopifnot(unique(meta_cell_BSDE[[var2test]]) %in% c(0,1))
    stopifnot(colnames(count_matrix) == meta_cell_BSDE$cell_id)
    
    sim_matrix_log = log2(1 + count_matrix)
    
    # if(exists("cl")){
    #   if(inherits(cl, "cluster")){
    if(!is.null(cl) & inherits(cl, "cluster")){
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
      
      if(recal){
        recal_index <- which(pval_mat["raw", ] < 0.01) # NA values will be excluded
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
    } else {
      for (i_g in 1:nrow(count_matrix)) {
        # if(i_g %% 1000 == 0)  cat(i_g, "\n")
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
      
      if(recal){
        recal_index <- which(pval_mat["raw", ] < 0.01) # NA values will be excluded
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
    pval_temp <- pval_mat["raw", ]
    pi0 = 2*sum(pval_temp > 0.5, na.rm = TRUE)/sum(!is.na(pval_temp))
    # really need this pi0, otherwise in some cases, error will occur
    pi0 = ifelse(pi0>1, 1, ifelse(pi0<0.5, 0.5, pi0)) # 0.5<=pi0<=1
    pval_mat["adjusted", ] <- qvalue(pval_temp, pi0=pi0)$qvalues
    return(pval_mat)
  }))
}

setups <- function(dat_list, ngene = 700, ncell = 100, nsample = 5){
  count_matrix <- dat_list$count_matrix
  gene_index <- dat_list$gene_index
  meta_ind <- dat_list$meta_ind
  meta_cell <- dat_list$meta_cell
  
  count_matrix <- count_matrix[1:ngene,]
  gene_index <- lapply(gene_index, function(gs) gs[gs <= ngene])
  meta_cell$cell_rd <- colSums(count_matrix)
  
  meta_cell <- meta_cell %>% group_by(individual) %>% slice_head(n=ncell)
  count_matrix <- count_matrix[, meta_cell$cell_id]
  
  meta_ind <- meta_ind %>% group_by(phenotype) %>% slice_head(n=nsample)
  meta_cell <- meta_cell %>% filter(individual %in% meta_ind$individual)
  count_matrix <- count_matrix[, meta_cell$cell_id]
  
  return(list(count_matrix = count_matrix, gene_index = gene_index,
              meta_ind = meta_ind, meta_cell = meta_cell))
}

system("lscpu | grep 'Model name:'")
system("lscpu | grep 'CPU MHz:'")

registerDoParallel(cores=1)
options(mc.cores=1)
RNGkind("L'Ecuyer-CMRG")

set.seed(seed = 123456)
load("/home/guanwh/zhan7474/mayo/res_to_rev/n200_200_d1.1_c1000_vFALSE_0811_b1.RData")

boot = 2
NGENES <- c(100, 200, 300, 400)
NCELL <- c(100, 200, 300, 400)
NSAMPLE <- c(5, 10, 15, 20)
parameters <- bind_rows(
  data.frame(nc = NCELL[1], ns = NSAMPLE[1], ng = NGENES),
  data.frame(nc = NCELL[1], ns = NSAMPLE, ng = NGENES[1]),
  data.frame(nc = NCELL, ns = NSAMPLE[1], ng = NGENES[1])
) %>% distinct(.keep_all = TRUE)

results <- bind_rows(sapply(1:nrow(parameters), simplify = FALSE, FUN = \(rr){
  cat("current iteration:", rr, "\n")
  nc = parameters[rr, "nc"]
  ns = parameters[rr, "ns"]
  ng = parameters[rr, "ng"]
  dat_list_temp <- setups(dat_list, ngene = ng, ncell = nc, nsample = ns)
  count_matrix = dat_list_temp$count_matrix
  meta_cell = dat_list_temp$meta_cell
  meta_cell$phenotype = factor(ifelse(meta_cell$phenotype == 1, "case", "control"))
  names(meta_cell)[which(names(meta_cell) == "cell_rd")] <- "read_depth"
  meta_ind = dat_list_temp$meta_ind
  meta_ind$phenotype = factor(ifelse(meta_ind$phenotype == 1, "case", "control"))
  var2test = "phenotype"
  var2adjust = NULL
  var2test_type = "binary"
  res <- microbenchmark(
    DiSC_anal(count_matrix, meta_cell, meta_ind, var2test, var2adjust),
    IDEAS_anal(count_matrix, meta_cell, meta_ind, var2test, var2adjust, var2test_type, recal = FALSE),
    DESeq2_anal(count_matrix, meta_cell, meta_ind, var2test, var2adjust),
    iDESC_anal(count_matrix, meta_cell, meta_ind, var2test, cl = NULL),
    BSDE_anal(count_matrix, meta_cell, meta_ind, var2test, cl = NULL, recal = FALSE),
    times = c(boot*100, boot, boot*100, boot, boot))
  res$ngene = ng
  res$ncell = nc
  res$nsample = ns
  return(as.data.frame(res))
}))
saveRDS(object = results, file = "time_results_0107.rds")

# typical datasets
dat_list_temp <- setups(dat_list, ngene = 2000, ncell = 375, nsample = 12)
count_matrix = dat_list_temp$count_matrix
meta_cell = dat_list_temp$meta_cell
meta_cell$phenotype = factor(ifelse(meta_cell$phenotype == 1, "case", "control"))
names(meta_cell)[which(names(meta_cell) == "cell_rd")] <- "read_depth"
meta_ind = dat_list_temp$meta_ind
meta_ind$phenotype = factor(ifelse(meta_ind$phenotype == 1, "case", "control"))
var2test = "phenotype"
var2adjust = NULL
var2test_type = "binary"
results <- microbenchmark(
  DiSC_anal(count_matrix, meta_cell, meta_ind, var2test, var2adjust, perm.no = 999),
  DiSC_anal(count_matrix, meta_cell, meta_ind, var2test, var2adjust, perm.no = 99),
  DESeq2_anal(count_matrix, meta_cell, meta_ind, var2test, var2adjust),
  IDEAS_anal(count_matrix, meta_cell, meta_ind, var2test, var2adjust, var2test_type, recal = FALSE),
  times = c(10,10,10,3))
saveRDS(object = as.data.frame(results), 
        file = "time_results_ng2000_nc375_ns24_0107.rds")

# dat_list_temp <- setups(dat_list, ngene = 8000, ncell = 375, nsample = 12)
# count_matrix = dat_list_temp$count_matrix
# meta_cell = dat_list_temp$meta_cell
# meta_cell$phenotype = factor(ifelse(meta_cell$phenotype == 1, "case", "control"))
# names(meta_cell)[which(names(meta_cell) == "cell_rd")] <- "read_depth"
# meta_ind = dat_list_temp$meta_ind
# meta_ind$phenotype = factor(ifelse(meta_ind$phenotype == 1, "case", "control"))
# var2test = "phenotype"
# var2adjust = NULL
# var2test_type = "binary"
# results <- microbenchmark(
#   DiSC_anal(count_matrix, meta_cell, meta_ind, var2test, var2adjust),
#   DESeq2_anal(count_matrix, meta_cell, meta_ind, var2test, var2adjust),
#   IDEAS_anal(count_matrix, meta_cell, meta_ind, var2test, var2adjust, var2test_type, recal = FALSE),
#   times = 3)
# saveRDS(object = as.data.frame(results), 
#         file = "time_results_ng8000_nc375_ns24_0107.rds")

# dat_list_temp <- setups(dat_list, ngene = 8000, ncell = 1000, nsample = 50)
# count_matrix = dat_list_temp$count_matrix
# meta_cell = dat_list_temp$meta_cell
# meta_cell$phenotype = factor(ifelse(meta_cell$phenotype == 1, "case", "control"))
# names(meta_cell)[which(names(meta_cell) == "cell_rd")] <- "read_depth"
# meta_ind = dat_list_temp$meta_ind
# meta_ind$phenotype = factor(ifelse(meta_ind$phenotype == 1, "case", "control"))
# var2test = "phenotype"
# var2adjust = NULL
# var2test_type = "binary"
# results <- microbenchmark(
#   DiSC_anal(count_matrix, meta_cell, meta_ind, var2test, var2adjust),
#   DESeq2_anal(count_matrix, meta_cell, meta_ind, var2test, var2adjust),
#   times = 3)
# saveRDS(object = as.data.frame(results), 
#         file = "time_results_ng8000_nc1000_ns100_0107_2.rds")

dat_list_temp <- dat_list
count_matrix = dat_list_temp$count_matrix
meta_cell = dat_list_temp$meta_cell
meta_cell$phenotype = factor(ifelse(meta_cell$phenotype == 1, "case", "control"))
names(meta_cell)[which(names(meta_cell) == "cell_rd")] <- "read_depth"
meta_ind = dat_list_temp$meta_ind
meta_ind$phenotype = factor(ifelse(meta_ind$phenotype == 1, "case", "control"))
var2test = "phenotype"
var2adjust = NULL
var2test_type = "binary"
results <- microbenchmark(
  DiSC_anal(count_matrix, meta_cell, meta_ind, var2test, var2adjust, perm.no = 999),
  # DiSC_anal(count_matrix, meta_cell, meta_ind, var2test, var2adjust, perm.no = 99),
  DESeq2_anal(count_matrix, meta_cell, meta_ind, var2test, var2adjust),
  # IDEAS_anal(count_matrix, meta_cell, meta_ind, var2test, var2adjust, var2test_type, recal = FALSE)
  times = 3)
saveRDS(object = as.data.frame(results), 
        file = "time_results_ng8000_nc1000_ns400_0107.rds")

stopImplicitCluster()