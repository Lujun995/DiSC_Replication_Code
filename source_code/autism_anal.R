##########################
#parameters to be specified:
#general parameters: WORKINGDIRECTORY, NCORE, DATE, DATADIR, DCADIR
##########################

library(tidyverse)
library(Matrix)
library(parallel)
library(doParallel)
library(DESeq2)
library(ideas)
library(qvalue)
library(data.table)
# rm(list=ls())


setwd("WORKINGDIRECTORY")
source("DiSC.R")

n.cores = NCORE
DDDD = "DATE"
dca_dir = "DCADIR"
data_dir = "DATADIR"

meta_raw = read_tsv(file.path(data_dir, "meta.tsv"), 
                    col_types = "cffcfnffffnnnndd")
SPARSITY_THRESHOLD = 0.8 #or NULL. Keep genes with sparsity < SPARSITY_THRESHOLD
ctypes = c("L2_3", "L4", "Microglia", "Endothelial", "IN-SV2C", "AST-PP", 
           "IN-VIP", "IN-SST", "IN-PV", "AST-FB", "Oligodendrocytes", 
           "L5_6", "L5_6-CC", "OPC", "Neu-NRGN-II", "Neu-NRGN-I", "Neu-mat")

pval_array_list <- vector(mode = "list", length = length(ctypes))
names(pval_array_list) <- ctypes
meths <- c("disc", "deseq", "ideas", "dca_ideas")
ptypes <- c("raw", "adjusted")
var2test <- "diagnosis"
var2adjust <- c("age", "sex", "Seqbatch", "RIN")

registerDoParallel(cores=n.cores)
options(mc.cores=n.cores)
RNGkind("L'Ecuyer-CMRG")

for(ctp in ctypes){
  tryCatch(
    {
      cat("Using cell type:", ctp, "\n")
      ## count_matrix, genes in rows and cells in columns
      count_matrix = readRDS(file.path(sprintf("%sct_mtx/PFC_%s.rds", data_dir, ctp)))
      ngene = nrow(count_matrix)
      ncell = ncol(count_matrix)
      read_depth <- colSums(count_matrix)
      if(ngene < ncell) 
        warning("the number of genes is smaller than the number of cells, please check if genes are in rows")
      cat("ngene: ", ngene, "ncell: ", ncell, "read_depth[1:10]: ", read_depth[1:10], "\n")
      #### Sparsity filtering
      if(!is.null(SPARSITY_THRESHOLD)){
        cat("sparsity filtering with threshold ", SPARSITY_THRESHOLD, " :\n")
        gene2keep = which(rowSums(count_matrix == 0) < SPARSITY_THRESHOLD*ncell)
        count_matrix = count_matrix[gene2keep,]
        read_depth = colSums(count_matrix)
        ngene = nrow(count_matrix)
        ncell = ncol(count_matrix)
        read_depth <- colSums(count_matrix)
        cat("ngene: ", ngene, "ncell: ", ncell, "read_depth[1:10]: ", read_depth[1:10], "\n")
      }
      ## meta data
      meta_cell <- meta_raw %>% filter(region == "PFC") %>%
        filter(cluster == str_replace(ctp, "_", "/")) %>%
        dplyr::rename(`cell_id` = "cell")
      if(sum(!(colnames(count_matrix) %in% meta_cell$cell_id)) >0)
        error("meta_cell and count_matrix dont match") else{
          meta_cell <- meta_cell %>% filter(cell_id %in% colnames(count_matrix))
          if(sum(meta_cell$cell_id != colnames(count_matrix)) >0 )
            error("meta_cell and count_matrix dont match2")
          meta_cell$read_depth <- read_depth
        }
      ## meta_ind
      meta_ind <- meta_cell %>% distinct(individual, .keep_all = TRUE) %>%
        dplyr::select(-`cell_id`, -`genes`, -`UMIs`, -`RNA mitochondr. percent`,
                      -`RNA ribosomal percent`, -`read_depth`) %>%
        dplyr::rename(RIN = "RNA Integrity Number")
      nindi = nrow(meta_ind)
      ## prepossessing
      count_matrix = as.matrix(count_matrix)
      meta_cell = as.data.frame(meta_cell) 
      meta_ind = as.data.frame(meta_ind)
      meta_ind$age <- scale(meta_ind$age)
      meta_ind$RIN <- scale(meta_ind$RIN)
      for(i in 1:ncol(meta_ind)){
        if(is.character(meta_ind[[i]])){
          meta_ind[[i]] = as.factor(meta_ind[[i]])
        }
      }
      meta_ind <- droplevels(meta_ind)
      meta_ind$diagnosis = factor(meta_ind$diagnosis, levels=c("Control", "ASD"))
      pval_array_list[[ctp]]<- 
        array(data = 1,
              dim = c(length(meths), length(ptypes), nrow(count_matrix)),
              dimnames = list(meths, ptypes, rownames(count_matrix)))
      
      ## DiSC
      # discard insufficiently sequenced cells
      depth = min(5000, round(quantile(read_depth, probs = 0.1)))
      cat("read depth filtering with threshold", depth, "\n")
      count_matrix_temp <- count_matrix[, which(read_depth >= depth)]
      meta_cell_temp <- meta_cell %>% 
        filter(cell_id %in% colnames(count_matrix_temp))
      rarefy_mat <- count_matrix_temp %>% t() %>% 
        GUniFrac::Rarefy(.) %>%
        magrittr::extract2(., "otu.tab.rff") %>%
        t()
      
      obj <- DiSC(ct.mat = rarefy_mat, 
                  cell.ind = meta_cell_temp, 
                  metadata = meta_ind,
                  outcome = var2test, 
                  covariates = var2adjust, 
                  cell.id = "cell_id",
                  individual.id = "individual",
                  verbose = FALSE)
      # list(meths, ptypes, ctypes)
      pval_array_list[[ctp]]["disc", "adjusted", ] <- obj$p.adj.fdr
      pval_array_list[[ctp]]["disc", "raw", ] <- obj$p.raw
      
      
      #ideas
      dist1 = ideas_dist(count_input = count_matrix, 
                         meta_cell = meta_cell, meta_ind = meta_ind, 
                         var_per_cell = "read_depth", var2test = var2test, 
                         var2test_type = "binary",
                         d_metric = "Was", fit_method = "nb")
      ideas_pval_raw <- permanova(dist_array = dist1, meta_ind = meta_ind, 
                                  var2test = var2test, var2adjust = var2adjust, 
                                  var2test_type = "binary", n_perm=999)
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
                                      var2test_type = "binary", n_perm=99999)
        ideas_pval_raw[recal_index] <- ideas_pval_raw_r
      }
      pval_temp <- ideas_pval_raw
      pi0 = 2*sum(pval_temp > 0.5, na.rm = TRUE)/sum(!is.na(pval_temp))
      # really need this pi0, otherwise in some cases, error will occur
      pi0 = ifelse(pi0>1, 1, ifelse(pi0<0.5, 0.5, pi0)) # 0.5<=pi0<=1
      ideas_pval_adj <- qvalue(pval_temp, pi0=pi0)$qvalues
      
      pval_array_list[[ctp]]["ideas", "adjusted", ] <- ideas_pval_adj
      pval_array_list[[ctp]]["ideas", "raw", ] <- ideas_pval_raw
      
      # dca-ideas
      f_mean = file.path(dca_dir, paste0(ctp, "_mean_norm.tsv"))
      f_disp = file.path(dca_dir, paste0(ctp, "_dispersion.tsv"))
      f_pi   = file.path(dca_dir, paste0(ctp, "_pi.tsv"))
      
      dca_mean = fread(f_mean, sep="\t", data.table = FALSE)
      dca_disp = fread(f_disp)
      dca_pi   = fread(f_pi)
      
      if(sum(!(meta_cell$cell_id == colnames(dca_mean)[-1])) != 0 | 
         sum(!(rownames(count_matrix) %in% dca_mean$V1)) != 0 )
        stop("dca_mean does not match meta_cell or count_matrix")
      
      w2kp = match(rownames(count_matrix), dca_mean$V1)
      dca_mean = dca_mean[w2kp,]
      dca_disp = dca_disp[w2kp,]
      dca_pi   = dca_pi[w2kp,]
      
      rownames(dca_mean) = dca_mean$V1
      rownames(dca_disp) = dca_disp$V1
      rownames(dca_pi)   = dca_pi$V1
      
      dca_mean = data.matrix(dca_mean[,-1])
      dca_disp = data.matrix(dca_disp[,-1, with=FALSE])
      dca_pi   = data.matrix(dca_pi[,-1,   with=FALSE])
      dca_par_list = list(dca_mean, dca_disp, dca_pi)
      
      dist1 = ideas_dist(count_input = dca_par_list, 
                         meta_cell = meta_cell, meta_ind = meta_ind, 
                         var_per_cell = "read_depth", var2test = var2test, 
                         var2test_type = "binary",
                         d_metric = "Was", fit_method = "dca_direct")
      dca_ideas_pval_raw <- permanova(dist_array = dist1, meta_ind = meta_ind, 
                                      var2test = var2test, var2adjust = var2adjust, 
                                      var2test_type = "binary", n_perm=999)
      recal_index <- which(dca_ideas_pval_raw < 0.01)# NA values will be excluded
      if(length(recal_index) > 0){
        if (length(recal_index) == 1){
          slice_dist = dist1[recal_index, , ]  
          recal_dist = array(dim = c(1, dim(slice_dist)))
          recal_dist[1, , ] = slice_dist
        }else{
          recal_dist = dist1[recal_index, , ]
        }
        dca_ideas_pval_raw_r <- permanova(dist_array = recal_dist, meta_ind = meta_ind, 
                                          var2test = var2test, var2adjust = var2adjust, 
                                          var2test_type = "binary", n_perm=99999)
        dca_ideas_pval_raw[recal_index] <- dca_ideas_pval_raw_r
      }
      pval_temp <- dca_ideas_pval_raw
      pi0 = 2*sum(pval_temp > 0.5, na.rm = TRUE)/sum(!is.na(pval_temp))
      # really need this pi0, otherwise in some cases, error will occur
      pi0 = ifelse(pi0>1, 1, ifelse(pi0<0.5, 0.5, pi0)) # 0.5<=pi0<=1
      dca_ideas_pval_adj <- qvalue(pval_temp, pi0=pi0)$qvalues
      
      pval_array_list[[ctp]]["dca_ideas", "adjusted", ] <- dca_ideas_pval_adj
      pval_array_list[[ctp]]["dca_ideas", "raw", ] <- dca_ideas_pval_raw
      
      #pseudo_bulk sample method, DEseq2
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
      
      dds <- DESeqDataSetFromMatrix(countData = count_matrix_bulk,
                                    colData = meta_ind,
                                    design = ~ age + sex + Seqbatch + RIN + diagnosis)
      dds <- DESeq(dds, quiet = TRUE)
      pval_array_list[[ctp]]["deseq", "adjusted", ] <- 
        results(dds, name = "diagnosis_ASD_vs_Control")$padj
      pval_array_list[[ctp]]["deseq", "raw", ] <- 
        results(dds, name = "diagnosis_ASD_vs_Control")$pvalue
    }, error = function(err) {
      print(err)
      save.image(str_glue("unfinished_autism_{ctp}_{DDDD}.RData"))
    }) #try-catch part
  # when error occurs, save the current working space and go to the next ctp
  
}# go to the next ctp

stopImplicitCluster()
save(SPARSITY_THRESHOLD, pval_array_list, 
     file = str_glue("results_autism_pfc_{DDDD}.RData"))



