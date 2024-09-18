##########################
#parameters to be specified:
#general parameters: /home/guanwh/zhan7474/mayo/para_data_sim/, 80, 0, 12.1, 
#simulation_parameters_0811.txt, 0811, T1
#parameters in simulations: DIFF, NCELL, NPG, VNCELL
##########################

setwd("/home/guanwh/zhan7474/mayo/para_data_sim/")
library(tidyverse)
library(parallel)

# Need to specify the method when creating this file
method = "deseq" # a single value
# be careful, the script assumes to save/return a value called "{method}_pval"
require(DESeq2) # specifying package dependencies for the "method"
# source("DiSC.R")
# library(doParallel)
# library(ideas)
# library(qvalue)

set.seed(seed = 12345181)
n_cores = 80
n_lim_1 = 0
n_lim_2 = 12.1
DATE = "0811"# the date when simulation datasets were created
para_out = "simulation_parameters_0811.txt"
catgr = "T1"

# read in parameters
parameters <- read_file(para_out) %>%
  str_split_1(pattern = ":") 
parameters <- parameters[parameters != ""]
parameters <- parameters[!duplicated(parameters)]
n_per_group_vec <- parameters %>%
  str_extract(pattern = "NPG[:blank:]*=[:blank:]*[:digit:]+") %>%
  str_remove(pattern = "NPG[:blank:]*=[:blank:]*") %>%
  as.numeric()
diff_vec <- parameters %>%
  str_extract(pattern = "DIFF[:blank:]*=[:blank:]*[[:digit:]]+\\.*[[:digit:]]*") %>%
  str_remove(pattern = "DIFF[:blank:]*=[:blank:]*") %>%
  as.numeric()
if(catgr == "T1"){
  parameters <- parameters[diff_vec == 1] # for type I error and FDR
}else{
  parameters <- parameters[n_per_group_vec >= n_lim_1 & 
                             n_per_group_vec < n_lim_2 &
                             diff_vec != 1]
}


tryCatch({
  for(paras in parameters){
    eval(parse(text=paras)) # VNCELL=FALSE;DIFF=1.1;NCELL=375;NPG=12;BOOT=50
    boot = BOOT
    diff = exp(abs(log(DIFF))) #promise diff >= 1
    n_cas = NPG
    n_ctr = NPG
    nall = n_cas + n_ctr
    neq_n_cell = VNCELL
    n_cell= NCELL
    
    conf <- str_glue("n{n_cas}_{n_ctr}_d{diff}_c{n_cell}_v{neq_n_cell}_{DATE}")
    if(file.exists(str_glue("results_{method}_{conf}.RData"))){
      cat(str_glue("finished results_{method}_{conf}.RData\n"))
      next
    }
    if(file.exists(str_glue("unfinished_{method}_{conf}.RData"))){
      # exception handling
      # to resume the previous work
      source("exam_sim_data_resume_disc_0805.R")
      save(list = str_glue("{method}_pval"),
           file = str_glue("results_{method}_{conf}.RData"))
      file.remove(str_glue("unfinished_{method}_{conf}.RData"))
      cat(str_glue("finished results_{method}_{conf}.RData\n"))
      next
    }
    
    ############pseudo-bulk method, DESeq2 starts
    cl <- makeCluster(n_cores)
    clusterExport(cl, c("conf"))
    clusterEvalQ(cl, {
      library(stringr)
      library(DESeq2)
      setwd("/home/guanwh/zhan7474/mayo/para_data_sim/")
    })
    deseq_pval <- parSapply(cl, 1:boot, simplify = FALSE, FUN = function(bbbb){
      load(str_glue("{conf}_b{bbbb}.RData"))
      #str(dat_list)
      # nGeneTotal = sum(unlist(lapply(dat_list[["gene_index"]], FUN = length)))
      meta_ind <- dat_list[["meta_ind"]]
      meta_cell <- dat_list[["meta_cell"]]
      count_matrix <- dat_list[["count_matrix"]]
      meta_ind$phenotype <- as.factor(meta_ind$phenotype)
      meta_cell$phenotype <- as.factor(meta_cell$phenotype)
      
      # prepossing
      count_matrix_bulk = matrix(NA, nrow = nrow(count_matrix), 
                                 ncol = nrow(meta_ind))
      rownames(count_matrix_bulk) = rownames(count_matrix)
      colnames(count_matrix_bulk) = meta_ind[["individual"]]
      if(sum(meta_cell[["cell_id"]] != colnames(count_matrix))==0)
        for (i_ind in 1:length(meta_ind[["individual"]])) {
          cur_ind   = meta_ind[["individual"]][i_ind]
          cur_ind_m = count_matrix[, meta_cell[["individual"]] == cur_ind] 
          count_matrix_bulk[, i_ind] = rowSums(cur_ind_m, na.rm = TRUE)
        }else
          stop("meta_cell$cell_id and colnames(count_matrix) dont match")
      
      dds <- DESeq2::DESeqDataSetFromMatrix(countData = count_matrix_bulk,
                                            colData = meta_ind,
                                            design = ~ RIN + phenotype)
      dds <- DESeq2::DESeq(dds, quiet = TRUE)
      deseq_pval_adj <- DESeq2::results(dds, name = "phenotype_1_vs_0")$padj
      deseq_pval_adj[is.na(deseq_pval_adj)] <- 1 
      # when allowing independentFiltering and cooksCutoff, p-values for some
      # genes may be returned as NA. Set them to an insignificant p-value.
      deseq_pval_raw <- DESeq2::results(dds, name = "phenotype_1_vs_0")$pvalue
      deseq_pval_raw[is.na(deseq_pval_raw)] <- 1
      # the default method for adjusting for multiplicity is BH in DESeq
      #When disabling independentFiltering and cooksCutoff:
      #identical(p.adjust(deseq_pval_raw, method = "BH"), deseq_pval_adj)#TRUE
      deseq_pval <- matrix(c(deseq_pval_raw, deseq_pval_adj),
                           ncol = 2, byrow = FALSE)
      colnames(deseq_pval) <- c("raw", "adjusted")
      
      return(deseq_pval)
    })
    stopCluster(cl)
    ############DESeq2 completed
    
    save(list = str_glue("{method}_pval"),
         file = str_glue("results_{method}_{conf}.RData"))
    cat(str_glue("finished results_{method}_{conf}.RData\n"))
  }
}, error=function(err){
  # exception handling. See the part above
  # if error occurs, save current work image and restart outside
  print(err)
  print("Need to be restarted...")
  save.image(file = str_glue("unfinished_{method}_{conf}.RData"))
  q(save = "no")
  return(0)
})# try-catach part

q(save = "no")


