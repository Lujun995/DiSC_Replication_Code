##########################
#parameters to be specified:
#general parameters: /home/guanwh/zhan7474/mayo/para_data_sim/, 30, 20.1, 35.1, 
#simulation_parameters_0811.txt, 0811, 3
#parameters in simulations: DIFF, NCELL, NPG, VNCELL
##########################

setwd("/home/guanwh/zhan7474/mayo/para_data_sim/")
library(tidyverse)
library(parallel)

# Need to specify the method when creating this file
method = "ideas" # a single value
# be careful, the script assumes to save/return a value called "{method}_pval"
library(ideas) # specifying package dependencies for the "method"
library(qvalue)
library(doParallel)
# source("DiSC.R")
# library(doParallel)
# library(ideas)
# library(qvalue)

set.seed(seed = 12345181)
n_cores = 30
n_lim_1 = 20.1
n_lim_2 = 35.1
DATE = "0811"# the date when simulation datasets were created
para_out = "simulation_parameters_0811.txt"
catgr = "3"

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
    
    ############ideas starts
    ideas_pval <- vector(mode = "list", length = boot)
    registerDoParallel(cores=n_cores)
    options(mc.cores=n_cores)
    RNGkind("L'Ecuyer-CMRG")
    for(bbbb in 1:boot){
      tryCatch(
        {
          load(str_glue("{conf}_b{bbbb}.RData"))
          #str(dat_list)
          # nGeneTotal = sum(unlist(lapply(dat_list[["gene_index"]], FUN = length)))
          meta_ind <- dat_list[["meta_ind"]]
          meta_cell <- dat_list[["meta_cell"]]
          count_matrix <- dat_list[["count_matrix"]]
          meta_ind$phenotype <- as.factor(meta_ind$phenotype)
          meta_cell$phenotype <- as.factor(meta_cell$phenotype)
          
          # rate-limiting step
          dist1 <- ideas_dist(count_input = count_matrix, 
                              meta_cell = meta_cell, meta_ind = meta_ind, 
                              var_per_cell = "cell_rd", var2test = "phenotype", 
                              var2test_type = "binary",
                              d_metric = "Was", fit_method = "nb")
          ideas_pval_raw <- permanova(dist_array = dist1, meta_ind = meta_ind, 
                                      var2test = "phenotype", var2adjust = "RIN", 
                                      var2test_type = "binary", n_perm=999)
          recal_index <- which(ideas_pval_raw < 0.01)# NA values will be excluded
          if(length(recal_index) > 0){
            if (length(recal_index) == 1){
              slice_dist <- dist1[recal_index, , ]  
              recal_dist <- array(dim = c(1, dim(slice_dist)))
              recal_dist[1, , ] <- slice_dist
            }else{
              recal_dist <- dist1[recal_index, , ]
            }
            # diff = exp(abs(log(DIFF))) #promise diff >= 1
            # n_cas = NPG
            # n_ctr = NPG
            # nall = n_cas + n_ctr
            # neq_n_cell = VNCELL
            # n_cell= NCELL
            ideas_pval_raw_r <- permanova(dist_array = recal_dist, 
                                          meta_ind = meta_ind, 
                                          var2test = "phenotype", 
                                          var2adjust = "RIN", 
                                          var2test_type = "binary", 
                                          n_perm=ifelse(diff<=1.2|nall<=10|n_cell<=250, 
                                                        99999, 9999))
            # when the power is too low, we need more permutations for q-value
            ideas_pval_raw[recal_index] <- ideas_pval_raw_r
          }
          
          pval_temp <- ideas_pval_raw
          pi0 = 2*sum(pval_temp > 0.5, na.rm = TRUE)/sum(!is.na(pval_temp))
          # really need this pi0, otherwise in some cases, error will occur
          pi0 = ifelse(pi0>1, 1, ifelse(pi0<0.5, 0.5, pi0)) # 0.5<=pi0<=1
          ideas_pval_adj <- qvalue(pval_temp, pi0=pi0)$qvalues
          ideas_pval[[bbbb]]  <- matrix(c(ideas_pval_raw, ideas_pval_adj),
                                        ncol = 2, byrow = FALSE, 
                                        dimnames = list(NULL, 
                                                        c("raw", "adjusted")))
          gc()
        },
        error = function(err){
          #due to the terrible garbage collection features in R, the pro-
          #gram may stop occasionally and irreproducibly.
          #in this case, try to restart the P-value calculation and resume the job
          print(err)
          print("Try to restart...")
          stopImplicitCluster()
          
          registerDoParallel(cores=n_cores)
          # it seems the nodes registered inside the function are still accessible
          #outside
          options(mc.cores=n_cores)
          RNGkind("L'Ecuyer-CMRG")
          
          load(str_glue("{conf}_b{bbbb}.RData"))
          #str(dat_list)
          # nGeneTotal = sum(unlist(lapply(dat_list[["gene_index"]], FUN = length)))
          meta_ind <- dat_list[["meta_ind"]]
          meta_cell <- dat_list[["meta_cell"]]
          count_matrix <- dat_list[["count_matrix"]]
          meta_ind$phenotype <- as.factor(meta_ind$phenotype)
          meta_cell$phenotype <- as.factor(meta_cell$phenotype)
          
          # rate-limiting step
          dist1 <- ideas_dist(count_input = count_matrix, 
                              meta_cell = meta_cell, meta_ind = meta_ind, 
                              var_per_cell = "cell_rd", var2test = "phenotype", 
                              var2test_type = "binary",
                              d_metric = "Was", fit_method = "nb")
          ideas_pval_raw <- permanova(dist_array = dist1, meta_ind = meta_ind, 
                                      var2test = "phenotype", var2adjust = "RIN", 
                                      var2test_type = "binary", n_perm=999)
          recal_index <- which(ideas_pval_raw < 0.01)# NA values will be excluded
          if(length(recal_index) > 0){
            if (length(recal_index) == 1){
              slice_dist <- dist1[recal_index, , ]  
              recal_dist <- array(dim = c(1, dim(slice_dist)))
              recal_dist[1, , ] <- slice_dist
            }else{
              recal_dist <- dist1[recal_index, , ]
            }
            ideas_pval_raw_r <- permanova(dist_array = recal_dist, 
                                          meta_ind = meta_ind, 
                                          var2test = "phenotype", 
                                          var2adjust = "RIN", 
                                          var2test_type = "binary", n_perm=9999)
            ideas_pval_raw[recal_index] <- ideas_pval_raw_r
          }
          
          pval_temp <- ideas_pval_raw
          pi0 = 2*sum(pval_temp > 0.5, na.rm = TRUE)/sum(!is.na(pval_temp))
          # really need this pi0, otherwise in some cases, error will occur
          pi0 = ifelse(pi0>1, 1, ifelse(pi0<0.5, 0.5, pi0)) # 0.5<=pi0<=1
          ideas_pval_adj <- qvalue(pval_temp, pi0=pi0)$qvalues
          
          ideas_pval[[bbbb]] <<- matrix(c(ideas_pval_raw, ideas_pval_adj),
                                        ncol = 2, byrow = FALSE, 
                                        dimnames = list(NULL, 
                                                        c("raw", "adjusted")))
          gc()
          
          print("Restarted successfully. Continue...")
          
          return(0)
        }
      )
    }
    stopImplicitCluster()
    ############IDEAS completed
    
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


