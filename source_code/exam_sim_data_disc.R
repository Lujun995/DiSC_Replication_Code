##########################
#parameters to be specified:
#general parameters: WORKINGDIRECTORY, NCORE, NLIM_1, NLIM_2, 
#PARA_OUT, SIMU_D, CAT
#parameters in simulations: DIFF, NCELL, NPG, VNCELL
##########################

setwd("WORKINGDIRECTORY")
library(tidyverse)
library(parallel)

# Need to specify the method when creating this file
method = "disc" # a single value
# be careful, the script assumes to save/return a value called "{method}_pval"
library(Matrix) # specifying package dependencies for the "method"
source("DiSC.R")
# library(doParallel)
# library(ideas)
# library(qvalue)

set.seed(seed = 12345181)
n_cores = NCORE
n_lim_1 = NLIM_1
n_lim_2 = NLIM_2
DATE = "SIMU_D"# the date when simulation datasets were created
para_out = "PARA_OUT"
catgr = "CAT"

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
    
    #############DiSC start
    cl <- makeCluster(n_cores)
    clusterExport(cl, c("na.pad", "perm.fdr.adj", "perm.fwer.adj",
                        "DiSC", "basic.func", "conf"))
    clusterEvalQ(cl, {
      library(matrixStats)
      library(tidyverse)
      setwd("WORKINGDIRECTORY")
    })
    disc_pval <- parSapply(cl, 1:boot, simplify = FALSE, FUN = function(bbbb){
      load(str_glue("{conf}_b{bbbb}.RData"))
      #str(dat_list)
      # nGeneTotal = sum(unlist(lapply(dat_list[["gene_index"]], FUN = length)))
      meta_ind <- dat_list[["meta_ind"]]
      meta_cell <- dat_list[["meta_cell"]]
      count_matrix <- dat_list[["count_matrix"]]
      meta_ind$phenotype <- as.factor(meta_ind$phenotype)
      meta_cell$phenotype <- as.factor(meta_cell$phenotype)
      
      
      obj <- DiSC(data.mat = count_matrix,
                  cell.ind = meta_cell,
                  metadata = meta_ind,
                  outcome = "phenotype", covariates = "RIN",
                  cell.id = "cell_id",
                  individual.id = "individual", 
                  verbose = FALSE, sequencing.data = TRUE)
      disc_pval_adj <- obj$p.adj.fdr
      disc_pval_raw <- obj$p.raw
      disc_pval <- matrix(c(disc_pval_raw, disc_pval_adj),
                          ncol = 2, byrow = FALSE)
      colnames(disc_pval) <- c("raw", "adjusted")
      return(disc_pval)
    })
    stopCluster(cl)
    ############DiSC completed
    
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


