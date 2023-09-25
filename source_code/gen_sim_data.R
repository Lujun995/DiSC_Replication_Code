##########################
#parameters to be specified:
#genral parameters: WORKINGDIRECTORY, NCORE, DATE, PARA_OUT, EST_OUT
#parameters in model: DIFF, NPG, NCELL, VNCELL, (also BOOT)
##########################
library(tidyverse)
# library(emdbook)
library(MASS)
library(doParallel)
library(foreach)
library(doRNG)

setwd("WORKINGDIRECTORY")
est_out = "EST_OUT"
para_out = "PARA_OUT"
dddd = "DATE"
load(est_out) #from parameters estimation step 
# save(list = c("gene_index", "cov_matrix_list", "beta_list",
#               "individual_log_mean", "individual_log_disp", 
#               "individual_logit_drop", "individual_log_mean_sd"), 
#      file = str_glue("para_data_sim/{est_out}"))
attach(gene_index)
n_cores = NCORE
nGeneTotal = sum(unlist(lapply(gene_index, FUN = length)))

# read in paramters
parameters <- read_file(para_out) %>% #from simulation design step 
  str_split_1(pattern = ":")
parameters <- parameters[parameters != ""]
parameters <- parameters[!duplicated(parameters)]

# global parallel computing settings
registerDoParallel(cores = n_cores)
options(mc.cores = n_cores)

tryCatch({
  for(paras in parameters){
    set.seed(seed = 12334235)
    eval(parse(text=paras)) # VNCELL=FALSE;DIFF=1.1;NCELL=375;NPG=12;BOOT=50
    boot = BOOT
    diff = exp(abs(log(DIFF))) #promise diff >= 1
    n_cas = NPG
    n_ctr = NPG
    nall = n_cas + n_ctr
    neq_n_cell = VNCELL
    n_cell= NCELL
    
    # genes that are modified should be the same across each bootstrap
    # (effectively the same for each scenario)
    random_idx_gene <- sample.int(nGeneTotal)
    flag.list <- lapply(gene_index, FUN = function(genes){
      return(sample.int(length(genes)) <= length(genes)/2)
    })
    
    for(b in 1:boot){ # for FDP calculation, no set.seed inside this loop
      # ------------------------------------------------------------------------
      # estimate 4 parameters for each gene **across individuals** based on 
      # a multivariate-normal distribution estimation for 
      # log_mean, log_disp, logit_drop, log of log_mean_sd. (Already done)
      # simulate the 4 parameters for each individual
      # ------------------------------------------------------------------------
      # at the population level, assume the log_mean, disp, mean_sd follows a MVN
      par.names <- c("mean", "dispersion", "dropout", "mean_sd")
      indi_para <- array(dim = c(nGeneTotal, nall, 4), 
                         dimnames = list(paste("gene", 1:nGeneTotal, sep = "_"), 
                                         paste("ind", 1:nall, sep = "_"), 
                                         par.names))
      
      # covariate at the individual level
      # first apply a normal quantile transformation to RIN, so that 
      # later we can simply simulate RIN from standard normal distribution
      RIN.simu <- rnorm(nall) 
      
      for (ig in 1:nGeneTotal) {# index of gene
        individual_data <- cbind(
          c(individual_log_mean[ig, ]),
          c(individual_log_disp[ig, ]),
          c(individual_logit_drop[ig, ]),
          c(log(individual_log_mean_sd[ig, ]))
        )# 23 rows with the 4 variables
        
        individual_data_mean <- apply(individual_data, 2, mean, na.rm = TRUE) # length of 4
        # cov_matrix <- cov(individual_data)# the covariance between individual_log_mean,
        #individual_log_disp, individual_logit_drop and individual_log_mean_sd
        cov_matrix <- cov_matrix_list[[ig]]
        
        # log_mean_ig <- individual_data[,1]
        # lmi  <- lm(log_mean_ig ~ RIN.qn)
        # beta <- lmi$coefficients # beta: how effective RIN (normal) is in log_mean_ig
        beta <- beta_list[[ig]]
        
        # add some extra variance for the mean parameter
        e1 <- rnorm(nall, mean=0, sd=sqrt(cov_matrix[1,1]))
        log_mean_ig_simu <- beta[1] + beta[2]*RIN.simu + e1 # simulated RIN covariate
        
        for(j in 1:nall){
          individual_data_mean_j    <- individual_data_mean
          individual_data_mean_j[1] <- log_mean_ig_simu[j]
          indi_para[ig, j, ]  <- exp(mvrnorm(1, mu = individual_data_mean_j, 
                                             Sigma = cov_matrix))  
        }# for each gene,
        # for each individual, 
        # the mean is correlated with a individual-level RIN + noise,
        # the dispersion, droupout and sd of mean are population-level averages
        # a MVN of the four parameters (one individual one time, looped for each gene).
      }
      
      # the dropout: expit function, see also the exp part above.
      indi_para[, , 3] <- indi_para[, , 3] / (1 + indi_para[, , 3])
      
      # shuffle genes and samples
      random_idx_sam <- sample.int(nall)
      RIN.simu <- RIN.simu[random_idx_sam]
      indi_para <- indi_para[random_idx_gene, random_idx_sam, ] 
      #random_idx_gene: outside the b in 1:boot loop
      
      # double check the p-value of the covariate RIN
      # pvals.check = rep(NA, nGeneTotal)
      # for (ig in 1:nGeneTotal) {
      #   y_ig = log(indi_para[ig, , 1]) # mean
      #   lm2  = lm(y_ig ~ RIN.simu)
      #   pvals.check[ig] = summary(lm2)$coef[2,4]
      # }
      # summary(pvals.check)# RIN is highly correlated with the mean
      
      # ------------------------------------------------------------------------
      # set up parameters for cases and controls
      # ------------------------------------------------------------------------
      # parameterization of ZINB:
      # mu: mean for the NB part; theta: dispersion for the NB part; pi: dropout
      # mean(NB part) = mu
      # var(NB part) = mu + (mu^2)/theta
      # in the "difference in var" scenario, the mean of NB will not be changed.
      # in the "difference in mean" scenario, the variance of NB will still be 
      # changed
      
      individual_param_case <- indi_para[, 1:n_cas, ] 
      individual_param_ctrl <- indi_para[, (n_cas + 1):nall, ]
      # # check the mean value parameter across genes
      # rd_case = (colSums(individual_param_case[, , 1]))
      # rd_ctrl = (colSums(individual_param_ctrl[, , 1]))
      # t.test(rd_case, rd_ctrl)# no difference, the read-depths are the same
      
      calc_par_var = function(mu, theta, r_v) {
        # the function calc_par_var returns the over-dispersion parameters 
        # for changing the variance of a negative binomial distribution. 
        # to make sure theta is larger than 0, r_v should be > 1. 
        theta2 = theta  * mu / (mu * r_v + (r_v - 1) * theta)
        if(theta2 < 0){ stop("negative theta2.") }
        theta2
      }
      
      # mean_index (attach(gene_index) at the start)
      # have half of the genes increased and the other half decreased for balancing
      flag <- flag.list[["mean_index"]]
      # flag[1:10]
      # table(flag)
      k = 0
      for (i in mean_index) {
        k = k + 1
        if(flag[k]){
          for (j in 1:n_cas) {
            individual_param_case[i, j, 1] = individual_param_case[i, j, 1]*diff# fold change
            # individual_param_case[1:5, 1:2, 1] # [genes, individuals, mean]
            # individual_param_ctrl[1:5, 1:2, 1]
          }
        }else{
          for (j in 1:n_ctr) {
            individual_param_ctrl[i, j, 1] = individual_param_ctrl[i, j, 1]*diff
          }
        }
      }
      # flag[1:5]
      # individual_param_case[1:5, 1:2, 1] # genes, individuals, mean
      # individual_param_ctrl[1:5, 1:2, 1]
      # 
      # rd_case    = colSums(individual_param_case[, , 1])
      # rd_control = colSums(individual_param_ctrl[, , 1])
      # t.test(rd_case, rd_control)# read depths are still balanced
      
      # var_index (attach(gene_index) at the start)
      flag <- flag.list[["var_index"]]
      # flag[1:10]
      # table(flag)
      k = 0
      for (i in var_index) {
        k = k + 1
        if(flag[k]){
          for (j in 1:n_cas) {
            x = individual_param_case[i, j, ]
            individual_param_case[i,j,2] = 
              calc_par_var(mu = x[1], theta = x[2], r_v = diff)
          }# 2 for dispersion
        }else{
          for (j in 1:n_ctr) {
            x = individual_param_ctrl[i, j, ]
            individual_param_ctrl[i,j,2] = 
              calc_par_var(mu = x[1], theta = x[2], r_v = diff)      
          }
        }
      }
      # flag[1:5]
      # individual_param_case[1:5, 1:2, 2] # genes, individuals, dispersion
      # individual_param_ctrl[1:5, 1:2, 2]
      # t.test(colSums(individual_param_case[, , 1]), 
      #        colSums(individual_param_case[, , 1])) # balanced
      
      # mean_var_index (attach(gene_index) at the start)
      flag <- flag.list[["mean_var_index"]]
      k = 0
      for (i in mean_var_index) {
        k = k + 1
        if(flag[k]){
          for (j in 1:n_cas) {
            individual_param_case[i, j, 1] = individual_param_case[i, j, 1]*diff# fold change
            x = individual_param_case[i, j, ]
            individual_param_case[i,j,2] = 
              calc_par_var(mu = x[1], theta = x[2], r_v = diff)
          }
        }else{
          for (j in 1:n_ctr) {
            individual_param_ctrl[i, j, 1] = individual_param_ctrl[i, j, 1]*diff
            x = individual_param_ctrl[i, j, ]
            individual_param_ctrl[i,j,2] = 
              calc_par_var(mu = x[1], theta = x[2], r_v = diff)
          }
        }
      }
      
      
      # ------------------------------------------------------------------------
      # simulate scRNAseq based on zinb parameters of cases and controls
      # ------------------------------------------------------------------------
      # Assume first n_cas individuals are cases, the remaining are controls
      # get the number of cells per individual
      n_cell_each <- rnorm(mean = n_cell,
                           sd = ifelse(neq_n_cell, n_cell_obs[["sd"]][1], 0),
                           n = nall) %>% round()
      n_cell_each[n_cell_each < 100] <- 100 # truncated at 100
      n_cell_tot = sum(n_cell_each)
      n_cell_cumsum = c(0, cumsum(n_cell_each))
      
      
      sim_matrix = matrix(0, nrow = nGeneTotal, ncol = n_cell_tot)
      for(i in 1:nall){#1:20
        # if(i %% 5 ==0){
        #   cat(i, "\n")
        # }
        idx_i = (n_cell_cumsum[i]+1):(n_cell_cumsum[i+1]) 
        # the index that these cells belonging to this individual
        if (i > n_cas) {
          mean_i = individual_param_ctrl[, (i - n_cas), 1] # [gene, indi, para]
          #> str(individual_param_ctrl)
          # num [1:3000, 1:10, 1:4] 2.83 7.11 7.99 2.36 5.35 ...
          # - attr(*, "dimnames")=List of 3
          # ..$ : chr [1:3000] "gene_930" "gene_59" "gene_866" "gene_90" ...
          # ..$ : chr [1:10] "ind_19" "ind_1" "ind_6" "ind_12" ...
          # ..$ : chr [1:4] "mean" "dispersion" "dropout" "mean_sd"
          disp_i = individual_param_ctrl[, (i - n_cas), 2]
          drop_i = individual_param_ctrl[, (i - n_cas), 3]
          individual_mean_sd_i = individual_param_ctrl[, (i - n_cas), 4]
        } else{
          mean_i = individual_param_case[, i, 1]
          disp_i = individual_param_case[, i, 2]
          drop_i = individual_param_case[, i, 3]
          individual_mean_sd_i = individual_param_case[, i, 4]
        }
        
        sim_matrix[,idx_i] = 
          foreach(k = 1:(n_cell_each[i]), .combine=cbind) %dorng% {
            #calcualte the each cell "mean" parameter first, using the mean_i and 
            #individual_mean_sd_i, for each of the 3000 genes.
            #generate a ZINB random number according to the three parameters 
            #(mean, dispersion and droupout)
            #each cell has the mean_i and some variance, and the same 
            #dispersion and dropout
            #looped for 3000 genes.
            individual_mean_k = exp(rnorm(n = nGeneTotal, mean = log(mean_i), 
                                          sd = individual_mean_sd_i))
            sim_vector_cell_k = rep(NA, nGeneTotal)
            for (ig in 1:nGeneTotal) {
              sim_vector_cell_k[ig] = emdbook::rzinbinom(1, individual_mean_k[ig], 
                                                         disp_i[ig], drop_i[ig])
            }
            sim_vector_cell_k
          }
      }
      # dim(sim_matrix)
      # sim_matrix[1:10, 1:5]
      # table(c(sim_matrix) == 0)/(nrow(sim_matrix)*ncol(sim_matrix))
      
      # create meta data
      
      phenotype  = c(rep(1, sum(n_cell_each[1:n_cas])), 
                     rep(0, sum(n_cell_each[(n_cas+1):nall])))
      individual = paste0("ind", c(rep(1:nall, times = n_cell_each)))
      # table(phenotype)
      # sum(n_cell_each[1:10])
      
      # Count info for matrix
      cell_id = paste0("cell", 1:ncol(sim_matrix))
      gene_id = paste0("gene", 1:nrow(sim_matrix))
      rownames(sim_matrix) <- gene_id
      colnames(sim_matrix) <- cell_id
      
      # Cell info for meta
      cell_rd = colSums(sim_matrix) # only this change for each boot
      # CDR     = colSums(sim_matrix > 0) / nrow(sim_matrix)
      meta    = data.frame(cell_id, individual, phenotype, cell_rd, #CDR, 
                           stringsAsFactors=FALSE)
      # dim(meta)
      # meta[1:2,]
      
      meta_ind = meta[, c("individual", "phenotype")]
      meta_ind = unique(meta_ind)
      rownames(meta_ind) = meta_ind$individual
      meta_ind$RIN = RIN.simu
      # dim(meta_ind)
      # meta_ind[1:2,]
      
      mode(sim_matrix) <- "integer"
      dat_list = list(count_matrix = sim_matrix, meta_cell = meta, 
                      meta_ind = meta_ind, gene_index = gene_index, 
                      random_idx_gene = random_idx_gene)
      
      save(dat_list, file = 
             str_glue("n{n_cas}_{n_ctr}_d{diff}_c{n_cell}_v{neq_n_cell}_{dddd}_b{b}.RData"))
      
    }# to the next b
  }# to the next parameter set
}, error=function(err){
  # exception handling. See the part above
  # if error occurs, save current work image and restart outside
  message(err)
  message("Need to be restarted...")
  save.image(file = str_glue("unfinished_gen_sim_data_{dddd}.RData"))
  q(save = "no")
  return(0)
})# try-catach part

stopImplicitCluster()
detach(gene_index)
q(save = "no")