#!/bin/bash
#SBATCH --job-name=t1dca
#SBATCH -o /home/guanwh/zhan7474/mayo/supp_info_anal/t1DCA_IDEAS.out
#SBATCH -e /home/guanwh/zhan7474/mayo/supp_info_anal/t1DCA_IDEAS.err
#SBATCH -p msismall
#SBATCH -c 20
#SBATCH --mem-per-cpu=4000
#SBATCH -t 96:00:00

cd ~
module load R/4.3.0-openblas
module load singularity

# Define the range of numbers to iterate over
start=1
end=20

# Iterate over the range and run R code for each number
for ((i=start; i<=end; i++)); do
    filename1="para_data_sim/n12_12_d1_c375_vFALSE_0811_b${i}.RData"
	filename2="supp_info_anal/null${i}.csv"
	filename3="supp_info_anal/null${i}"
	filename4="supp_info_anal/null${i}.RData"
	
    R --vanilla <<EOF
    setwd("/home/guanwh/zhan7474/mayo/")
    load("$filename1")
    write.csv(dat_list[["count_matrix"]], "$filename2")
    q(save = "no")
EOF
    
	singularity exec --bind /home/guanwh/zhan7474/mayo:/container/input \
                 --bind /home/guanwh/zhan7474/mayo:/container/output \
                 /home/guanwh/zhan7474/tensorflow_2.4.3.sif \
                 /home/guanwh/zhan7474/.local/bin/dca /container/input/$filename2 /container/output/$filename3 --threads 20 --type zinb-conddisp
				 
	R --vanilla <<EOF
    setwd("/home/guanwh/zhan7474/mayo/")
	library(data.table)
    library(ideas)
    library(Matrix)
    library(parallel)
    library(doParallel)
    load("$filename1")
	
	# normalized mean
	meta_cell <- dat_list[["meta_cell"]]
    meta_ind <- dat_list[["meta_ind"]]
    count_matrix <- dat_list[["count_matrix"]]
	
	all_counts = dat_list[["count_matrix"]]
    col_sums = colSums(all_counts)
    med_rd = median(col_sums)
    size_factors = as.numeric(col_sums)/med_rd
	
	mean_mat = fread("$filename3/mean.tsv", header = TRUE, sep = "\t")
    mean_norm_mat = t(t(mean_mat[, 2:ncol(mean_mat)])/size_factors)
    mean_norm_dt = as.data.table(mean_norm_mat)
    mean_norm_dt[, V1:= mean_mat[, 1]]
    setcolorder(mean_norm_dt, "V1")
	
	# DCA-IDEAS
	dca_mean = mean_norm_dt
    dca_disp = fread("$filename3/dispersion.tsv")
    dca_pi   = fread("$filename3/dropout.tsv")
	
	rownames(dca_mean) = dca_mean[["V1"]]
    rownames(dca_disp) = dca_disp[["V1"]]
    rownames(dca_pi)   = dca_pi[["V1"]]
    dca_mean = data.matrix(dca_mean[,-1])
    dca_disp = data.matrix(dca_disp[,-1, with=FALSE])
    dca_pi   = data.matrix(dca_pi[,-1,   with=FALSE])
    dca_par_list = list(dca_mean, dca_disp, dca_pi)
	
	n.cores = 20
    registerDoParallel(cores=n.cores)
    options(mc.cores=n.cores)
    RNGkind("L'Ecuyer-CMRG")
	
	dist1 = ideas_dist(count_input = dca_par_list, 
                       meta_cell = meta_cell, meta_ind = meta_ind, 
                       var_per_cell = "cell_rd", var2test = "phenotype", 
                       var2test_type = "binary",
                       d_metric = "Was", fit_method = "dca_direct")
					   
	dca_ideas_pval_raw <- permanova(dist_array = dist1, meta_ind = meta_ind, 
                                    var2test = "phenotype", var2adjust = "RIN", 
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
                                        var2test = "phenotype", var2adjust = "RIN", 
                                        var2test_type = "binary", n_perm=99999)
      dca_ideas_pval_raw[recal_index] <- dca_ideas_pval_raw_r
    }
    pval_temp <- dca_ideas_pval_raw
    pi0 = 2*sum(pval_temp > 0.5, na.rm = TRUE)/sum(!is.na(pval_temp))
    # really need this pi0, otherwise in some cases, error will occur
    pi0 = ifelse(pi0>1, 1, ifelse(pi0<0.5, 0.5, pi0)) # 0.5<=pi0<=1
    dca_ideas_pval_adj <- qvalue::qvalue(pval_temp, pi0=pi0)[["qvalues"]]
	
	save(dca_ideas_pval_adj, dca_ideas_pval_raw, file = "$filename4")
    q(save = "no")
EOF

# avoid using rm -r
rm /home/guanwh/zhan7474/mayo/$filename3/mean.tsv
rm /home/guanwh/zhan7474/mayo/$filename3/dispersion.tsv
rm /home/guanwh/zhan7474/mayo/$filename3/dropout.tsv
rm /home/guanwh/zhan7474/mayo/$filename3/latent.tsv
rm /home/guanwh/zhan7474/mayo/$filename2
done


