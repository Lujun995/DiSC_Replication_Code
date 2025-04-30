library(tidyverse)
library(Matrix)
library(Seurat)

# workdir = "/home/guanwh/zhan7474/mayo/real_data_anal/AD"

thre_n_cell_per_individual = 30
thre_mito_pct = 5
min_qt_UMI_per_cell = 0.05 # 5% quantile
max_qt_UMI_per_cell = 0.98 # 99% quantile
SPARSITY_THRESHOLD = 0.6 #or NULL. Keep genes with sparsity < SPARSITY_THRESHOLD
options(timeout=600)

# more filtering to reduce computational time (IDEAS will take too much time)
ngene_2_keep = 1000
ncell_2_keep = 250
get_n_cell <- function(read_depth, ncell_2_keep = 300){
  n_cell <- length(read_depth)
  if(ncell_2_keep >= n_cell) 
    return(rep(TRUE, n_cell))
  HDI <- HDInterval::hdi(read_depth, credMass= ncell_2_keep/n_cell, allowSplit = FALSE)
  return(read_depth >= HDI["lower"] & read_depth < HDI["upper"])
}

setwd("/home/guanwh/zhan7474/mayo/res_to_rev/")

cell_type_fp = c(
  Lamp5_Lhx6 = "https://datasets.cellxgene.cziscience.com/3d690bcf-c9d3-4fcf-b7e1-e0e622bbf958.rds",
  Lamp5 = "https://datasets.cellxgene.cziscience.com/abc02bdf-8bb7-42dd-8d5f-eff48fb8225b.rds",
  Pax6 = "https://datasets.cellxgene.cziscience.com/a2ee847d-c5b4-4318-b6ba-0d36ed57d769.rds",
  Sncg = "https://datasets.cellxgene.cziscience.com/ee226a77-6ec1-4a16-b653-8cbacd3876bc.rds",
  Vip = "https://datasets.cellxgene.cziscience.com/08f013eb-107a-4f25-a8b2-4ebb30f6e7ea.rds",
  Chandelier = "https://datasets.cellxgene.cziscience.com/7bb8238f-b5a7-4bbd-9c00-244e2b72e140.rds",
  Pvalb = "https://datasets.cellxgene.cziscience.com/bcdf1ef7-3129-4de8-87a2-066978199d5f.rds",
  Sst = "https://datasets.cellxgene.cziscience.com/b52e6423-1687-49cc-bf31-810c055b6656.rds",
  Sst_Chodl = "https://datasets.cellxgene.cziscience.com/1ac8ed44-27f6-4ade-96e6-e703cb171cfc.rds",
  # L2_3_IT = "https://datasets.cellxgene.cziscience.com/f9e42464-5984-4c77-99ea-9c6d6c130897.h5ad",
  L4_IT = "https://datasets.cellxgene.cziscience.com/94717711-bdeb-4002-94bf-4098b969b061.rds",
  L5_IT = "https://datasets.cellxgene.cziscience.com/2af5392b-7579-43be-bacc-40f67e921fc7.rds",
  L6_IT = "https://datasets.cellxgene.cziscience.com/a37f85f1-4ef4-4d8f-8a65-63e45441b006.rds",
  L6_IT_Car3 = "https://datasets.cellxgene.cziscience.com/5f78b6bc-57ff-469d-ab6a-91787159e88f.rds",
  L5_ET = "https://datasets.cellxgene.cziscience.com/9d53f7bb-dc23-4c05-b2a6-4afa9a6e3be0.rds",
  L5_6_NP = "https://datasets.cellxgene.cziscience.com/97112340-5db0-423a-8077-9386c1b815ff.rds",
  L6b = "https://datasets.cellxgene.cziscience.com/922bd0fa-edd7-44e8-bc7c-59de700a686b.rds",
  L6_CT = "https://datasets.cellxgene.cziscience.com/8406fa2b-a671-4f7f-aaeb-0877e61a2a04.rds",
  Astrocyte = "https://datasets.cellxgene.cziscience.com/eb9ac67f-3041-4cb8-85d4-8616931cc218.rds",
  Oligodendrocyte = "https://datasets.cellxgene.cziscience.com/32b747ce-7afb-4342-89d2-ca5d274164ed.rds",
  OPC = "https://datasets.cellxgene.cziscience.com/8172e13e-dfa0-44fc-aec5-3988ace41f37.rds",
  Endothelial = "https://datasets.cellxgene.cziscience.com/d4f4d813-aa81-4fd9-a12d-b4a60788b604.rds",
  VLMC = "https://datasets.cellxgene.cziscience.com/0fd1a3ed-51fd-44ba-808c-9b94f3ff88da.rds",
  Microglia_PVM = "https://datasets.cellxgene.cziscience.com/32b32549-9c68-453d-91c6-74acad9df928.rds"
)

cell_types <- names(cell_type_fp)

# variables should be included in the meta
var_interested_cell <- c("cell", "Number of UMIs", "Fraction mitochrondrial UMIs", "assay")
var_interested_indi <- c("donor_id", "disease", 
                         "Age at death", "sex", "self_reported_ethnicity")
var_interested <- c(var_interested_cell, var_interested_indi)

for(ctp in cell_types){ #ctp = cell_types[1]
  cat("#########################################################\n")
  cat("Cleaning data for cell type:", ctp, "\n")
  cat("#########################################################\n")
  
  DATA <- readRDS(url(cell_type_fp[ctp]))
  ct_mat <- LayerData(DATA, assay = "RNA", layer = "counts") # looks good
  meta <- DATA[[]]
  meta$cell <- rownames(meta)
  
  #######################################
  # check meta and ct_mat
  #######################################
  meta = droplevels(meta)
  print(str(meta))
  for(vv in var_interested)
    stopifnot(vv %in% names(meta)) #an error occurs when any of var_interested is not in meta
  stopifnot(inherits(meta[["cell"]], "character"))
  stopifnot(inherits(meta[["assay"]], "factor"))
  stopifnot(inherits(meta[["Number of UMIs"]], "numeric"))
  stopifnot(inherits(meta[["Fraction mitochrondrial UMIs"]], "numeric"))
  stopifnot(inherits(meta[["donor_id"]], "factor")) # this should be converted to character
  stopifnot(inherits(meta[["disease"]], "factor"))
  stopifnot(inherits(meta[["Age at death"]], "factor")) # this should be converted to numeric
  stopifnot(inherits(meta[["sex"]], "factor"))
  stopifnot(inherits(meta[["self_reported_ethnicity"]], "factor"))
  
  stopifnot(inherits(ct_mat, "dgCMatrix"))
  stopifnot(sum(colnames(ct_mat) != meta[["cell"]]) == 0)
  read_depth <- colSums(ct_mat)
  if(sum(read_depth != meta[["Number of UMIs"]]) != 0){
    cat("colSums(ct_mat) != Number of UMIs\n")
    if(sum(read_depth > meta[["Number of UMIs"]]) == 0){
      cat("colSums(ct_mat) <= Number of UMIs\n") 
      # this may be because some sequences have been discarded in QC steps
      
      index <- 
        (read_depth <= quantile(read_depth, max_qt_UMI_per_cell)) &
        (read_depth >= quantile(read_depth, min_qt_UMI_per_cell))
      diff = read_depth - meta[["Number of UMIs"]]
      cat("summary of the difference:\n")
      print(summary(diff[index]))
      # We will remove insufficiently sequenced cells and overly deeply sequenced cells later
      # As long as the Number of UMIs in the remaining cells is close to read_depth
      # it is good to replace Number of UMIs with colSums(ct_mat), or read_depth
      cat("summary of the difference in percentage:\n")
      frac = diff / meta[["Number of UMIs"]]
      print(summary(frac[index]))
      
      if(max(abs(frac[index])) <= 0.03){
        cat("Number of UMIs * 0.97 <= colSums(ct_mat) <= Number of UMIs\n")
        cat("Replacing Number of UMIs with colSums(ct_mat)\n")
        # when the difference is small, "Fraction mitochrondrial UMIs" are still considerably accurate
        # Fraction mitochrondrial UMIs = 5% -> 5.15%
        meta[["Number of UMIs"]] <- read_depth
      } else {
        stop("for some colSums(ct_mat) < 0.97 * Number of UMIs! It may be inaccurate to replace Number of UMIs")
      }
    } else {
      stop("for some colSums(ct_mat) > Number of UMIs!!!!!")
    }
  }
  stopifnot(sum(colSums(ct_mat) != meta[["Number of UMIs"]]) == 0)
  
  ############################################
  # clean the dataset
  ############################################
  meta_cleaned <- meta
  ##outcome of interest: disease (dementia vs normal)
  table(meta_cleaned[["Cognitive status"]], meta_cleaned[["disease"]]) 
  ##one donor can have multiple specimens
  # table(meta_cleaned[["donor_id"]], meta_cleaned[["Specimen ID"]]) 
  ##check if each individual will have a unique value for outcome of interest and covariates
  cat("# of individuals:", length(unique(meta_cleaned[["donor_id"]])), "\n")
  stopifnot(length(unique(meta_cleaned[["donor_id"]])) == # unique donors
              nrow(distinct(meta_cleaned, pick(all_of(var_interested_indi)))))
  meta_cleaned %>% 
    distinct(donor_id, disease, 
             sex, `Age at death`, self_reported_ethnicity) %>%
    tableone::CreateTableOne(data =., vars = setdiff(names(.), c("disease", "donor_id")),
                             strata = c("disease"),
                             includeNA = TRUE, test = FALSE, addOverall = TRUE)%>%
    print()# balanced
  
  ######################################################
  # filtering and datatype conversion
  ######################################################
  # datatype conversion
  meta_cleaned <- meta_cleaned %>%
    select(all_of(var_interested)) %>%
    rename(individual = "donor_id", cell_id = "cell") %>%
    mutate(individual = as.character(individual)) %>%
    mutate(age_lvl = case_match(`Age at death`, 
                                "Less than 65 years old" ~ 1, 
                                "65 to 77 years old" ~ 2,
                                "78 to 89 years old" ~ 3,
                                "90+ years old" ~ 4))
  
  # self_reported_ethnicity: subjects are dominantly European descendants and
  # the sample sizes for other ethnicity are too small (1 or 3 people)
  meta_cleaned <- meta_cleaned %>% 
    filter(self_reported_ethnicity == "European") %>%
    filter(assay == "10x 3' v3") # the number squenced by "10x 3' v3" was 9x as many as "10x multiome"
  
  
  # mitochrondrial RNA fraction
  print(summary(meta_cleaned[["Fraction mitochrondrial UMIs"]]))
  if(sum(meta_cleaned[["Fraction mitochrondrial UMIs"]] > 1) == 0){
    cat("Fraction mitochrondrial UMIs are decimals. Convert them to percentage\n")
    meta_cleaned[["Fraction mitochrondrial UMIs"]] <- 
      meta_cleaned[["Fraction mitochrondrial UMIs"]]*100
  }
  cat(str_glue("number of cells with Fraction mitochrondrial UMIs greater than {thre_mito_pct}: {sum(meta_cleaned[['Fraction mitochrondrial UMIs']]>thre_mito_pct)}\n\n"))
  meta_cleaned <- meta_cleaned %>% 
    filter(`Fraction mitochrondrial UMIs` <= thre_mito_pct)
  
  # number of UMIs
  max_count_UMI_per_cell <- 
    quantile(meta_cleaned[["Number of UMIs"]], probs = max_qt_UMI_per_cell)
  min_count_UMI_per_cell <- 
    quantile(meta_cleaned[["Number of UMIs"]], probs = min_qt_UMI_per_cell)
  cat("filtering threshold of Number of UMIs: min:", min_count_UMI_per_cell,
      ", max:", max_count_UMI_per_cell, "\n")
  meta_cleaned <- meta_cleaned %>%
    filter(`Number of UMIs` <= max_count_UMI_per_cell &
             `Number of UMIs` >= min_count_UMI_per_cell)
  
  # filtering insufficient cell per individual
  # this step was less explicit in COVID-PBMC pipeline, where filtering starts 
  # with sample_ids -> columns in mat_filtered -> rows in meta_cell
  ncell_tbl <- meta_cleaned %>% 
    summarise(ncell = n(), .by = individual)%>%
    filter(ncell < thre_n_cell_per_individual)
  cat(str_glue("individuals with ncell < {thre_n_cell_per_individual}:"), 
      ncell_tbl[['individual']], "\n")
  meta_cleaned <- meta_cleaned %>%
    filter(!(individual %in% ncell_tbl[["individual"]]))
  
  ######################################################
  # creating the required dataset for DE analyses
  ######################################################
  var2test = "disease"
  #var2adjust = c("age_lvl", "sex", "self_reported_ethnicity") 
  var2adjust = c("age_lvl", "sex") 
  var2test_type = "binary"
  
  meta_ind <- meta_cleaned %>% 
    distinct(pick(all_of(c("individual", var2test, var2adjust)))) %>%
    droplevels()
  if(nrow(meta_ind) < 6){
    cat(str_glue("Sample size too small for {ctp} (<6 individuals). {ctp} will be excluded!\n\n"))
    next
  }
  if(sum(meta_ind[["disease"]] == "normal") < 2){
    cat(str_glue("For {ctp}, we only have 0 or 1 healthy control. {ctp} will be excluded!\n\n"))
    next
  }
  if(sum(meta_ind[["disease"]] == "dementia") < 2){
    cat(str_glue("For {ctp}, we only have 0 or 1 dementia case. {ctp} will be excluded!\n\n"))
    next
  }
  
  count_matrix <- ct_mat[, colnames(ct_mat) %in% meta_cleaned[["cell_id"]]]
  ## count_matrix, genes in rows and cells in columns
  # no need to convert to dense matrix here. Do it after sparsity filtering
  ngene = nrow(count_matrix)
  ncell = ncol(count_matrix)
  read_depth <- colSums(count_matrix)
  if(ngene < ncell) 
    message("the number of genes is smaller than the number of cells, please check if genes are in rows")
  cat("ngene: ", ngene, "ncell: ", ncell, "read_depth[1:10]: ", read_depth[1:10], "\n")
  #### Sparsity filtering
  if(!is.null(SPARSITY_THRESHOLD)){
    cat("sparsity filtering with threshold ", SPARSITY_THRESHOLD, " :\n")
    gene2keep = which(rowSums(count_matrix == 0) < SPARSITY_THRESHOLD*ncell) # all zero rows (100% sparsity) will be removed here
    count_matrix = count_matrix[gene2keep,]
    read_depth = colSums(count_matrix)
    ngene = nrow(count_matrix)
    ncell = ncol(count_matrix)
    read_depth <- colSums(count_matrix)
    cat("ngene: ", ngene, "ncell: ", ncell, "read_depth[1:10]: ", read_depth[1:10], "\n")
  }
  count_matrix <- as.matrix(count_matrix)
  
  meta_cell <- meta_cleaned
  stopifnot(sum(colnames(count_matrix) != meta_cell[["cell_id"]]) == 0)
  meta_cell$read_depth <- read_depth
  meta_cell <- droplevels(meta_cell)
  
  ## prepossessing
  meta_cell <- as.data.frame(meta_cell) 
  meta_cell$individual <- as.character(meta_cell$individual)
  meta_cell$cell_id <- as.character(meta_cell$cell_id)
  rownames(meta_cell) <- NULL
  meta_cell <- droplevels(meta_cell)
  str(meta_cell)
  meta_ind <- as.data.frame(meta_ind)
  meta_ind$age_lvl <- scale(meta_ind$age_lvl)
  meta_ind$sex <- as.factor(meta_ind$sex)
  # meta_ind$self_reported_ethnicity <- as.factor(meta_ind$self_reported_ethnicity)
  meta_ind$individual <- as.character(meta_ind$individual)
  rownames(meta_ind) <- NULL
  meta_ind <- droplevels(meta_ind)
  str(meta_ind)

  # more filtering to reduce the computational time 
  # We have to reduce the computational time otherwise IDEAS will take a VERY long time to complete the DE analsyses
  if(!is.null(ngene_2_keep)){
    ngene <- nrow(count_matrix)
    if(ngene > ngene_2_keep) {
      cat("More sparsity filtering...\n")
      ncell <- ncol(count_matrix)
      sparsity <- rowSums(count_matrix == 0)/ncell
      gene2keep <- (rank(sparsity) <= ngene_2_keep)
      count_matrix <- count_matrix[gene2keep,]
      meta_cell$read_depth <- colSums(count_matrix)

      ngene <- nrow(count_matrix)
      ncell <- ncol(count_matrix)
      cat("ngene: ", ngene, "ncell: ", ncell, "read_depth[1:10]: ", meta_cell$read_depth[1:10], "\n")
    }
  }
  if(!is.null(ncell_2_keep)){
    stopifnot(meta_cell$cell_id == colnames(count_matrix))
    cell_include <- meta_cell %>%
      group_by(individual) %>%
      mutate(cell_include = get_n_cell(read_depth, ncell_2_keep)) %>%
      pull(cell_include)
    count_matrix <- count_matrix[, cell_include]
    meta_cell <- meta_cell[cell_include, ]
    meta_cell %>% group_by(individual) %>% 
      summarise(min= min(read_depth), max = max(read_depth), n_cell =n()) %>% 
      print(n=100)
    ngene <- nrow(count_matrix)
    ncell <- ncol(count_matrix)
    cat("ngene: ", ngene, "ncell: ", ncell, "read_depth[1:10]: ", meta_cell$read_depth[1:10], "\n")
  }
  
  var2adjust = NULL
  save(list = c("count_matrix", "meta_cell", "meta_ind",
                "var2test", "var2adjust", "var2test_type"),
       file = str_glue("/home/guanwh/zhan7474/mayo/res_to_rev/AD2/{ctp}.RData"))
  # this dataset should have:
  # count_matrix: genes in rows, cells in cols, rowname: gene_name, colname: cell_id
  # meta_cell: with variables called "individual" (string), "cell_id" (string) and "read_depth" (numeric)
  # meta_ind: with variables called "individual", the outcome and all the covariates (numeric or factor)
  # var2test, var2adjust, var2test_type: strings
  
  rm(list = setdiff(ls(), c("ngene_2_keep", "ncell_2_keep", "get_n_cell",
                            "thre_n_cell_per_individual",
                            "thre_mito_pct", "min_qt_UMI_per_cell",
                            "max_qt_UMI_per_cell", "SPARSITY_THRESHOLD",
                            "cell_types", "cell_type_fp", "var_interested",
                            "var_interested_cell", "var_interested_indi",
                            "ctp")))
  gc()
  gc()
}




