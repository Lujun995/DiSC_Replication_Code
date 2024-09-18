####################################################
# setup, data input and exploratory analysis
####################################################
library(tidyverse)
library(Matrix)
library(Seurat)
library(SeuratDisk)
options(timeout=180)

workdir = "/home/guanwh/zhan7474/mayo/real_data_anal"
fp = "https://covid19.cog.sanger.ac.uk/submissions/release1/haniffa21.processed.h5ad"

DELETE = FALSE # delete the downloaded data files
intermediate_results = FALSE # save intermediate results
thre_n_cell_per_individual = 30
thre_mito_pct = 5
thre_count_per_cell = 3800 # 5% quantile
max_count_per_cell  = 37000 # 99% quantile
SPARSITY_THRESHOLD = 0.8 #or NULL. Keep genes with sparsity < SPARSITY_THRESHOLD

setwd(workdir)

if(!file.exists("haniffa21.processed.h5ad"))
  curl::curl_download(fp, "haniffa21.processed.h5ad", quiet = FALSE)
if(!file.exists("haniffa21.processed.h5seurat"))
  Convert("haniffa21.processed.h5ad", dest = "h5seurat", verbose = TRUE)
covid_pbmc <- LoadH5Seurat("haniffa21.processed.h5seurat")
covid_pbmc
meta <- covid_pbmc[[]]
dim(meta)

#fp = "http://cells.ucsc.edu/covid19-pbmc"
#meta2 = read_tsv(str_glue("{fp}/meta.tsv"))
#sum(abs(meta$pct_counts_mt - meta2$pct_counts_mt)) #0.02773324, two datasets are the same

meta[["cellId"]] <- rownames(meta)
meta$sample_id <- as.character(meta$sample_id)
covid_pbmc[["cellId"]] <- rownames(meta)
#sum(abs(meta$nCount_raw - colSums(covid_pbmc[["raw"]]$counts))) #0
covars_of_interest <- 
  c("Site", "Sex", "Age_interval", "Smoker",
    "Status", "Worst_Clinical_Status", "Swab_result", "Outcome", 
    "Collection_Day","Days_from_onset",
    "Status_on_day_collection", "Status_on_day_collection_summary")
unique(meta$sample_id) %>% length() # 143
persons_covar <- meta %>% 
  select("sample_id", all_of(covars_of_interest)) %>%
  distinct() 
dim(persons_covar) # 143 rows, individual-level covariates have only one value for each subject across all cells

persons_cell_info <- meta %>% 
  group_by(sample_id) %>%
  summarise(n_cells = n(), 
            n_cells_mito_pct_l5 = sum(pct_counts_mt<5),
            mean_genes_per_cell = mean(n_genes),
            sd_genes_per_cell = sd(n_genes),
            median_seq_depth = median(nCount_raw),
            min_seq_depth = min(nCount_raw),
            IQR_seq_depth = IQR(nCount_raw),
            .groups = "drop")

persons_cellgrp_info <- meta %>%
  group_by(sample_id, initial_clustering) %>%
  summarise(n_cell = n(), .groups = "drop") %>%
  pivot_wider(id_cols = sample_id, names_from = initial_clustering,
              values_from = n_cell, names_prefix = "n_")

persons <- left_join(
  left_join(persons_covar, persons_cell_info, by = "sample_id"),
  persons_cellgrp_info, by = "sample_id") %>%
  rename("Status_on_collection" = "Status_on_day_collection_summary")

persons %>%
  mutate(across(starts_with("n_"), 
                \(x)
                cut(x = x,
                    breaks = c(0, 49.9, 199.9, 999.9, 4999.9, Inf),
                    labels = c("< 50", "[50, 200)", "[200, 1000)",
                               "[1000, 5000)", ">= 5000"),
                    ordered_result = TRUE, include.lowest = TRUE) )) %>%
  select("Site", "Status_on_collection",
         "Sex", "Age_interval", "Smoker",
         starts_with("n_"), 
         "mean_genes_per_cell", "min_seq_depth", "median_seq_depth") %>%
  tableone::CreateTableOne(data =.,
                           vars = setdiff(names(.), "Site"),
                           # strata = c("Site", "Status_on_collection"),
                           strata = c("Site"),
                           includeNA = TRUE, test = FALSE, addOverall = TRUE)

####################################################
# Filtering start
####################################################
# person_filtered: filtered individual-level covariate table
# filtering based on Status_on_collection of interest and site
person_filtered <- persons %>%
  filter(Status_on_collection %in% c("Healthy", "Mild", "Moderate", 
                                     "Severe", "Critical")) %>%
  mutate(Status_on_collection = 
           ordered(Status_on_collection, 
                   levels = c("Healthy", "Mild", "Moderate", 
                              "Severe", "Critical"))) %>%
  mutate(Status_on_collection = 
           fct_collapse(Status_on_collection,
                        "Severe/Critical" = c("Severe", "Critical"))) %>%
  filter(Site == "Cambridge") %>% 
  droplevels()
# More balanced
person_filtered %>% 
  select("Site", "Status_on_collection",
         "Sex", "Age_interval", "Smoker") %>%
  tableone::CreateTableOne(data =.,
                           includeNA = TRUE, test = FALSE, addOverall = TRUE)

# cell_filtered: filtered cell-level covariate table 
#Excluding insufficiently sequenced cells, high pct_counts_mt cells, 
#and overly deeply sequenced cells (in the 99th percentile). 
meta %>% filter(sample_id %in% person_filtered$sample_id) %>%
  summarise(nCount_q5 = quantile(nCount_raw, prob = 0.05),
            nCount_q99 = quantile(nCount_raw, prob = 0.99))
cell_filtered <- meta %>% 
  select("cellId", "sample_id", "n_genes", "n_genes_by_counts",
         "nCount_raw", "total_counts_mt", "pct_counts_mt", 
         "full_clustering", "initial_clustering") %>%
  filter(sample_id %in% person_filtered$sample_id) %>%
  filter(pct_counts_mt < thre_mito_pct) %>%
  filter(nCount_raw >= thre_count_per_cell & nCount_raw <= max_count_per_cell)

####################################################
# Splitting, sparsity filtering, and Prepare the data for DE methods
####################################################
#This is to count how many cells are available for analysis in each individual 
#for each cell type and to exclude individuals with an insufficient number of 
#cells (set to NA and excluded).
persons_cellgrp_info_filtered <- cell_filtered %>%
  group_by(sample_id, initial_clustering) %>%
  summarise(n_cell = n(), .groups = "drop") %>%
  pivot_wider(id_cols = sample_id, names_from = initial_clustering,
              values_from = n_cell, names_prefix = "n_") %>%
  mutate(across(starts_with("n_"), 
                \(x) if_else(x < thre_n_cell_per_individual, NA, x))) %>%
  pivot_longer(cols = starts_with("n_"), 
               values_to = "n_cell", values_drop_na = TRUE,
               names_to = "cell_type", names_prefix = "^n_")
cell_types <- as.character(unique(meta$initial_clustering))
# Prepare the data for DE methods (DiSC, DESeq2, IDEAS ...)
meta_cell_all <- cell_filtered %>% 
  dplyr::rename(`cell_id` = "cellId", # IDEAS requires the variable names to be exactly "cell_id" and "individual"
                `individual` = "sample_id") %>%
  mutate(individual = as.character(individual),
         cell_id = as.character(cell_id)) %>%		
  droplevels()
meta_ind_all <- person_filtered %>% 
  select(sample_id, Site, Sex, Age_interval, 
         Status, Status_on_collection) %>%
  mutate(Age = case_match(Age_interval, "(20, 29]"~25,
                          "(30, 39]"~35, "(40, 49]"~45,
                          "(50, 59]"~55, "(60, 69]"~65,
                          "(70, 79]"~75, "(80, 89]"~85),
         Status_on_collection = 
           case_match(Status_on_collection, "Healthy" ~ 0,
                      "Mild" ~ 1, "Moderate" ~ 2, 
                      "Severe/Critical" ~ 3)) %>%
  select(- Age_interval) %>%
  dplyr::rename(`individual` = "sample_id") %>%
  mutate(individual = as.character(individual)) %>%
  droplevels()
var2test <- "Status_on_collection"
var2adjust <- c("Age", "Sex")
var2test_type <- "continuous"

for(ct in cell_types){
  cat("Splitting cell type:", ct, "\n")
  # sample_ids: for each cell type, the samples/individuals who have a status 
  # and site of interest and a sufficient number of cells 
  sample_ids <- persons_cellgrp_info_filtered %>% 
    filter(cell_type == ct) %>%
    pull(sample_id)
  
  # if the number of samples/individuals for this cell type is too small 
  # (in total or in either treatment arm), the cell type will not be analyzed
  if(length(sample_ids) < 10){
    cat(str_glue("Sample size too small for {ct} (<10 individuals)."))
    next
  }
  status <- person_filtered %>% 
    filter(sample_id %in% sample_ids) %>%
    pull(Status)
  if(sum(status == "Healthy") < 2){
    cat(str_glue("For {ct}, we only have 0 or 1 Healthy sample. {ct} will be excluded!"))
    next
  }
  if(sum(status == "Covid") < 2){
    cat(str_glue("For {ct}, we only have 0 or 1 COVID sample. {ct} will be excluded!"))
    next
  }
  
  # split the count data into different cell types
  # Seurat v5 assays store data in layers. These layers can store raw, 
  # un-normalized counts (layer='counts'), normalized data (layer='data'), or 
  # z-scored/variance-stabilized data (layer='scale.data').
  # The layer "data" == layer "counts", but assay "RNA" comprises normalized 
  # counts (having negative values) while assay "raw" comprises raw counts
  mat_filtered <- subset(covid_pbmc, 
                         subset = initial_clustering == ct & # this cell type
                           sample_id %in% sample_ids & # individuals with enough cells and with status and sites of interest
                           pct_counts_mt < thre_mito_pct & # low mito_pct
                           nCount_raw >= thre_count_per_cell & # sufficiently sequenced
                           nCount_raw <= max_count_per_cell) %>% # not overly deep-sequenced
    LayerData(assay = "raw", layer = "counts")
  mat_filtered <- mat_filtered[rowSums(mat_filtered) != 0, ]
  mat_filtered[1:50, 1:2]
  
  # mat_filtered2 <- subset(covid_pbmc, 
  # subset = initial_clustering == ct & # this cell type
  # sample_id %in% sample_ids & # individuals with enough cells and with status and sites of interest
  # pct_counts_mt < thre_mito_pct & # low mito_pct
  # total_counts >= thre_count_per_cell) %>% # sufficiently sequenced
  # LayerData(assay = "raw", layer = "data")
  # mat_filtered2 <- mat_filtered2[rowSums(mat_filtered2) != 0, ]
  # mat_filtered2[1:50, 1:2]
  # sum(abs(mat_filtered - mat_filtered2)) # 0
  
  if(intermediate_results)
    saveRDS(mat_filtered, str_glue("{workdir}/COVID_PBMC/{ct}.rds"), 
            compress = "xz")
  
  # Prepare the data for DE methods (DiSC, DESeq2, IDEAS ...)
  ## count_matrix, genes in rows and cells in columns
  count_matrix = mat_filtered # no need to convert to dense matrix here. Do it after sparsity filtering
  ngene = nrow(count_matrix)
  ncell = ncol(count_matrix)
  read_depth <- colSums(count_matrix)
  if(ngene < ncell) 
    message("the number of genes is smaller than the number of cells, please check if genes are in rows")
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
  count_matrix <- as.matrix(count_matrix)
  ## meta data
  meta_cell <- meta_cell_all
  if(sum(!(colnames(count_matrix) %in% meta_cell$cell_id)) >0)
    error("meta_cell and count_matrix dont match") else{
      meta_cell <- meta_cell %>% filter(cell_id %in% colnames(count_matrix))
      if(sum(meta_cell$cell_id != colnames(count_matrix)) >0 )
        error("meta_cell and count_matrix dont match2")
      meta_cell$read_depth <- read_depth
    }
  meta_cell <- droplevels(meta_cell)
  ## meta_ind
  meta_ind <- meta_ind_all %>% 
    filter(individual %in% meta_cell$individual)
  meta_ind <- droplevels(meta_ind)
  nindi = nrow(meta_ind)
  ## prepossessing
  meta_cell <- as.data.frame(meta_cell) 
  meta_cell$individual <- as.character(meta_cell$individual)
  meta_cell$cell_id <- as.character(meta_cell$cell_id)
  meta_cell <- droplevels(meta_cell)
  str(meta_cell)
  meta_ind <- as.data.frame(meta_ind)
  meta_ind$Age <- scale(meta_ind$Age)
  meta_ind$Sex <- as.factor(meta_ind$Sex)
  meta_ind$individual <- as.character(meta_ind$individual)
  meta_ind <- droplevels(meta_ind)
  str(meta_ind)
  
  save(list = c("count_matrix", "meta_cell", "meta_ind",
                "var2test", "var2adjust", "var2test_type"),
       file = str_glue("{workdir}/COVID_PBMC/{ct}.RData"),
       compress = "xz")
  
  gc() 

}

if(intermediate_results)
  save(list = c("person_filtered", "cell_filtered"), 
       file = str_glue("{workdir}/COVID_PBMC/meta.RData"))

if (DELETE) {
  #Delete file if it exists
  file.remove("haniffa21.processed.h5ad")
  file.remove("haniffa21.processed.h5seurat")
}