zp_prediction <- function (mat, norm_factor, span) 
{
  if (class(mat)[1] != "dgCMatrix") {
    mat <- Matrix::Matrix(mat, sparse = TRUE)
  }
  prop_zero <- 1 - tabulate(mat@i + 1, nbins = nrow(mat))/ncol(mat)
  names(prop_zero) <- rownames(mat)
  prop_zero <- prop_zero[prop_zero != 1]
  gl = names(prop_zero)
  norm_gene_all <- try(apply(mat[gl, ], 1, function(gene) {
    mean(log1p(gene/norm_factor * 10000), na.rm = T)
  }))
  if (class(norm_gene_all) == "try-error") {
    tmp_seurat <- Seurat::CreateSeuratObject(counts = mat[gl, 
      ], min.cells = 0, min.features = 0)
    tmp_seurat <- Seurat::NormalizeData(tmp_seurat)
    norm_gene_all <- Matrix::rowMeans(tmp_seurat@assays$RNA@data)
  }
  zp_plot = data.frame(zp = prop_zero, log_mu = log(norm_gene_all))
  loessMod05 <- stats::loess(zp ~ log_mu, data = zp_plot, 
    span = span)
  smoothed05 <- stats::predict(loessMod05)
  names(smoothed05) <- gl
  plot(zp ~ log_mu, data = zp_plot, ylab = "Proportion of zeros", 
    xlab = "log(avg.norm.expression)", main = paste0("span=", 
      span), cex = 0.3)
  graphics::lines(smoothed05[order(zp_plot$log_mu)], x = zp_plot$log_mu[order(zp_plot$log_mu)], 
    col = "red")
  return(smoothed05)
}