perm.fdr.adj <- function (F0, Fp) {
  ord <- order(F0, decreasing = T)
  F0 <- F0[ord]
  perm.no <- ncol(Fp)
  Fp <- as.vector(Fp)
  Fp <- Fp[!is.na(Fp)]
  Fp <- sort(c(Fp, F0), decreasing = F)
  n <- length(Fp)
  m <- length(F0)
  FPN <- (n + 1) - match(F0, Fp) - 1:m
  p.adj.fdr <- FPN / perm.no / (1:m)
  
  # Impute 0s - pseudo-ct 0.5
  p.adj.fdr[p.adj.fdr == 0] <- 0.5 / perm.no / which(p.adj.fdr == 0)
  
  p.adj.fdr <- pmin(1, rev(cummin(rev(p.adj.fdr))))[order(ord)]
  
  return(p.adj.fdr)
}

perm.fwer.adj <- function (F0, Fp) {
  ord <- order(F0, decreasing = T)
  m <- length(F0)
  F0 <- F0[ord]
  Fp <- Fp[ord, , drop=FALSE]
  col.max <- Fp[m, ]
  p.adj.fwer <- sapply(m:1, function(i) {
    x <- F0[i]
    y <- Fp[i, ]
    col.max <<- ifelse(y > col.max, y, col.max)
    # Impute 0s
    mean(c(1, col.max >= x))
  })
  p.adj.fwer <- rev(p.adj.fwer)
  
  p.adj.fwer <- pmin(1, rev(cummin(rev(p.adj.fwer))))[order(ord)]
  return(p.adj.fwer)
}

na.pad <- function (vec, ind) {
  vec0 <- numeric(length(ind))
  vec0[!ind] <- vec
  vec0[ind] <- NA
  return(vec0)
}

basic.func <-
  function (meta.dat,
            feature.list,
            grp.name,
            adj.name = NULL,
            stats.combine.func = mean,
            perm.no = 999,
            strata = NULL,
            is.fwer = FALSE,
            verbose = TRUE)
  {
    this.call <- match.call()
    
    
    feature.no <- length(feature.list)
    
    sample.no <- ncol(feature.list[[1]])
    otu.no <- nrow(feature.list[[1]])
    row.names <- rownames(feature.list[[1]])
    
    if (verbose)
      cat("The data has ",
          sample.no,
          " samples and ",
          otu.no,
          " genes will be tested!\n")
    
    
    if (!is.null(strata)) {
      strata <- factor(strata)
    }
    if (is.null(adj.name)) {
      M0 <- model.matrix( ~ 1, meta.dat)
    }
    else {
      data0 <- meta.dat[, c(adj.name), drop = FALSE]
      if (sum(is.na(data0)) != 0) {
        stop("Please remove or impute NAs in the variables to be adjusted!\n")
      }
      M0 <- model.matrix( ~ ., data0)
    }
    data1 <- meta.dat[, c(grp.name), drop = FALSE]
    if (sum(is.na(data1)) != 0) {
      stop("Please remove or impute NAs in the variable of interest!\n")
    }
    M1 <- model.matrix( ~ ., data1)[,-1, drop = FALSE]
    M01 <- cbind(M0, M1)
    qrX0 <- qr(M0, tol = 1e-07)
    Q0 <- qr.Q(qrX0)
    Q0 <- Q0[, 1:qrX0$rank, drop = FALSE]
    H0 <- (Q0 %*% t(Q0))
    qrX1 <- qr(M1, tol = 1e-07)
    Q1 <- qr.Q(qrX1)
    Q1 <- Q1[, 1:qrX1$rank, drop = FALSE]
    qrX01 <- qr(M01, tol = 1e-07)
    Q01 <- qr.Q(qrX01)
    Q01 <- Q01[, 1:qrX01$rank, drop = FALSE]
    R0 <- as.matrix(resid(lm(Q1 ~ Q0 - 1)))
    pX0 <- ncol(Q0)
    pX1 <- ncol(Q1)
    pX01 <- ncol(Q01)
    df.model <- pX01 - pX0
    df.residual <- sample.no - pX01
    

    if (verbose)
      cat("Permutation testing ...\n")

    Y <- matrix(NA, sample.no, feature.no * otu.no)
    
    for (j in 1:length(feature.list)) {
      Y[, feature.no * (0:(otu.no - 1)) + j] <- t(feature.list[[j]])
    }
    
    Y <- t(Y)
    TSS <- rowSums(Y ^ 2)
    MSS01 <- rowSums((Y %*% Q01) ^ 2)
    MSS0 <- rowSums((Y %*% Q0) ^ 2)
    MSS <- (MSS01 - MSS0)
    RSS <- (TSS - MSS01)
    getPermuteMatrix <- getFromNamespace("getPermuteMatrix",
                                         "vegan")
    perm.ind <-
      getPermuteMatrix(perm.no, sample.no, strata = strata)
    perm.no <- nrow(perm.ind)
    MRSSp <- sapply(1:perm.no, function(ii) {
      if (verbose) {
        if (ii %% 10 == 0)
          cat(".")
      }
      Rp <- R0[perm.ind[ii,], , drop = FALSE]
      Rp <- Rp - H0 %*% Rp
      qrRp <- qr(Rp, tol = 1e-07)
      Q1p <- qr.Q(qrRp)
      Q1p <- Q1p[, 1:qrRp$rank, drop = FALSE]
      MSS01p <- MSS0 + rowSums((Y %*% Q1p) ^ 2)
      MSSp <- (MSS01p - MSS0)
      RSSp <- (TSS - MSS01p)
      c(MSSp, RSSp)
    })

    unit <- feature.no * otu.no
    MSSp <- MRSSp[1:unit,]
    RSSp <- MRSSp[(unit + 1):(2 * unit),]
    RSS.m <- t(array(RSS, c(feature.no, otu.no)))

    F0.m <- t(array((MSS / df.model) / (RSS / df.residual), c(feature.no,
                                                            otu.no)))

    R2.m <- t(array(MSS / TSS, c(feature.no, otu.no)))

    F0 <- (MSS / df.model) / (RSS / df.residual)
    Fp <- (MSSp / df.model) / (RSSp / df.residual)
    F0 <- array(F0, c(feature.no, otu.no))
    Fp <- array(Fp, c(feature.no, otu.no, perm.no))
    F0 <- apply(F0, 2, stats.combine.func)
    Fp <- apply(Fp, c(2, 3), stats.combine.func)

    if (verbose)
      cat("\n")
    if (mean(is.na(F0)) >= 0.1) {
      warning("More than 10% observed F stats have NA! Please check! \n")
    }
    if (mean(is.na(Fp)) >= 0.1) {
      warning("More than 10% permuted F stats have NA! Please check! \n")
    }
    na.ind <- is.na(F0)
    F0 <- F0[!na.ind]
    Fp <- Fp[!na.ind,]
    which.nan.ind <- which(!na.ind)
    p.raw <- rowMeans(cbind(Fp, F0) >= F0)
    
    p.adj.fdr <- perm.fdr.adj(F0, Fp)
    p.raw <- na.pad(p.raw, na.ind)
    p.adj.fdr <- na.pad(p.adj.fdr, na.ind)

    names(p.raw) <- names(p.adj.fdr) <- rownames(R2.m) <- rownames(RSS.m) <- rownames(F0.m) <- row.names
    colnames(R2.m) <- colnames(F0.m) <- colnames(RSS.m) <- paste0("F", 1:feature.no)
    
    if (is.fwer) {
      p.adj.fwer <- perm.fwer.adj(F0, Fp)
      p.adj.fwer <- na.pad(p.adj.fwer, na.ind)
      names(p.adj.fwer) <- row.names
    }
    else {
      p.adj.fwer <- NULL
    }
    cat('*')
    coef.list <- NULL
    for (j in 1:feature.no) {
      coef.list[[j]] <- solve(t(M01) %*% M01) %*% t(M01) %*%
        t(feature.list[[j]])
    }
    
    if (verbose)
      cat("Completed!\n")
    return(
      list(
        call = this.call,
        R2 = R2.m,
        F0 = F0.m,
        RSS = RSS.m,
        df.model = df.model,
        df.residual = df.residual,
        coef.list = coef.list,
        p.raw = p.raw,
        p.adj.fdr = p.adj.fdr,
        p.adj.fwer = p.adj.fwer
      )
    )
  }


DiSC <- function (data.mat, cell.ind, metadata, outcome, covariates = NULL,
                  cell.id = "cell_id", individual.id = "individual", 
                  perm.no = 999, features = c('prev', 'nzm', 'nzsd'), 
                  verbose = TRUE, sequencing.data = TRUE) {
  #data.mat: data matrix. Genes/variables are in rows, cells are in columns, 
  # colnames are `cell_id`
  #cell.ind: a data frame linking `cell_id` to `individual_id`
  #metadata: a dataframe including an 'individual_id', an outcome of interest, 
  # and covariates
  #cell.id: the variable name of cell ids in `cell.ind`
  #individual.id: the variable name of the individual id variables in `cell.ind`
  # and `metadata`
  #perm.no: number of permutations
  #features: feature extracted from the data matrix for each individual
  #verbose: do we need to print the progresses?
  #sequencing.data: is the data.mat a sequencing data matrix 
  # (e.g. scRNA sequencing data/count data)?
  
  require(matrixStats)
  
  inds <- unique(cell.ind[[individual.id]])
  #match(target, df$name), order metadata according to inds 
  metadata <- metadata[base::match(inds, metadata[[individual.id]]), ]
  
  cell.list <- list()
  for (ind in inds) {
    cell.list[[ind]] <- cell.ind[cell.ind[[individual.id]] == ind, ][[cell.id]]
  }

  if(sequencing.data){
    depth <- colSums(data.mat)
    data.mat <- t(t(data.mat) / depth)
    log_md_depth <- numeric(length = length(inds))
    names(log_md_depth) <- inds
    for(ind in inds) 
      log_md_depth[ind] <- log(median(depth[cell.list[[ind]]]))
    metadata$log_md_depth <- log_md_depth
  } else {
    depth <- log_md_depth <- NULL
  }
  
  yy.list <- list()
  for (feature in features) {
    yy.list[[feature]] <- matrix(NA, nrow(data.mat), length(inds), 
                                 dimnames = list(rownames(data.mat), inds))
  }
  # feature extraction
  for (ind in inds) {
    expr <- data.mat[, cell.list[[ind]]]
    if('prev' %in% features){
      prv <- rowMeans(expr == 0)
      prv <- log(prv / (1 - prv))
      prv[prv > 7] <- 7 
      prv[prv < -7] <- -7 #prv < 0.001
      yy.list[['prev']][, ind] <- prv 
    }
    if('nzm' %in% features){
      yy.list[['nzm']][, ind] <- sqrt(rowMeans(expr))
    }
    if('nzsd' %in% features){
      yy.list[['nzsd']][, ind] <- 
        sqrt(apply(expr, MARGIN = 1, 
                   FUN = function(x){
                     temp = sd(x[x > 0])
                     if(is.na(temp)) return(0) else
                       return(temp)}
      ))}
    if('sd' %in% features){
      yy.list[['sd']][, ind]   <- sqrt(rowSds(expr)) 
    }
    if('nzm^1' %in% features){
      yy.list[['nzm^1']][, ind] <- rowMeans(expr)
    }
    if('nzsd^1' %in% features){
      yy.list[['nzsd^1']][, ind] <- 
        apply(expr, MARGIN = 1, 
              FUN = function(x){
                temp = sd(x[x > 0])
                if(is.na(temp)) return(0) else
                  return(temp)}
        )
    }
    # pretty balanced. Some genes have highest Fs in "prev", while others in 
    #"nzm" or "nzsd" 
  }
  
  obj <- basic.func(metadata, yy.list, 
                    grp.name = outcome,
                    adj.name = c(covariates, ifelse(sequencing.data, 
                                                    "log_md_depth", NULL)),
                    stats.combine.func = base::max,
                    perm.no = perm.no, verbose = verbose)
  
  return(obj)
}