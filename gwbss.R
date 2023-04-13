gpower <- function(A, power = -1){
  
  # init tol
  tol = sqrt(.Machine$double.eps)
  
  # eigen decomposition
  eig.decompt = eigen(A, symmetric = TRUE)
  eig.values = eig.decompt$values
  thresh = max(tol * eig.values[1L], 0)
  positives = eig.values > thresh
  
  ddiag = eig.values[positives]^power
  val = eig.decompt$vectors[, positives] %*% (t(eig.decompt$vectors[, positives]) * ddiag)
  
  return(list(mat = val, vals = eig.values, valspower = ddiag, orthmat = eig.decompt$vectors[, positives])) 
}

# spatial version of first scatter - dat is spatial field, wt the current set of weights 
S1 <- function(dat, wt, eps){
  
  # only positive weights
  use <- abs(wt) > eps
  wt  <- wt[use]
  dat_w <- dat[use, ]
  dat_w <- sweep(dat_w, 1, sqrt(wt / sum(wt)), "*")
  
  # calculate scatter
  s <- t(dat_w) %*% dat_w
  
  return(s)
}




# sbss scatter with square root weights transformed data 
S2_wcov <- function(dat, k_mat, wt, eps){
  
  # only positive weights
  n <- nrow(dat)
  use <- abs(wt) > eps
  wt  <- wt[use]
  dat_w <- dat[use, ]
  dat_w <- sweep(dat_w, 1, sqrt(wt / sum(wt)), "*")
  
  # compute spatial scatter
  s <- t(dat_w) %*% k_mat[use, use] %*% dat_w 
  
  return(s)
}

# sbss scatter with square root weights transformed data 
S2_wvario <- function(dat, k_mat, wt, eps){
  
  # only positive weights
  n <- nrow(dat)
  use <- abs(wt) > eps
  rhot  <- k_mat[use, use] %*% sqrt(wt / sum(wt))
  dat_w <- dat[use, ]
  dat_w <- sweep(dat_w, 1, rhot * sqrt(wt / sum(wt)), "*")
  
  # compute spatial scatter
  s <- 2 * t(dat_w) %*% dat_w -
    2 * S2_wcov(dat, k_mat, wt, eps)

  return(s)
}

# sbss scatter with square root weights transformed data 
S2_wgraph <- function(dat, k_mat, xi, wt, eps){
  
  # only positive weights
  n <- nrow(dat)
  use <- abs(wt) > eps
  wt  <- wt[use]
  dat_w <- dat[use, ]
  dat_w <- sweep(dat_w, 2, xi, "-")
  dat_w <- sweep(dat_w, 1, sqrt(wt / sum(wt)), "*")
  
  # compute spatial scatter
  s <- t(dat_w) %*% k_mat[use, use] %*% dat_w 
  
  return(s)
}

# spatial version of matrix of fourth moments 
S2_sfobi <- function(dat, wt, eps){
  
  # only positive weights
  use <- abs(wt) > eps
  wt  <- wt[use]
  dat <- dat[use, ]
  dat_w <- sweep(dat, 1, sqrt(rowSums(dat ^ 2) * wt / sum(wt)), "*")
  
  # calculate scatter
  s <- t(dat_w) %*% dat_w

  return(s)
}

gwbss <- function(x, coords, bw, spatial_mean = TRUE, S2_type = c("wcov", "wvario", "wgraph", "sfobi"),
                  field_order = c("gwpca", "gwica"), kernel_type = "ball", kernel_parameters = bw) {
  # init
  S2_type <- match.arg(S2_type)
  field_order <- match.arg(field_order) 
  n <- nrow(coords)
  p <- ncol(x)
  eps <- .Machine$double.eps
  x_cen <- matrix(0, n, p)
  
  # weights
  weights <- GWmodel::gw.weight(as.matrix(dist(coords)), bw = bw, kernel = "gaussian", adaptive = FALSE)

  # center data 
  if(spatial_mean){
    x_cen <- lapply(seq_len(n),function(idx){ 
                  wt  <- weights[idx, ]
                  use <- abs(wt) > eps
                  wt  <- wt[use]
                  xwc <- colSums(sweep(x[use, ], 1, wt / sum(wt) , "*"))
                  x[idx,] - xwc
              })
    x_cen <- do.call("rbind", x_cen)
  } else {
    x_cen <- scale(x, center = colMeans(x), scale = FALSE)  # global
  }
  cen <- x - x_cen
  
  # compute S1
  S1_invsq <- lapply(seq_len(n), function(idx) gpower(S1(x_cen, weights[idx, ], eps), - 1 / 2))
  x_wh <- lapply(seq_len(n), function(idx) x_cen %*% S1_invsq[[idx]]$mat)
  
  # compute and diagonalize second scatter
  k_mat <- SpatialBSS::spatial_kernel_matrix(coords, kernel_type = kernel_type,
                                             kernel_parameters = kernel_parameters)[[1]]
  if (S2_type == "wcov") {
    S2_diago <- lapply(seq_len(n), function(idx) { 
      s2 <- list(S2_wcov(x_wh[[idx]], k_mat, weights[idx, ], eps)) 
      attr(s2, "lcov") <- "lcov"
      SpatialBSS:::diag_scatters(s2, ordered = TRUE)
    })
  } else if (S2_type == "wvario") {
    S2_diago <- lapply(seq_len(n), function(idx) { 
      s2 <- list(S2_wvario(x_wh[[idx]], k_mat, weights[idx, ], eps)) 
      attr(s2, "lcov") <- "lcov"
      SpatialBSS:::diag_scatters(s2, ordered = FALSE)
    })
    S2_diago <- lapply(S2_diago, function(S2) { 
      ord <- order(diag(S2$d), decreasing = FALSE)
      S2$u <- S2$u[, ord]
      S2$d <- S2$d[ord, ord]
      S2
    })
  } else if (S2_type == "wgraph") {
    S2_diago <- lapply(seq_len(n), function(idx) { 
      s2 <- list(S2_wgraph(x_wh[[idx]], k_mat, x_wh[[idx]][idx, ], weights[idx, ], eps)) 
      attr(s2, "lcov") <- "lcov"
      SpatialBSS:::diag_scatters(s2, ordered = FALSE)
    })
    S2_diago <- lapply(S2_diago, function(S2) { 
      ord <- order(diag(S2$d), decreasing = FALSE)
      S2$u <- S2$u[, ord]
      S2$d <- S2$d[ord, ord]
      S2
    })
  } else if (S2_type == "sfobi") {
    S2_diago <- lapply(seq_len(n), function(idx) { 
      s2 <- list(S2_sfobi(x_wh[[idx]], weights[idx, ], eps)) 
      attr(s2, "lcov") <- "lcov"
      SpatialBSS:::diag_scatters(s2, ordered = TRUE)
    })
  }
 
  # order of latent components
  if (field_order == "gwpca") {
    S2_diago <- lapply(seq_len(n), function(idx){ 
      ord <- order(S1_invsq[[idx]]$vals, decreasing = TRUE)
      S2_diago[[idx]]$u <- S2_diago[[idx]]$u[, ord]
      S2_diago[[idx]]$d <- S2_diago[[idx]]$d[ord, ord]
      S2_diago[[idx]]
    })
  }

  # compute loadings 
  loadings <- lapply(seq_len(n), function(idx){ 
                     lo <- crossprod(S2_diago[[idx]]$u, S1_invsq[[idx]]$mat)
                     lo <- sweep(lo, 1, sign(lo[, 1]), "*")
                     lo
               })
  loadings <- abind(loadings,along = 3)
  loadings <- aperm(loadings,perm = c(3,1,2))
  
  # scores
  scores <- lapply(seq_len(n), function(idx) x_cen[idx, ] %*% t(loadings[idx, , ]))
  scores <- do.call("rbind", scores)

  # eigvalues
  d <- lapply(seq_len(n), function(idx) as.vector(diag(S2_diago[[idx]]$d)))
  d <- do.call("rbind", d)
  
  return(list(s = scores, loadings = loadings, d = d, gwmeans = cen))
}


