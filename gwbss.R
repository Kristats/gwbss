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


# covariance scatter to s'
S2_scov <- function(dat, xi, wt, eps){
  
  # only positive weights
  use <- abs(wt) > eps
  wt  <- wt[use]
  dat_w <- dat[use, ]
  dat_w <- sweep(dat_w, 1, wt / sum(wt), "*")
  
  # calculate scatter
  s <- do.call(rbind,
                lapply(xi, function(xi) colSums(xi * dat_w)))
  s <- (s + t(s)) / 2

  return(s)
}

# variogram scatter with center s'
S2_vario <- function(dat, xi, wt, eps){
  
  # only positive weights
  use <- abs(wt) > eps
  wt  <- wt[use]
  dat_w <- dat[use, ]
  dat_w <- - sweep(dat_w, 2, xi, "-")
  dat_w <- sweep(dat_w, 1, sqrt(wt / sum(wt)), "*")
  
  # calculate scatter
  s <- t(dat_w) %*% dat_w 
  
  return(s)
}

# sbss scatter with square root weights transformed data 
S2_sbssw <- function(dat, k_mat, wt, eps){
  
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

# spatial version of matrix of fourth moments 
S2_sfobi <- function(dat, wt, eps){
  
  # only positive weights
  use <- abs(wt) > eps
  wt  <- wt[use]
  dat <- dat[use, ]
  dat_w <- sweep(dat, 1, sqrt(wt / sum(wt)), "*")
  
  # calculate scatter
  s <- crossprod(sweep(dat_w, 1, rowSums(dat ^ 2), "*"), 
                dat_w)

  return(s)
}

gwbss <- function(x, coords, bw, kernel_type = "ball", kernel_parameters = bw,
                   spatial_mean = TRUE, S2_type = c("scov", "vario", "sbssw", "sfobi"),
                  field_order = c("gwpca", "gwica")) {
  # init
  S2_type <- match.arg(S2_type)
  field_order <- match.arg(field_order) 
  n <- nrow(coords)
  p <- ncol(x)
  eps <- .Machine$double.eps
  x_cen <- matrix(0,n,p)
  
  # weights
  weights <- GWmodel::gw.weight(as.matrix(dist(coords)), bw = bw, kernel = "gaussian", adaptive = FALSE)

  # center data 
  if(spatial_type$spatial_mean){
    x_cen <- lapply(seq_len(n),function(n.idx){ 
                  wt  <- weights[n.idx, ]
                  use <- abs(wt) > eps
                  wt  <- wt[use]
                  xwc <- colSums(sweep(x[use, ], 1, wt / sum(wt) , "*"))
                  x[n.idx,] - xwc
              })
    x_cen <- do.call("rbind",x_cen)
  } else {
    x_cen <- scale(x, center = colMeans(x), scale = FALSE)  # global
  }
  cen <- x - x_cen
  
  # compute S1
  S1_invsq <- lapply(seq_len(n), function(n.idx) gpower(S1(x_cen,weights[n.idx, ]), - 1 / 2))

  # compute and diagonalize second scatter
  if (S2_type == "scov") {
    stop("Not implemented!")
  } else if(S2_type == "vario") {
    stop("Not implemented!")
  } else if(S2_type == "sbssw") {
    k_mat <- as.matrix(SpatialBSS::spatial_kernel_matrix(coords, kernel_type = kernel_type, 
                                                         kernel_parameters = kernel_parameters)[[1]])
    S2s <- lapply(seq_len(n), function(n.idx){ 
      x_whtd <- x_whtd_s[[n.idx]]
      s2l    <- list(S2spatial(x_whtd, k_mat, weights[n.idx,])) 
      attr(s2l, "lcov") <- "lcov"
      SpatialBSS:::diag_scatters(s2l, ordered = FALSE)
    })
  } else if (S2_type == "sfobi") {
    stop("Not implemented!")
  }
 
  # compute loadings 
  loadings <- lapply(seq_len(n), function(n.idx){ 
                     lo <- crossprod(S2s[[n.idx]]$u, S1s_sqinv[[n.idx]]$mat)
                     lo <- lo[, order(S1s_sqinv[[n.idx]]$vals, decreasing = TRUE)]
                     lo <- sweep(lo, 1, sign(lo[1, ]), "*")
                     lo
               })
  loadings <- abind(loadings,along = 3)
  loadings <- aperm(loadings,perm = c(3,1,2))
  
  # scores
  scores <- lapply(seq_len(n), function(n.idx){   x_cen[n.idx,] %*% t(loadings[n.idx,,])  })
  scores <- do.call("rbind",scores)

  # eigvalues
  d <- lapply(seq_len(n),function(n.idx){ as.vector(diag(S2s[[n.idx]]$d)) })
  d <- do.call("rbind",d)
  
  return(list(s = scores, loadings = loadings, d = d, gwmeans = cen))
}


