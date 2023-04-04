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


# gwbss <- function(x, coords, bw, scatter = c("1", "2"), kernel_type = "ball", kernel_parameters = bw) {
#   scatter <- match.arg(scatter)
#   
#   n <- nrow(coords)
#   p <- ncol(x)
#   
#   weights <- GWmodel::gw.weight(as.matrix(dist(coords)), bw = bw, kernel = "gaussian", adaptive = FALSE)
#   
#   loadings <- array(0, dim = c(n, p, p))
#   d <- matrix(0, nrow = n, ncol = p)
#   scores <- matrix(0, nrow = n, ncol = p)
#   gwmeans <- matrix(0, nrow = n, ncol = p)
#   
#   
#   if (scatter == "2") {
#     k_mat <- SpatialBSS::spatial_kernel_matrix(coords, kernel_type = kernel_type, 
#                                                kernel_parameters = kernel_parameters)
#   }
#   
#   
#   # cen <- lapply(seq_len(n),function(n.idx){ 
#   #   wt  <- weights[n.idx, ]
#   #   use <- abs(wt) > 0
#   #   wt  <- wt[use]
#   #   xwc <- colSums(sweep(x[use, ], 1, wt / sum(wt) , "*"))
#   #   xwc
#   # })
#   # cen = do.call("rbind",cen)
#   # 
#   
#   # local bss
#   for (idx in 1:n) {
#     # prepare weights and data
#     wt <- weights[idx, ]
#     use <- wt > 0
#     wt <- wt[use]
#     dat <- x[use, ]
# 
#     # center data with gwmean
#     gwmean <- colSums(sweep(dat, 1, wt / sum(wt) , "*")) #cen[idx,] #
#     #dat_cen <- sweep(dat, 2, gwmean, "-")
#     #dat_cen = dat - cen
#     dat_cen_w <- sweep(dat_cen, 1, sqrt(wt), "*")
#     
#     # white data with gwcov
#     s1 <- crossprod(dat_cen_w) / sum(wt)
#     s1_inv_sqrt <- gpower(s1, - 1 / 2)
# 
#     # compute scatters 
#     if (scatter == "1") {
#       dat_cen_w <- sweep(dat_cen, 1, wt, "*")
#       s2 <- do.call(rbind,
#                     lapply(x[idx, ] - gwmean, function(xi) colSums(xi * dat_cen_w)))
#       s2 <- (s2 + t(s2)) / 2 / sum(wt)
#       s2 <- list(s1_inv_sqrt$mat %*% s2 %*% s1_inv_sqrt$mat)
#     } else if (scatter == "2") {
#       k <- list(k_mat[[1]][use, use])
#       s2 <- SpatialBSS::local_covariance_matrix(dat_cen_w, k, center = FALSE)
#       s2 <- list(s1_inv_sqrt$mat %*% (s2[[1]] * n / sum(wt)) %*% s1_inv_sqrt$mat)
#     }
#     attr(s2, "lcov") <- "lcov"
#     
#     # diag scatter
#     s_diag <- SpatialBSS:::diag_scatters(s2, ordered = FALSE)
# 
#     # scores, loadings, d
#     gwmeans[idx, ] <- gwmean
#     d[idx, ] <- diag(s_diag$d)
#     l <- crossprod(s_diag$u, s1_inv_sqrt$mat)
#     l <- l[, order(s1_inv_sqrt$vals, decreasing = TRUE)]
#     l <- sweep(l, 1, sign(l[1, ]), "*")
#     loadings[idx,, ] <- l
#     scores[idx, ] <- (x[idx, ] - gwmean) %*% t(loadings[idx,, ]) 
#   } 
#   
#   return(list(s = scores, loadings = loadings, d = d, gwmeans = gwmeans))
# }



gwbss <- function(x, coords, bw, kernel_type = "ball", kernel_parameters = bw,
                   spatial_type = list(spatial_mean = TRUE, spatial_scores = TRUE)) {

  # spatial version of first scatter - dat is spatial field, wt the current set of weights 
  S1spatial <- function(dat,wt){
    
    # only positive weights
    use <- abs(wt) > eps
    wt  <- wt[use]
    dat_w <- dat[use, ]
    dat_w <- sweep(dat_w, 1, sqrt(wt / sum(wt)), "*")
    
    # calculate scatter
    s1  = t(dat_w) %*% dat_w
    
    return(s1)
  }
  
  
  # spatial version of second scatter - dat is spatial field, wt the current set of weights 
  S2spatial <- function(dat,k_mat,wt){
    
    # only positive weights
    n <- nrow(dat)
    use <- abs(wt) > eps
    wt  <- wt[use]
    dat_w <- dat[use, ]
    dat_w <- sweep(dat_w, 1, sqrt(wt / sum(wt)), "*")

    # compute spatial scatter
    #s2 <- SpatialBSS::local_covariance_matrix(x = dat_w, kernel_list = list(k_mat[use, use]), lcov  = 'lcov', center = FALSE)
    #s2 <- s2[[1]] * n
    s2 <- t(dat_w) %*% k_mat[use,use] %*% dat_w # the same as the above 
    
    return(s2)
  }

  
  # weights
  weights <- GWmodel::gw.weight(as.matrix(dist(coords)), bw = bw, kernel = "gaussian", adaptive = FALSE)
  
  # init
  n <- nrow(coords)
  p <- ncol(x)
  eps <- 0#.Machine$double.eps
  x_cen <- matrix(0,n,p)
  
  # initialize spatial kernel
  k_mat <- as.matrix(SpatialBSS::spatial_kernel_matrix(coords, kernel_type = kernel_type, 
                                               kernel_parameters = kernel_parameters)[[1]])
  
  
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
  }else if(!spatial_type$spatial_mean){
    x_cen <- scale(x, center = colMeans(x), scale = FALSE)  # global
  }
  cen <- x - x_cen
  
  
  # compute first scatters with field and take sqrt inverse
  S1s_sqinv <- lapply(seq_len(n), function(n.idx){ gpower(S1spatial(x_cen,weights[n.idx,]),-1/2) })
  
  # whiten data 
  x_whtd_s <- lapply(seq_len(n), function(n.idx){ x_cen %*% S1s_sqinv[[n.idx]]$mat })
  
  # compute second scatters with field and 
  if(spatial_type$spatial_scores){
     S2s <- lapply(seq_len(n), function(n.idx){ 
                   x_whtd <- x_whtd_s[[n.idx]]
                   s2l    <- list(S2spatial(x_whtd, k_mat, weights[n.idx,])) 
                   attr(s2l, "lcov") <- "lcov"
                   SpatialBSS:::diag_scatters(s2l, ordered = FALSE)
             })
  }else{
    stop("Not implemented!")
    x_whtd <- lapply(seq_len(n), function(n.idx) x_whtd_s[[n.idx]][n.idx, ])
    S2s <- S2spatial(x_whtd, k_mat, rep(1,n))
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


