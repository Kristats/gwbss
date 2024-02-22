library("plyr")


# function to return maximum index of absolutes of each column of a matrix
rowMaxAbs <- function(A, vrn)  apply(A, 1, function(u){ vrn[which.max(abs(u))] })
rowOrderAbs <- function(A, vrn) t(apply(A, 1, function(u){ vrn[order(abs(u), decreasing = TRUE)] }))

# generalized inverse function 
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

# function that works with arrays - data should be centered already!
S1_arr <-  function(dat, wt, eps){
  
  # only positive weights
  nl <- dim(dat)[1]
  pl <- dim(dat)[2]
  rel <- dim(dat)[3]
  dat_w <- sweep(dat, 1, sqrt(wt / sum(wt)), "*")
  dat_w <- apply(dat_w, 2, c)
  
  # calculate scatter
  s <- 1 / rel * t(dat_w) %*% dat_w
  
  return(s)
}


# sbss scatter with square root weights transformed data 
S2_wcov <- function(dat, k_mat, wt, eps){
  
  # only positive weights
  n <- nrow(dat)
  use <- abs(wt) > eps
  wt  <- wt[use]
  dat_w <- dat[use, ]
  wt_s  <- wt / sum(wt)
  dat_w <- sweep(dat_w, 1, sqrt(wt_s), "*")
  
  # compute spatial scatter
  s <- t(dat_w) %*% k_mat[use, use] %*% dat_w 
  
  return(s)
}

# function that works with arrays - data should be centered already!
S2_wcov_arr <-  function(dat, k_mat, wt, eps){

  # only positive weights
  nl <- dim(dat)[1]
  pl <- dim(dat)[2]
  rel <- dim(dat)[3]
  wt_s  <- wt / sum(wt)
  dat_w <- lapply(1:rel, function(idx){ sweep(dat[,,idx],1,sqrt(wt_s), "*") }) #sweep(dat, 1, sqrt(wt_s), "*")
  s     <- lapply(1:rel, function(idx){ t(dat_w[[idx]]) %*% k_mat %*% dat_w[[idx]]  })
  s <- 1 / rel * do.call("+",s)

  return(s)
}



# sbss scatter with square root weights transformed data 
S2_wvario <- function(dat, k_mat, wt, eps){
  
  # only positive weights
  n <- nrow(dat)
  wt_s  <- wt / sum(wt)
  rhot  <- as.vector(k_mat %*% sqrt(wt_s))
  dd    <- as.vector(sqrt(wt_s)) * rhot
  sc    <- sum(rhot * sqrt(wt_s)) # scale
  dat_w <- dat

  # compute spatial scatter
  s <- 2 / sc * (t(dat_w) %*% sweep(dat_w, 1, dd, "*") - 
      S2_wcov(dat, k_mat, wt, eps))

  return(s)
}

# function that works with arrays - data should be centered already!
S2_wvario_arr <- function(dat, k_mat, wt, eps){
  
  # only positive weights
  rel <- dim(dat)[3]
  pl <- dim(dat)[2]
  dat_w <- dat
  wt_s  <- wt / sum(wt)
  rhot  <- as.vector(k_mat %*% sqrt(wt_s))
  dd    <- as.vector(sqrt(wt_s)) * rhot
  sc    <- rel * sum(rhot * sqrt(wt_s))  # rel * sqrt(wt_s) %*% k_mat %*% sqrt(wt_s)
  if(sc == 0){ return( matrix(0,pl,pl))}

  #s     <- lapply(1:rel, function(idx){  
  #                         dat_loc <- sweep(dat_w[,,idx], 1, dd, "*")
  #                         t(dat_loc) %*% dat_loc })
  #s     <-  do.call("+",s)
  
  s     <- lapply(1:rel, function(idx){ t(dat_w[,,idx]) %*%  sweep(dat_w[,,idx], 1, dd, "*") })
  s     <-  do.call("+",s)
  
  # compute spatial scatter
  s <- 2 / sc * (s - rel * S2_wcov_arr(dat, k_mat, wt, eps))

  return(s)
}


# sbss scatter with square root weights transformed data 
S2_wgraph <- function(dat, k_mat, xi, wt, eps){
  
  # only positive weights
  n <- nrow(dat)
  #use <- abs(wt) > eps
  wt_s  <- wt / sum(wt)
  sc    <- sum(sqrt(wt_s) %*% k_mat %*% sqrt(wt_s))
  dat_w <- dat
  dat_w <- sweep(dat_w, 2, xi, "-")
  dat_w <- sweep(dat_w, 1, sqrt(wt_s), "*")
  
  # compute spatial scatter
  s <- 1 / sc * t(dat_w) %*% k_mat %*% dat_w 
  
  return(s)
}

# function that works with arrays - data should be centered already!
S2_wgraph_arr <- function(dat, k_mat, s_idx, wt, eps){
  
  # only positive weights
  rel   <- dim(dat)[3]
  pl   <- dim(dat)[2]
  wt_s  <- wt / sum(wt)
  U_cen <- dat[s_idx,,]
  sc    <- rel * sum(sqrt(wt_s) %*% k_mat %*% sqrt(wt_s))
  if(sc == 0){ return( matrix(0,pl,pl))}
  
  # subtract "center" 
  dat_w <- sweep(dat, MARGIN = 2:3, STATS = U_cen, FUN = "-")
 
  # scatter
  s <- lapply(1:rel, function(idx){ 
                                datloc <- sweep(dat_w[,,idx], 1, sqrt(wt_s), "*") # multiply by sqrt of weights
                                t(datloc) %*% k_mat %*% datloc })                 # make scatter
  s     <- 1 / sc * do.call("+",s)

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

gwbss <- function(x, coords, spatial_mean = TRUE, weights_adj, S2_type = c("wcov", "wvario", "wgraph", "sfobi"),
                  field_order = c("gwpca", "gwica"), kernel_type = "ball", kernel_parameters = bw) {
  # init
  S2_type <- match.arg(S2_type)
  field_order <- match.arg(field_order) 
  n <- nrow(coords)
  p <- ncol(x)
  eps <- .Machine$double.eps
  x_cen <- matrix(0, n, p)
  
  # names
  varnames <- if(!is.null(colnames(x))){ colnames(x) }else{ as.character(1:p) }

  # center data 
  if(spatial_mean){
    cen <- lapply(seq_len(n),function(idx){ 
                  wt  <- weights_adj[idx, ]
                  use <- abs(wt) > eps
                  wt  <- wt[use]
                  xwc <- colSums(sweep(x[use, ], 1, wt / sum(wt) , "*"))
                  xwc
              })
    cen <- do.call("rbind", cen)  
  } else {
    cen <- t(matrix(colMeans(x), ncol = n, nrow = p))  # global
  }
  
  # center data
  x_cen <- lapply(seq_len(n), function(idx){ scale(x, center = cen[idx,], scale = FALSE)  })
  
  ### compute S1
  # sqrt inverse of first scatter
  S1_invsq <- lapply(seq_len(n), function(idx){ gpower(S1(x_cen[[idx]], weights_adj[idx, ], eps), - 1 / 2)  })
  
  ## whiten data 
  x_wh <- lapply(seq_len(n), function(idx) x_cen[[idx]] %*% S1_invsq[[idx]]$mat)
  
  # compute and diagonalize second scatter
  k_mat <- SpatialBSS::spatial_kernel_matrix(coords, kernel_type = kernel_type,
                                             kernel_parameters = kernel_parameters)[[1]]
  
  
  

   if (S2_type == "wcov") {
     S2_diago <- lapply(seq_len(n), function(idx) {
       s2 <- list(S2_wcov(x_wh[[idx]], k_mat, weights_adj[idx, ], eps))
       attr(s2, "lcov") <- "lcov"
       SpatialBSS:::diag_scatters(s2, ordered = TRUE)
     })
   } else if (S2_type == "wvario") {
     S2_diago <- lapply(seq_len(n), function(idx) {
       s2 <- list(S2_wvario(x_wh[[idx]], k_mat, weights_adj[idx, ], eps))
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
       s2 <- list(S2_wgraph(x_wh[[idx]], k_mat, x_wh[[idx]][idx, ], weights_adj[idx, ], eps))
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
       s2 <- list(S2_sfobi(x_wh[[idx]], weights_adj[idx, ], eps))
       attr(s2, "lcov") <- "lcov"
       SpatialBSS:::diag_scatters(s2, ordered = TRUE)
     })
   }

  # # order of latent components
  # if (field_order == "gwpca") {
  #   S2_diago <- lapply(seq_len(n), function(idx){
  #     ord <- order(S1_invsq[[idx]]$vals, decreasing = TRUE)
  #     S2_diago[[idx]]$u <- S2_diago[[idx]]$u[, ord]
  #     S2_diago[[idx]]$d <- S2_diago[[idx]]$d[ord, ord]
  #     S2_diago[[idx]]
  #   })
  # }

   # compute unmixing matrices 
  unmixing <- lapply(seq_len(n), function(idx){ 
                     lo <- crossprod(S2_diago[[idx]]$u, S1_invsq[[idx]]$mat)
                     lo <- sweep(lo, 1, sign(lo[, 1]), "*")
                     lo
               })
  
  # scores
  scores <- lapply(seq_len(n), function(idx) x_cen[[idx]][idx,] %*% t(unmixing[[idx]]))
  scores <- do.call("rbind", scores)

  # compute scores of edges 
  scores_edges <- lapply(seq_len(n), function(idx)  x_cen[[idx]] %*% t(unmixing[[idx]]) )
  scores_edges <- abind(scores_edges,along = 3)
  scores_edges <- aperm(scores_edges,  perm = c(3,1,2))     # first dim is s, second s_tilde, third the variables
  
  
  # eigvalues
  d <- lapply(seq_len(n), function(idx) as.vector(diag(S2_diago[[idx]]$d)))
  d <- do.call("rbind", d)
  
  # loadings
  #unmixing <- abind(unmixing,along = 3)
  #unmixing <- aperm(unmixing,perm = c(3,1,2))
  
  # return winning loading variable name
  winningloading  <- lapply(1:n, function(idx){ rowMaxAbs(unmixing[[idx]], varnames) })
  winningloading <- do.call("rbind", winningloading)
    
  orderedloadings <- lapply(1:n, function(idx){ rowOrderAbs(unmixing[[idx]], varnames) })
  orderedloadings <- abind(orderedloadings,along = 3)
  orderedloadings <- aperm(orderedloadings,perm = c(3,1,2))
  
  return(list(s        = scores, 
              unmixing = loadings, 
              d        = d, 
              gwmeans = cen, 
              winningloading = winningloading,
              orderedloadings = orderedloadings))
}


# function to calculate global center
glb_center <- function(x){
   n <- nrow(x); p <- ncol(x)
   arm <- apply(x, 2, mean)
   cen <- t(matrix(arm, nrow = p, ncol = n))
   cen
}

# function to calculate global covariance of centered data
glb_cov <- function(x){
   n <- nrow(x)
   cen <- glb_center(x)[1,]
   xc <- sweep(x, MARGIN = 2, STATS = cen, FUN = "-") 
   S1_arr(xc, rep(1,n), eps)
}


spatial_center <- function(x, weights_adj){
  n <- nrow(x)
  cen <- lapply(seq_len(n),function(idx){ 
                  wt  <- weights_adj[idx, ]
                  xwc <- apply(x, c(2,3), function(u) u * wt / sum(wt))
                  xwc <- apply(xwc, c(2,3), sum)
                  xwc <- rowMeans(xwc)
              })
  cen <- do.call("rbind", cen)
  cen
}

# graph weighted bss 
# if k_mat is provided then Mittinen scatter is used - weights are overwritten 
graphbss <- function(x, weights_adj, dist_neighbors = 1,  
                     spatial_mean = TRUE, spatial_s1 = TRUE, S2_type = c("wcov", "wvario", "wgraph"),  k_mat = NULL, mtf = FALSE,
                     geographic = list(flag = FALSE, coords = NULL, kernel_type = "ball", diam = NULL),
                     MDIs = TRUE){
  
  
  # check if x is in array form
  ns <- if(length(dim(x)) == 2){ 
    ns = c(dim(x),1) 
  }else if(length(dim(x)) == 3){
    ns = dim(x)
  }else{ 
    stop("array must be supplied") 
  }
  x <-  array(x, dim = ns)
  
  if(spatial_s1){
    if(!spatial_mean){ 
      warning("setting spatial mean to true as spatial_s1 true")
      spatial_mean = TRUE
      }
  }
  
  
  # init
  S2_type <- match.arg(S2_type)
  n <- dim(x)[1]
  p <- dim(x)[2]
  re <- dim(x)[3]
  S1s <- NULL
  diag(weights_adj) <- 0
  eps <- .Machine$double.eps
  x_cen <- array(0, dim = c(n, p, re))
  
  # names
  varnames <- if(!is.null(colnames(x))){ colnames(x) }else{ as.character(1:p) }

  # center data 
  if(spatial_mean){
    cen <- spatial_center(x, weights_adj)
  } else{
    cen <- glb_center(x)
  }

  # centered signals - for every node we return again an array
  x_cen <- lapply(seq_len(n), function(idx){ 
    sweep(x, MARGIN = 2, STATS = cen[idx,], FUN = "-") 
  })   

  ### compute S1
  # sqrt inverse of first scatter
  if(spatial_s1){
      S1s      <- lapply(seq_len(n), function(idx) S1_arr(x_cen[[idx]], weights_adj[idx, ], eps))
      S1_invsq <- lapply(seq_len(n), function(idx) gpower(S1s[[idx]], - 1 / 2))  
  }else{
      s1    <- glb_cov(x)
      s1inv <- gpower(s1, -1/2)
      S1_invsq <- lapply(seq_len(n), function(idx) s1inv)  
  }

  ## whiten data 
  x_wh <- lapply(seq_len(n), function(idx){
    Aloc <- aaply(x_cen[[idx]], 3 , function(A){ A %*% S1_invsq[[idx]]$mat })
    Aloc <- if(re != 1){ aperm(Aloc,perm = c(2,3,1)) }else{ array(Aloc, dim = dim(x)) }
    Aloc
  })
  
  # compute gamma
  if(mtf & !is.null(k_mat)){
    weights_adj  <- matrix(1,n,n)
    S2_type      <- "wgraph"
    spatial_mean <- FALSE
    
  }else if(is.null(k_mat)){
    if(!geographic$flag){
      if(dist_neighbors < 1){stop("dist_neighbors must be at least 1")}
      graph.ob <- graph_from_adjacency_matrix(weights_adj, weighted = FALSE)
      dists <- distances(graph.ob, mode = "all")
      dists <- t(apply(dists, 1, function(u){
                  v <- rep(0,length(u))
                  ind <- which(u <= dist_neighbors) 
                  v[ind] <- 1
                  v
                }))
      k_mat <- dists
    }else if(geographic$flag){
      if(is.null( geographic$coords) | is.null(geographic$diam)){ stop("coordinates and diameter must be supplied")}
      k_mat <- SpatialBSS::spatial_kernel_matrix(geographic$coords, kernel_type = geographic$kernel_type,
                                                 kernel_parameters = geographic$kernel_parameters)[[1]]
    }else{
      stop("geographic must be either TRUE or FALSE")
    }
  }


  # compute second scatter
  if (S2_type == "wcov") {
    S2_diago <- lapply(seq_len(n), function(idx) { 
      s2 <- list(S2_wcov_arr(x_wh[[idx]], k_mat, weights_adj[idx, ], eps)) 
      attr(s2, "lcov") <- "lcov"
      SpatialBSS:::diag_scatters(s2, ordered = TRUE)
    })
  } else if (S2_type == "wvario") {
    S2_diago <- lapply(seq_len(n), function(idx) { 
      s2 <- list(S2_wvario_arr(x_wh[[idx]], k_mat, weights_adj[idx, ], eps)) 
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
      s2 <- list(S2_wgraph_arr(x_wh[[idx]], k_mat, s_idx = idx, weights_adj[idx, ], eps)) 
      attr(s2, "lcov") <- "lcov"
      SpatialBSS:::diag_scatters(s2, ordered = FALSE)
    })
    S2_diago <- lapply(S2_diago, function(S2) { 
      ord <- order(diag(S2$d), decreasing = FALSE)
      S2$u <- S2$u[, ord]
      S2$d <- S2$d[ord, ord]
      S2
    })
  }

  # compute loadings 
  unmixing <- lapply(seq_len(n), function(idx){ 
                     lo <- crossprod(S2_diago[[idx]]$u, S1_invsq[[idx]]$mat)
                     lo <- sweep(lo, 1, sign(lo[, 1]), "*")
                     lo })

  # compute scores of edges 
  scores_edges <- lapply(seq_len(n), function(idx){  
                           scl <- aaply(x_cen[[idx]], 3, function(A){ A %*% t(unmixing[[idx]]) })
                           scl <- if(re != 1){ aperm(scl,  perm = c(2,3,1)) }else{ scl }
                           scl })  # first dim is s, second s_tilde, third the variables

  # scores
  scores <- lapply(seq_len(n), function(idx) t(x_cen[[idx]][idx,,]) %*%  t(unmixing[[idx]])  )
  scores <- if(re != 1){ 
    scores <- abind(scores,along = 3)
    scores <- aperm(scores,perm = c(3,2,1))
  }else{ 
     scores <- do.call("rbind", scores)
     scores <- array(scores, dim = ns)
  }

  # eigvalues
  d <- lapply(seq_len(n), function(idx) as.vector(diag(S2_diago[[idx]]$d)))
  d <- do.call("rbind", d)
    
  # winning loading variable name
  winningloading  <- lapply(1:n, function(idx){ rowMaxAbs(unmixing[[idx]], varnames) })
  winningloading <- do.call("rbind", winningloading)
    
  # also ordered loadings 
  orderedloadings <- lapply(1:n, function(idx){ rowOrderAbs(unmixing[[idx]], varnames) })
  orderedloadings <- abind(orderedloadings,along = 3)
  orderedloadings <- aperm(orderedloadings,perm = c(3,1,2))
  
  # calculate deviation from global center
  glcen <- glb_center(x)
  if(!spatial_mean){
    dist_center <- apply(glcen - spatial_center(x, weights_adj), 1, function(u) sqrt(sum(u*u))) /
      (sqrt(sum(glcen*glcen))+1)
  }else{
    dist_center <- apply(glcen - cen, 1, function(u) sqrt(sum(u*u))) /
      (sqrt(sum(glcen*glcen))+1)
  }
  
  
  # deviation from global covariance
  if(!spatial_s1){
    S1s      <- lapply(seq_len(n), function(idx) S1_arr(x_cen[[idx]], weights_adj[idx, ], eps))
  }  
  S1_inv   <- lapply(seq_len(n), function(idx) gpower(S1s[[idx]], - 1))  
  glc <- glb_cov(x)
  MDss1 <- lapply(1:n, function(idx){ 
    MD(W.hat = S1_inv[[idx]]$mat, A = glc) 
    norm(diag(rep(1,p)) - S1_inv[[idx]]$mat %*% glc,"F") / (norm(glc, "F") + 1)
    })
  MDss1 <- do.call("c", MDss1) 
  
  # calculate minimum distance index for global versus local estimates
  if(MDIs){
    if(!( S2_type %in% c("wvario", "wgraph"))){ 
      stop("MD not available for S2_type wcov. Choose wvario or wgraph")}
    # global unmixing matrix
    weights_ones <- matrix(1, ncol = n, nrow = n)
    global  <- graphbss(x = x, 
                       weights_adj = weights_ones, 
                       mtf = mtf,
                       dist_neighbors  = dist_neighbors ,
                       spatial_mean = spatial_mean, 
                       spatial_s1 = spatial_s1,
                       S2_type = S2_type ,
                       geographic = geographic,
                       MDIs = FALSE)
    mixingmat <- gpower(global$unmixing[[1]], power = -1)$mat
    
    MDss2 <- lapply(1:n, function(idx){ 
      MD(W.hat = unmixing[[idx]], A = mixingmat) 
      norm(diag(rep(1,p)) - unmixing[[idx]] %*% mixingmat,"F") / (norm(mixingmat, "F") + 1)
      })
    MDss2 <- do.call("c", MDss2) 
  }else{ MDss2 = NULL }
  
  return(list(scores  = scores,
              scores_edges = scores_edges,
              unmixing = unmixing, 
              eigenval = d, 
              gwmeans  = cen,
              winningloading = winningloading,
              orderedloadings = orderedloadings,
              dist_center = dist_center,
              MDss1 = MDss1,
              MDss2 = MDss2))
}

# function to do eigendecomposition and apply a function to its eigenvalues
# so g(A)
gA <- function(A, g){
  g <- Vectorize(g)
  eig <- eigen(A)
  eig$vectors %*% sweep(t(eig$vectors), 1, g(eig$values), FUN = "*")
}


# graph weighted bss 
# if k_mat is provided then Mittinen scatter is used - weights are overwritten 
grabss <- function(x,  
                  weights_adj,
                  nbr_nns = 2,
                  spatial_mean = TRUE,
                  spatial_cov = TRUE, 
                  laplacian = TRUE,
                  spatial_laplacian = TRUE,
                  gfun = NULL){
  # init
  n <- nrow(x)
  p <- ncol(x)
  eps <- .Machine$double.eps
  x_cen <- matrix(0, n, p)
  
  # compute distances
  graph.ob <- graph_from_adjacency_matrix(weights_adj, weighted = TRUE)
  dists <- distances(graph.ob, mode = "all")
  rho_mat <- dists / max(dists)
  sig   <- stats::quantile(rho_mat, probs = 0.9)
  rho_mat <- exp(- rho_mat^2 / sig^2 * 1 / 2)
  
  # names
  varnames <- if(!is.null(colnames(x))){ colnames(x) }else{ as.character(1:p) }

  # center data 
  if(spatial_mean){
    cen <- lapply(seq_len(n),function(idx){ 
                  wt  <- rho_mat[idx, ]
                  xwc <- colSums(sweep(x, 1, wt / sum(wt) , "*"))
                  xwc
              })
    cen <- do.call("rbind", cen)  
  } else if(!spatial_mean | spatial_cov){
    cen <- t(matrix(colMeans(x), ncol = n, nrow = p))  # global
  }
  
  # center data
  x_cen <- lapply(seq_len(n), function(idx){ scale(x, center = cen[idx,], scale = FALSE)  })
  
  ### compute S1s
  if(spatial_cov){
    S1s <- lapply(seq_len(n), function(idx){ S1(x_cen[[idx]], weights_adj[idx, ], eps) })
  }else{
    cen_gl <- glb_center(x)
    s1  <- S1(x - cen_gl, rep(1,n), eps)
    S1s <- lapply(seq_len(n), function(idx){ s1 })
  }
  
  
  # sqrt inverse of first scatter
  S1_invsq <- lapply(seq_len(n), function(idx){ gpower(S1s[[idx]], - 1 / 2)  })
  
  ## whiten data 
  x_wh <- lapply(seq_len(n), function(idx) x_cen[[idx]] %*% S1_invsq[[idx]]$mat)
  
  # gamma matrix - how far we go from s  
  gamma_mat <- Matrix(1 * (dists <= nbr_nns))
  
  # compute L(s) or W(s)
  Localmats <-  if(spatial_laplacian){
                     if(!is.null(gfun)){ message("g can only be used with sptial_laplacian = FALSE")}
                     lapply(seq_len(n), function(idx){ 
                                W_loc <-  Matrix(outer(rho_mat[idx,], rho_mat[idx,], FUN = "*")) * gamma_mat
                                A_loc <- if(laplacian){ diag(rowSums(W_loc)) - W_loc }else{ W_loc }
                                Matrix(A_loc)
                            })  
                }else{
                  W_gl <-  rho_mat * gamma_mat
                  A_gl <- if(laplacian){ diag(rowSums(W_gl)) - W_gl }else{ W_gl }
                  A_gl <- if(!is.null(gfun)){ gA(A_gl, function(x){gfun(x)})}else{ A_gl }
                  lapply(seq_len(n), function(idx){ A_gl })
                }
  
  # compute scatters
  S2_diago <- lapply(1:n, function(idx){ 
                   s2 <- 1 / n * t(x_wh[[idx]]) %*% Localmats[[idx]] %*% x_wh[[idx]]
                   eig <- eigen(s2, symmetric = TRUE)
                   if(laplacian){ inds <- (p-1):1 }else{ inds <- p:1 }
                   list(u = eig$vectors[,inds], d = eig$values[inds])
                })

  # compute unmixing matrices 
  unmixing_local <- lapply(seq_len(n), function(idx){ 
                     lo <- crossprod(S2_diago[[idx]]$u, S1_invsq[[idx]]$mat)
                     lo <- sweep(lo, 1, sign(lo[, 1]), "*")
                     colnames(lo) <- varnames
                     lo
               })
  
  # scores
  scores_ss <- lapply(seq_len(n), function(idx) x_cen[[idx]][idx,] %*% t(unmixing_local[[idx]]))
  scores_ss <- do.call("rbind", scores_ss)

  # compute scores of edges 
  scores_edges  <- lapply(seq_len(n), function(idx)  x_cen[[idx]] %*% t(unmixing_local[[idx]]) )
  scores_global <- 1 / n * Reduce("+", scores_edges) #1 / n * aaply(scores_edges, c(2:3), sum)
  
  # bind
  scores_edges <- abind(scores_edges, along = 3)
  scores_edges <- aperm(scores_edges,  perm = c(3,1,2))     # first dim is s, second s_tilde, third the variables
  
  # global scores 
  #scores_global <- 1 / n * aaply(scores_edges, c(2:3), sum)
  
  # global unmixing
  unmixing_global <-  1 / n * Reduce("+", unmixing_local)
  colnames(unmixing_global) <- varnames
  
  # eigvalues
  d <- lapply(seq_len(n), function(idx) as.vector(diag(S2_diago[[idx]]$d)))
  d <- do.call("rbind", d)

  
  # return winning loading variable name
  winningloading  <- lapply(1:n, function(idx){ rowMaxAbs(unmixing_local[[idx]], varnames) })
  winningloading <- do.call("rbind", winningloading)
    
  orderedloadings <- lapply(1:n, function(idx){ rowOrderAbs(unmixing_local[[idx]], varnames) })
  orderedloadings <- abind(orderedloadings,along = 3)
  orderedloadings <- aperm(orderedloadings,perm = c(3,1,2))
  
  return(list(scores_ss  = scores_ss, 
              scores_edges = scores_edges,
              scores_global = scores_global,
              unmixing_local = unmixing_local, 
              unmixing_global = unmixing_global,
              eigenvalues = d, 
              centers = cen, 
              winningloading = winningloading,
              orderedloadings = orderedloadings))
}


