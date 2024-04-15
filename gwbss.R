library("plyr")


# function to return maximum index of absolutes of each column of a matrix
rowMaxAbs <- function(A, vrn)  apply(A, 1, function(u){ vrn[which.max(abs(u))] })
rowOrderAbs <- function(A, vrn) t(apply(A, 1, function(u){ vrn[order(abs(u), decreasing = TRUE)] }))

# functions to return just negative and positive loadings 
rowOrderPos <- function(A, vrn){t(apply(A, 1, function(u){ 
  #ind <- !is.na(ifelse(u>0,u,NA)); vrn[!ind] <- NA
  #vrn[ind] <- vrn[ind][order(u[ind], decreasing = TRUE)]
  ord <- order(abs(u), decreasing = TRUE)
  ind <- which(u[ord] > 0)
  vrn <- vrn[ord]
  vrn[-ind] <- NA
  vrn }))}
rowOrderNeg <- function(A, vrn){t(apply(A, 1, function(u){ 
  #ind <- !is.na(ifelse(u<0,u,NA)); vrn[!ind] <- NA
  #vrn[ind] <- vrn[ind][order(u[ind], decreasing = TRUE)] 
  ord <- order(abs(u), decreasing = TRUE)
  ind <- which(u[ord] < 0)
  vrn <- vrn[ord]
  vrn[-ind] <- NA
  vrn }))}



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

# function to do eigendecomposition and apply a function to its eigenvalues
# so g(A)
gA <- function(A, g){
  g <- Vectorize(g)
  eig <- eigen(A)
  eig$vectors %*% sweep(t(eig$vectors), 1, g(eig$values), FUN = "*")
}

# function to compute scatters and return orthogonal vectors as well as eigenvalues
scatter2 <- function(xw, Lloc, laplacian){
  n <- nrow(xw); p <- ncol(xw)
  s2 <- (1 / n) * crossprod(xw, Lloc) %*% xw 
  eig <- eigen(s2, symmetric = TRUE)
  if(laplacian){ inds <- p:1 }else{ inds <- 1:p } #
  list(u = eig$vectors[,inds], d = eig$values[inds])
}


# graph weighted bss 
# if k_mat is provided then Mittinen scatter is used - weights are overwritten 
grabss <- function(x,  
                  weights_adj,
                  nbr_nns = 2,
                  local_cen = TRUE,
                  local_cov = TRUE, 
                  laplacian = TRUE,
                  local_laplacian = TRUE,
                  gammam = c("rho","unit"),
                  gfun = NULL){
  # init
  n <- nrow(x)
  p <- ncol(x)
  eps <- .Machine$double.eps
  x_cen <- matrix(0, n, p)
  
  # compute distances
  graph.ob <- graph_from_adjacency_matrix(weights_adj, weighted = TRUE)
  dists <- distances(graph.ob, mode = "all", algorithm = "dijkstra")
  dists[is.infinite(dists)] <- 0
  rho_mat <- dists / max(dists)
  #sig     <- stats::quantile(rho_mat, probs = 0.9)
  sig     <- mad(rho_mat)#
  rho_mat <- exp(- rho_mat^2 / sig^2 * 1 / 2)
  
  # names
  varnames <- if(!is.null(colnames(x))){ colnames(x) }else{ as.character(1:p) }

  # center data 
  if(local_cen){
    cen <- lapply(seq_len(n),function(idx){ 
                  wt  <- rho_mat[idx, ]
                  xwc <- colSums(sweep(x, 1, wt / sum(wt) , "*"))
                  xwc })
    cen <- do.call("rbind", cen)  
  } else if(!local_cen ){
    cen <- t(matrix(colMeans(x), ncol = n, nrow = p))      # global
  }
  
  # when local cov FALSE then use global scale to scale 
  if(local_cov){
    sc <- lapply(seq_len(n),function(idx){ 
                  wt  <- rho_mat[idx, ]
                  xwc <- colSums(sweep((x - cen)^2, 1, wt / sum(wt) , "*"))
                  xwc })
    sc <- do.call("rbind", sc)  
  }else{
    sc  <- t(matrix(apply(x, 2, sd), ncol = n, nrow = p))  # global
  }
  
  
  # center data
  x_cen <- lapply(seq_len(n), function(idx){ scale(x, center = cen[idx,], scale = sc[idx,])  })
  
  ### compute S1s
  if(local_cov){
    S1s <- lapply(seq_len(n), function(idx){ S1(x_cen[[idx]], rho_mat[idx, ], eps) })
  }else{
    cen_gl <- glb_center(x)
    s1  <- S1(x - cen_gl, rep(1,n), eps)
    S1s <- lapply(seq_len(n), function(idx){ s1 })
  }

    
  
  # sqrt inverse of first scatter
  S1_invsq <- lapply(seq_len(n), function(idx){ gpower(S1s[[idx]], - 1 / 2)  })
  
    # loadings gwpca
  S1ortho <- lapply(seq_len(n), function(idx){ t(S1_invsq[[idx]]$orthmat) })
  
  ## whiten data 
  x_wh <- lapply(seq_len(n), function(idx) x_cen[[idx]] %*% S1_invsq[[idx]]$mat)
  
  # gamma matrix - how far we go from s
  if(gammam == "rho"){ 
    gamma_mat <- as(rho_mat, "sparseMatrix")
  }else if(gammam == "unit"){
    gamma_mat <- as((1 * (dists <= nbr_nns)), "sparseMatrix")
  }

  
  S2_diago <- 
  if(local_laplacian){
           if(!is.null(gfun)){ message("g can only be used with sptial_laplacian = FALSE")}
           lapply(seq_len(n), function(idx){ 
                      #W_loc <- gamma_mat * rho_mat[idx,]
                      W_loc <-  as(outer(rho_mat[idx,], rho_mat[idx,], FUN = "*") * gamma_mat,"sparseMatrix")
                      A_loc <- if(laplacian){ diag(rowSums(W_loc)) - W_loc }else{ W_loc }
                      scatter2(x_wh[[idx]], A_loc, laplacian)
                  })  
      }else{
        W_gl <- as(rho_mat * gamma_mat, "sparseMatrix")
        A_gl <- if(laplacian){ diag(rowSums(W_gl)) - W_gl }else{ W_gl }
        A_gl <- if(!is.null(gfun)){ gA(A_gl, function(x){gfun(x)})}else{ A_gl }
        lapply(seq_len(n), function(idx){ scatter2(x_wh[[idx]], A_gl, laplacian) })
  }
  
  
  
  # remove object
  rm(x_wh)

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
  
  # unscaled unmixing 
  unmixing_local_unscaled <- 
    lapply(seq_len(n), function(idx){ sweep(unmixing_local[[idx]], 2, 1 / sc[idx,], "*")  }) # back to original scale

  # global unmixing
  unmixing_global <-  1 / n * Reduce("+", unmixing_local)
  colnames(unmixing_global) <- varnames
  
  # unscaled uglobal unmixing 
  unmixing_global_unscaled <-  1 / n * Reduce("+", unmixing_local_unscaled)
  colnames(unmixing_global_unscaled) <- varnames
  
  
  # eigvalues
  d <- lapply(seq_len(n), function(idx) as.vector(S2_diago[[idx]]$d))
  d <- do.call("rbind", d)
  dcummulated <- t(apply(d, 1, function(u){ cumsum(u)/ sum(u) }))
  
  # return winning loading variable name
  winningloading  <- lapply(1:n, function(idx){ rowMaxAbs(unmixing_local[[idx]], varnames) })
  winningloading <- do.call("rbind", winningloading)
    
  orderedloadings <- lapply(1:n, function(idx){ rowOrderAbs(unmixing_local[[idx]], varnames) })
  orderedloadings <- abind(orderedloadings,along = 3)
  orderedloadings <- aperm(orderedloadings,perm = c(3,1,2))
  
  # unscaled
  winningloading_unscaled  <- lapply(1:n, function(idx){ rowMaxAbs(unmixing_local_unscaled[[idx]], varnames) })
  winningloading_unscaled <- do.call("rbind", winningloading_unscaled)
    
  orderedloadings_unscaled <- lapply(1:n, function(idx){ rowOrderAbs(unmixing_local_unscaled[[idx]], varnames) })
  orderedloadings_unscaled <- abind(orderedloadings_unscaled,along = 3)
  orderedloadings_unscaled <- aperm(orderedloadings_unscaled,perm = c(3,1,2))
  

  
  # garbage collector
  gc()
  
  return(list(scores_ss  = scores_ss, 
              scores_edges = scores_edges,
              scores_global = scores_global,
              unmixing_local = unmixing_local, 
              unmixing_global = unmixing_global,
              S1ortho = S1ortho,
              unmixing_local_unscaled = unmixing_local_unscaled,
              unmixing_global_unscaled = unmixing_global_unscaled,
              eigenvalues = d, 
              eigenvaluescum = dcummulated,
              centers = cen, 
              scales = sc,
              winningloading = winningloading,
              orderedloadings = orderedloadings,
              winningloading_unscaled = winningloading_unscaled,
              orderedloadings_unscaled = orderedloadings_unscaled,
              dists = dists,
              gamma_mat = gamma_mat))
}


