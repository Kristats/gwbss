gwbss <- function(x, coords, bw, scatter = c("1", "2"), kernel_type = "ball", kernel_parameters = bw) {
  scatter <- match.arg(scatter)

  n <- nrow(coords)
  p <- ncol(x)

  weights <- GWmodel::gw.weight(as.matrix(dist(coords)), bw = bw, kernel = "gaussian", adaptive = FALSE)

  loadings <- array(0, dim = c(n, p, p))
  d <- matrix(0, nrow = n, ncol = p)
  scores <- matrix(0, nrow = n, ncol = p)
  gwmeans <- matrix(0, nrow = n, ncol = p)


  if (scatter == "2") {
    k_mat <- SpatialBSS::spatial_kernel_matrix(coords, kernel_type = kernel_type,
                                               kernel_parameters = kernel_parameters)
  }


  # cen <- lapply(seq_len(n),function(n.idx){
  #   wt  <- weights[n.idx, ]
  #   use <- abs(wt) > 0
  #   wt  <- wt[use]
  #   xwc <- colSums(sweep(x[use, ], 1, wt / sum(wt) , "*"))
  #   xwc
  # })
  # cen = do.call("rbind",cen)
  #

  # local bss
  for (idx in 1:n) {
    # prepare weights and data
    wt <- weights[idx, ]
    use <- wt > 0
    wt <- wt[use]
    dat <- x[use, ]

    # center data with gwmean
    gwmean <- colSums(sweep(dat, 1, wt / sum(wt) , "*")) #cen[idx,] #
    #dat_cen <- sweep(dat, 2, gwmean, "-")
    #dat_cen = dat - cen
    dat_cen_w <- sweep(dat_cen, 1, sqrt(wt), "*")

    # white data with gwcov
    s1 <- crossprod(dat_cen_w) / sum(wt)
    s1_inv_sqrt <- gpower(s1, - 1 / 2)

    # compute scatters
    if (scatter == "1") {
      dat_cen_w <- sweep(dat_cen, 1, wt, "*")
      s2 <- do.call(rbind,
                    lapply(x[idx, ] - gwmean, function(xi) colSums(xi * dat_cen_w)))
      s2 <- (s2 + t(s2)) / 2 / sum(wt)
      s2 <- list(s1_inv_sqrt$mat %*% s2 %*% s1_inv_sqrt$mat)
    } else if (scatter == "2") {
      k <- list(k_mat[[1]][use, use])
      s2 <- SpatialBSS::local_covariance_matrix(dat_cen_w, k, center = FALSE)
      s2 <- list(s1_inv_sqrt$mat %*% (s2[[1]] * n / sum(wt)) %*% s1_inv_sqrt$mat)
    }
    attr(s2, "lcov") <- "lcov"

    # diag scatter
    s_diag <- SpatialBSS:::diag_scatters(s2, ordered = FALSE)

    # scores, loadings, d
    gwmeans[idx, ] <- gwmean
    d[idx, ] <- diag(s_diag$d)
    l <- crossprod(s_diag$u, s1_inv_sqrt$mat)
    l <- l[, order(s1_inv_sqrt$vals, decreasing = TRUE)]
    l <- sweep(l, 1, sign(l[1, ]), "*")
    loadings[idx,, ] <- l
    scores[idx, ] <- (x[idx, ] - gwmean) %*% t(loadings[idx,, ])
  }

  return(list(s = scores, loadings = loadings, d = d, gwmeans = gwmeans))
}