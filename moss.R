# imports
library('SpatialBSS')
library('StatDA')
library('ggmap')
library('sf')
library("GWmodel")
library('abind')
library("gridExtra")
library("osmdata")
library("cccd")
library("pheatmap")
library("JADE")
source("utils.R")
source("gwbss.R")

# data 
data("moss")
data("chorizon")

#x <- moss
x <- chorizon


coords <- as.matrix(x[, c('XCOO', 'YCOO')])
row.names(coords) <- NULL

coords_sf <- sf::st_as_sf(data.frame(coords), 
                          coords = c(1,2))
st_crs(coords_sf) <- '+proj=utm +zone=35 +ellps=GRS80 +datum=NAD83'
coords_sf <- st_transform(coords_sf, '+proj=longlat +datum=WGS84')
coords_ll <- st_coordinates(coords_sf)

#kola_map <- get_stamenmap(bbox = unname(st_bbox(coords_sf)),
#                         maptype = 'terrain-background',
 #                        zoom = 7,
#                        crop = TRUE)

#save(kola_map, file = "kolamap.RData")
load(file = "kolamap.RData")


# kola moss
#vars <- c("Ag","Al","As","B","Ba","Bi","Ca","Cd","Co","Cr","Cu","Fe","Hg","K","Mg","Mn","Mo",
#          "Na","Ni","P","Pb","Rb","S","Sb","Si","Sr","Th","Tl","U","V","Zn")

# chorizon
vars <- c("Al_XRF","Ca_XRF", "Fe_XRF","K_XRF","Mg_XRF","Mn_XRF","Na_XRF","P_XRF",
          "Si_XRF","Ti_XRF")

elems <- as.matrix(x[, vars])


#row.names(elems) <- NULL

# if in ilr wanted uncomment
# c_mat <- contras_mat(elems)
# elems_ilr <- log(elems) %*% c_mat


elems_log <- log(elems)
colnames(elems_log) <- vars
#elems_log  <- t(scale(t(elems_log), center = TRUE, scale = FALSE))
elems_log  <- scale(elems_log, center = TRUE, scale = TRUE)
#colnames(elems_log) <- 1:ncol(elems_log)

# apply gwbss
bw <- 100000
graph.obj <- nng(x = as.matrix(stats::dist(coords)), k = 40, mutual = FALSE)
weights_adj <- as.matrix(as_adjacency_matrix(graph.obj))
weights_adj <- weights_adj + t(weights_adj)
weights_adj[weights_adj!=0] = 1

# dist_mat <- as.matrix(dist(coords))
# dist_mat[dist_mat > bw] <- 0
# weights_adj <- matrix(0, nrow(dist_mat), nrow(dist_mat))
# weights_adj[dist_mat != 0] <- 1 #/ dist_mat[dist_mat != 0]
# weights_adj <- weights_adj / max(weights_adj)
  
#SpatialBSS::spatial_kernel_matrix(coords, kernel_type = "ball",
#                                  kernel_parameters = bw)[[1]]
# weights_adj <- GWmodel::gw.weight(as.matrix(dist(coords)), bw = bw, kernel = "gaussian", adaptive = FALSE)
# 
# weights_adj <- t(apply(weights_adj, 1, function(u){ 
#   u[u < quantile(u,0.75)] <-0
#   return(u)}))
# weights_adj = 0.5*(weights_adj + t(weights_adj))


gwbss_res_1 <- gwbss(x = elems_log, coords = coords, 
                     spatial_mean = TRUE, weights_adj = weights_adj, 
                     S2_type = "wcov",
                     field_order = "gwica",
                     kernel_type = "ball", kernel_parameters = bw)
gwbss_res_2 <- gwbss(x = elems_log, coords = coords, 
                     spatial_mean = TRUE, weights_adj = weights_adj, 
                     S2_type = "wvario",
                     field_order = "gwica",
                     kernel_type = "ball", kernel_parameters = bw)
gwbss_res_3 <- gwbss(x = elems_log, coords = coords, 
                     spatial_mean = TRUE, weights_adj = weights_adj, 
                     S2_type = "wgraph",
                     field_order = "gwica",
                     kernel_type = "ball", kernel_parameters = bw)
gwbss_res_4 <- gwbss(x = elems_log, coords = coords, 
                     spatial_mean = TRUE, weights_adj = weights_adj, 
                     S2_type = "sfobi",
                     field_order = "gwica",
                     kernel_type = "ball", kernel_parameters = bw)



  


# apply sbss
sbss_res <- SpatialBSS:::sbss.default(x = elems_log, coords = coords, 
                                      kernel_type = "ring", kernel_parameters = c(0, bw))

# scores plots
idx <- 1
g_scov <- plot_map(coords_ll, gwbss_res_1$s[, idx], map = kola_map, quant = TRUE, title = "wcov")
g_vario <- plot_map(coords_ll, gwbss_res_2$s[, idx], map = kola_map, quant = TRUE, title = "wvario")
g_sbssw <- plot_map(coords_ll, gwbss_res_3$s[, idx], map = kola_map, quant = TRUE, title = "wgraph")
g_sfobi <- plot_map(coords_ll, sbss_res$s[, idx], map = kola_map, quant = TRUE, title = "sbss")
#grid.arrange(g_scov, g_vario, g_sbssw, g_sfobi, nrow = 2)
p1 <- arrangeGrob(g_scov, g_vario, g_sbssw, g_sfobi, nrow = 2)
plot(p1)

idx <- 2
g_scov <- plot_map(coords_ll, gwbss_res_1$s[, idx], map = kola_map, quant = TRUE, title = "wcov")
g_vario <- plot_map(coords_ll, gwbss_res_2$s[, idx], map = kola_map, quant = TRUE, title = "wvario")
g_sbssw <- plot_map(coords_ll, gwbss_res_3$s[, idx], map = kola_map, quant = TRUE, title = "wgraph")
g_sfobi <- plot_map(coords_ll, sbss_res$s[, idx], map = kola_map, quant = TRUE, title = "sbss")
#grid.arrange(g_scov, g_vario, g_sbssw, g_sfobi, nrow = 2)
p2 <- arrangeGrob(g_scov, g_vario, g_sbssw, g_sfobi, nrow = 2)

idx <- 3
g_scov <- plot_map(coords_ll, gwbss_res_1$s[, idx], map = kola_map, quant = TRUE, title = "wcov")
g_vario <- plot_map(coords_ll, gwbss_res_2$s[, idx], map = kola_map, quant = TRUE, title = "wvario")
g_sbssw <- plot_map(coords_ll, gwbss_res_3$s[, idx], map = kola_map, quant = TRUE, title = "wgraph")
g_sfobi <- plot_map(coords_ll, sbss_res$s[, idx], map = kola_map, quant = TRUE, title = "sbss")
p3 <- arrangeGrob(g_scov, g_vario, g_sbssw, g_sfobi, nrow = 2)

idx <- 4
g_scov <- plot_map(coords_ll, gwbss_res_1$s[, idx], map = kola_map, quant = TRUE, title = "wcov")
g_vario <- plot_map(coords_ll, gwbss_res_2$s[, idx], map = kola_map, quant = TRUE, title = "wvario")
g_sbssw <- plot_map(coords_ll, gwbss_res_3$s[, idx], map = kola_map, quant = TRUE, title = "wgraph")
g_sfobi <- plot_map(coords_ll, sbss_res$s[, idx], map = kola_map, quant = TRUE, title = "sbss")
p4 <- arrangeGrob(g_scov, g_vario, g_sbssw, g_sfobi, nrow = 2)


pdf("scores.pdf", width = 8, height = 12)
plot(p1)
plot(p2)
plot(p3)
plot(p4)
dev.off()


# loadings plots
idx <- 1
g_scov  <- plot_map(coords_ll, gwbss_res_1$winningloading[,idx], map = kola_map, quant = FALSE, title = "wcov")
g_vario <- plot_map(coords_ll, gwbss_res_2$winningloading[,idx], map = kola_map, quant = FALSE, title = "wvario")
g_sbssw <- plot_map(coords_ll, gwbss_res_3$winningloading[,idx], map = kola_map, quant = FALSE, title = "wgraph")
#g_sfobi <- plot_map(coords_ll, sbss_res$s[, idx], map = kola_map, quant = FALSE, title = "sbss")
p1 <- arrangeGrob(g_scov, g_vario, g_sbssw, nrow = 2)
plot(p1)


# second biggest
nbrord <- 2
g_scov  <- plot_map(coords_ll, gwbss_res_1$orderedloadings[,idx,nbrord], map = kola_map, quant = FALSE, title = "wcov")
g_vario <- plot_map(coords_ll, gwbss_res_2$orderedloadings[,idx,nbrord], map = kola_map, quant = FALSE, title = "wcov")
g_sbssw <- plot_map(coords_ll, gwbss_res_3$orderedloadings[,idx,nbrord], map = kola_map, quant = FALSE, title = "wcov")
p1 <- arrangeGrob(g_scov, g_vario, g_sbssw, nrow = 2)
plot(p1)







## with graph 
#k_mat <- SpatialBSS::spatial_kernel_matrix(coords, kernel_type = "ball", kernel_parameters = bw)[[1]]
gwbss_res_1 <- graphbss(x = elems_log, 
                        dist_neighbors = 1,
                        k_mat = NULL, 
                        mtf = FALSE,
                        weights_adj = weights_adj, 
                        spatial_mean = FALSE, 
                        spatial_s1 = FALSE,
                        S2_type = "wcov",
                        geographic = list(flag = FALSE, coords = NULL, kernel_type = "ball", diam = NULL),
                        MDIs = FALSE)
gwbss_res_2 <- graphbss(x = elems_log, 
                        dist_neighbors = 1,
                        k_mat = NULL, 
                        mtf =FALSE,
                        weights_adj = weights_adj,
                        spatial_mean = FALSE,
                        spatial_s1 = FALSE,
                        S2_type = "wvario",
                        geographic = list(flag = FALSE, coords = NULL, kernel_type = "ball", diam = NULL),
                        MDIs = TRUE)
gwbss_res_3 <- graphbss(x = elems_log,
                        dist_neighbors = 1,
                        k_mat = NULL, 
                        mtf =FALSE,
                        weights_adj = weights_adj,
                        spatial_mean = FALSE, 
                        spatial_s1 = FALSE,
                        S2_type = "wgraph",
                        geographic = list(flag = FALSE, coords = NULL, kernel_type = "ball", diam = NULL),
                        MDIs = TRUE)
#gwbss_res_4 <- gwbss(x = elems_log, 
#                     spatial_mean = TRUE, weights_adj = weights_adj, 
#                     S2_type = "sfobi",
#                     field_order = "gwica",)

# apply sbss
sbss_res <- SpatialBSS:::sbss.default(x = elems_log, coords = coords, 
                                      kernel_type = "ring", kernel_parameters = c(0, bw))



# plot center deviations
g_scov <- plot_map(coords_ll, gwbss_res_1$dist_center, map = kola_map, quant = TRUE, title = "wcov")
g_vario <- plot_map(coords_ll, gwbss_res_2$dist_center, map = kola_map, quant = TRUE, title = "wvario")
g_sbssw <- plot_map(coords_ll, gwbss_res_3$dist_center, map = kola_map, quant = TRUE, title = "wgraph")
plot(arrangeGrob(g_scov, g_vario, g_sbssw, nrow = 2))

# plot scatter deviations
g_scov <- plot_map(coords_ll, gwbss_res_1$MDss1, map = kola_map, quant = FALSE, title = "wcov")
g_vario <- plot_map(coords_ll, gwbss_res_2$MDss1, map = kola_map, quant = FALSE, title = "wvario")
g_sbssw <- plot_map(coords_ll, gwbss_res_3$MDss1, map = kola_map, quant = FALSE, title = "wgraph")
plot(arrangeGrob(g_scov, g_vario, g_sbssw, nrow = 2))

# plot scores
idx <- 1
g_scov <- plot_map(coords_ll, gwbss_res_1$scores[,idx,1], map = kola_map, quant = TRUE, title = "wcov")
g_vario <- plot_map(coords_ll, gwbss_res_2$scores[,idx,1], map = kola_map, quant = TRUE, title = "wvario")
g_sbssw <- plot_map(coords_ll, gwbss_res_3$scores[,idx,1], map = kola_map, quant = TRUE, title = "wgraph")
g_sfobi <- plot_map(coords_ll, sbss_res$s[, idx], map = kola_map, quant = TRUE, title = "sbss")
#grid.arrange(g_scov, g_vario, g_sbssw, g_sfobi, nrow = 2)
p1 <- arrangeGrob(g_scov, g_vario, g_sbssw, g_sfobi, nrow = 2)



idx <- 2
g_scov <- plot_map(coords_ll, gwbss_res_1$scores[,idx,1], map = kola_map, quant = TRUE, title = "wcov")
g_vario <- plot_map(coords_ll, gwbss_res_2$scores[,idx,1], map = kola_map, quant = TRUE, title = "wvario")
g_sbssw <- plot_map(coords_ll, gwbss_res_3$scores[,idx,1], map = kola_map, quant = TRUE, title = "wgraph")
g_sfobi <- plot_map(coords_ll, sbss_res$s[, idx], map = kola_map, quant = TRUE, title = "sbss")
#grid.arrange(g_scov, g_vario, g_sbssw, g_sfobi, nrow = 2)
p2 <- arrangeGrob(g_scov, g_vario, g_sbssw, g_sfobi, nrow = 2)


idx <- 3
g_scov <- plot_map(coords_ll, gwbss_res_1$scores[,idx,1], map = kola_map, quant = TRUE, title = "wcov")
g_vario <- plot_map(coords_ll, gwbss_res_2$scores[,idx,1], map = kola_map, quant = TRUE, title = "wvario")
g_sbssw <- plot_map(coords_ll, gwbss_res_3$scores[,idx,1], map = kola_map, quant = TRUE, title = "wgraph")
g_sfobi <- plot_map(coords_ll, sbss_res$s[, idx], map = kola_map, quant = TRUE, title = "sbss")
p3 <- arrangeGrob(g_scov, g_vario, g_sbssw, g_sfobi, nrow = 2)


pdf("scores.pdf", width = 8, height = 12)
plot(p1)
plot(p2)
plot(p3)
dev.off()


# plot loadings
idx <- 1
g_scov  <- plot_map(coords_ll, gwbss_res_1$winningloading[,idx], map = kola_map, quant = FALSE, title = "wcov")
g_vario <- plot_map(coords_ll, gwbss_res_2$winningloading[,idx], map = kola_map, quant = FALSE, title = "wvario")
g_sbssw <- plot_map(coords_ll, gwbss_res_3$winningloading[,idx], map = kola_map, quant = FALSE, title = "wgraph")
#g_sfobi <- plot_map(coords_ll, sbss_res$s[, idx], map = kola_map, quant = FALSE, title = "sbss")
p1 <- arrangeGrob(g_scov, g_vario, g_sbssw, nrow = 2)


# second biggest
nbrord <- 2
g_scov  <- plot_map(coords_ll, gwbss_res_1$orderedloadings[,idx,nbrord], map = kola_map, quant = FALSE, title = "wcov")
g_vario <- plot_map(coords_ll, gwbss_res_2$orderedloadings[,idx,nbrord], map = kola_map, quant = FALSE, title = "wcov")
g_sbssw <- plot_map(coords_ll, gwbss_res_3$orderedloadings[,idx,nbrord], map = kola_map, quant = FALSE, title = "wcov")
p2 <- arrangeGrob(g_scov, g_vario, g_sbssw, nrow = 2)

# third biggest
nbrord <- 3
g_scov  <- plot_map(coords_ll, gwbss_res_1$orderedloadings[,idx,nbrord], map = kola_map, quant = FALSE, title = "wcov")
g_vario <- plot_map(coords_ll, gwbss_res_2$orderedloadings[,idx,nbrord], map = kola_map, quant = FALSE, title = "wcov")
g_sbssw <- plot_map(coords_ll, gwbss_res_3$orderedloadings[,idx,nbrord], map = kola_map, quant = FALSE, title = "wcov")
p3 <- arrangeGrob(g_scov, g_vario, g_sbssw, nrow = 2)


pdf("scores.pdf", width = 8, height = 12)
plot(p1)
plot(p2)
plot(p3)
dev.off()



# plot z(i,j) for first loading
idx <- 1
zij <- lapply(1:nrow(x), function(nidx){gwbss_res_1$scores_edges[[nidx]][,idx]})
zij <- do.call("rbind",zij)
zij[weights_adj == 0] = 0
pheatmap(zij,  cluster_rows = FALSE, cluster_cols = TRUE)
pheatmap(zij,  cluster_rows = FALSE, cluster_cols = FALSE)













######## with graph function

nbr_nns <- 2
gwbss_res_1 <- grabss(x = elems_log, 
                      weights_adj,
                      nbr_nns = nbr_nns,
                      spatial_mean = TRUE,
                      spatial_cov = TRUE, 
                      laplacian = TRUE)
gwbss_res_2 <- grabss(x = elems_log, 
                      weights_adj,
                      nbr_nns = nbr_nns,
                      spatial_mean = TRUE,
                      spatial_cov = TRUE, 
                      laplacian = FALSE)
gwbss_res_3 <- grabss(x = elems_log, 
                      weights_adj,
                      nbr_nns = nbr_nns,
                      spatial_mean = FALSE,
                      spatial_cov = FALSE, 
                      laplacian = TRUE)
gwbss_res_4 <- grabss(x = elems_log, 
                      weights_adj,
                      nbr_nns = nbr_nns,
                      spatial_mean = FALSE,
                      spatial_cov = FALSE, 
                      laplacian = FALSE)



# apply sbss
sbss_res <- SpatialBSS:::sbss.default(x = elems_log, coords = coords, 
                                      kernel_type = "ring", kernel_parameters = c(0, bw))



# plot center deviations
# g_scov <- plot_map(coords_ll, gwbss_res_1$dist_center, map = kola_map, quant = TRUE, title = "wcov")
# g_vario <- plot_map(coords_ll, gwbss_res_2$dist_center, map = kola_map, quant = TRUE, title = "wvario")
# plot(arrangeGrob(g_scov, g_vario, g_sbssw, nrow = 2))
# 
# # plot scatter deviations
# g_scov <- plot_map(coords_ll, gwbss_res_1$MDss1, map = kola_map, quant = FALSE, title = "wcov")
# g_vario <- plot_map(coords_ll, gwbss_res_2$MDss1, map = kola_map, quant = FALSE, title = "wvario")
# plot(arrangeGrob(g_scov, g_vario, g_sbssw, nrow = 2))

# plot scores
idx <- 1
g1 <- plot_map(coords_ll, gwbss_res_1$scores[,idx], map = kola_map, quant = TRUE, title = "L")
g2 <- plot_map(coords_ll, gwbss_res_2$scores[,idx], map = kola_map, quant = TRUE, title = "W")
g3 <- plot_map(coords_ll, gwbss_res_3$scores[,idx], map = kola_map, quant = TRUE, title = "L-const")
g4 <- plot_map(coords_ll, gwbss_res_4$scores[,idx], map = kola_map, quant = TRUE, title = "W-const")
g_sfobi <- plot_map(coords_ll, sbss_res$s[, idx], map = kola_map, quant = TRUE, title = "sbss")
plot(arrangeGrob(g1, g2,g3,g4, g_sfobi, ncol = 3))

plots = list()
for(i in 1:4){
  plots[[i]] <- plot_map(coords_ll, gwbss_res_1$scores_edges[,i,idx], map = kola_map, quant = TRUE, title = "L")
}
do.call("grid.arrange", c(plots, ncol=4))


idx <- 2
g1 <- plot_map(coords_ll, gwbss_res_1$scores[,idx], map = kola_map, quant = TRUE, title = "L")
g2 <- plot_map(coords_ll, gwbss_res_2$scores[,idx], map = kola_map, quant = TRUE, title = "W")
g3 <- plot_map(coords_ll, gwbss_res_3$scores[,idx], map = kola_map, quant = TRUE, title = "L-const")
g4 <- plot_map(coords_ll, gwbss_res_4$scores[,idx], map = kola_map, quant = TRUE, title = "W-const")
g_sfobi <- plot_map(coords_ll, sbss_res$s[, idx], map = kola_map, quant = TRUE, title = "sbss")
plot(arrangeGrob(g1, g2,g3,g4, g_sfobi, ncol = 3))




idx <- 3
g1 <- plot_map(coords_ll, gwbss_res_1$scores[,idx], map = kola_map, quant = TRUE, title = "L")
g2 <- plot_map(coords_ll, gwbss_res_2$scores[,idx], map = kola_map, quant = TRUE, title = "W")
g3 <- plot_map(coords_ll, gwbss_res_3$scores[,idx], map = kola_map, quant = TRUE, title = "L-const")
g4 <- plot_map(coords_ll, gwbss_res_4$scores[,idx], map = kola_map, quant = TRUE, title = "W-const")
g_sfobi <- plot_map(coords_ll, sbss_res$s[, idx], map = kola_map, quant = TRUE, title = "sbss")
plot(arrangeGrob(g1, g2,g3,g4, g_sfobi, ncol = 3))





# plot loadings
idx <- 1
g1  <- plot_map(coords_ll, gwbss_res_1$winningloading[,idx], map = kola_map, quant = FALSE, title = "L")
g2 <- plot_map(coords_ll, gwbss_res_2$winningloading[,idx], map = kola_map, quant = FALSE, title = "W")
g3 <- plot_map(coords_ll, gwbss_res_3$winningloading[,idx], map = kola_map, quant = FALSE, title = "L-const")
g4 <- plot_map(coords_ll, gwbss_res_4$winningloading[,idx], map = kola_map, quant = FALSE, title = "W-const")
plot(arrangeGrob(g1,g2,g3, g4, ncol = 2))


# second biggest
nbr_lo <- 1
idx <- 2
g1  <- plot_map(coords_ll, gwbss_res_1$orderedloadings[,idx,nbr_lo], map = kola_map, quant = FALSE, title = "L")
g2 <- plot_map(coords_ll, gwbss_res_2$orderedloadings[,idx,nbr_lo], map = kola_map, quant = FALSE, title = "W")
g3 <- plot_map(coords_ll, gwbss_res_3$orderedloadings[,idx,nbr_lo], map = kola_map, quant = FALSE, title = "L-const")
g4 <- plot_map(coords_ll, gwbss_res_4$orderedloadings[,idx,nbr_lo], map = kola_map, quant = FALSE, title = "W-const")
plot(arrangeGrob(g1,g2,g3, g4, ncol = 2))


# plot z(i,j) for first loading
idx <- 2
zij <- lapply(1:nrow(x), function(nidx){gwbss_res_1$scores_edges[,nidx,idx]})
zij <- do.call("rbind",zij)
zij[weights_adj == 0] = 0
pheatmap(zij,  cluster_rows = FALSE, cluster_cols = TRUE)
pheatmap(zij,  cluster_rows = FALSE, cluster_cols = FALSE)



# how it changes according to different g
# graph weighted bss 
# if k_mat is provided then Mittinen scatter is used - weights are overwritten 

gwbss_res_1 <- grabss(x = elems_log, 
                      weights_adj,
                      nbr_nns = nbr_nns,
                      spatial_mean = FALSE,
                      spatial_cov = FALSE, 
                      laplacian = TRUE,
                      spatial_laplacian = TRUE,
                      gfun = function(x){(x)^2})
plot_map(coords_ll, gwbss_res_1$scores_global[,2], map = kola_map, quant = TRUE, title = "L")
plot_map(coords_ll, gwbss_res_1$scores_edges[,20,1], map = kola_map, quant = TRUE, title = "L")
plot_map(coords_ll, gwbss_res_1$scores_edges[20,,1], map = kola_map, quant = TRUE, title = "L")


