# imports
library('SpatialBSS')
library('StatDA')
library('ggmap')
library('sf')
library("GWmodel")
source("utils.R")
source("gwbss.R")
library('abind')

# data 
data("moss")

coords <- as.matrix(moss[, c('XCOO', 'YCOO')])
row.names(coords) <- NULL

coords_sf <- sf::st_as_sf(data.frame(coords), 
                          coords = c(1,2))
st_crs(coords_sf) <- '+proj=utm +zone=35 +ellps=GRS80 +datum=NAD83'
coords_sf <- st_transform(coords_sf, '+proj=longlat +datum=WGS84')
coords_ll <- st_coordinates(coords_sf)

kola_map <- get_stamenmap(bbox = unname(st_bbox(coords_sf)), 
                          maptype = 'terrain-background', 
                          zoom = 7,
                          crop = TRUE)

vars <- c("Ag","Al","As","B","Ba","Bi","Ca","Cd","Co","Cr","Cu","Fe","Hg","K","Mg","Mn","Mo",
          "Na","Ni","P","Pb","Rb","S","Sb","Si","Sr","Th","Tl","U","V","Zn")

elems <- as.matrix(moss[, vars])
row.names(elems) <- NULL

c_mat <- contras_mat(elems)

elems_ilr <- log(elems) %*% c_mat

# kernel definition
sbss_res <- SpatialBSS:::sbss.default(x = elems_ilr, coords = coords, 
                                      kernel_type = "ring", kernel_parameters = c(0, 50000))



d1s <- SpatialPointsDataFrame(coords, as.data.frame(scale(elems_ilr, center = TRUE, scale = FALSE)))
bw <- GWmodel::bw.gwpca(d1s, vars=colnames(d1s@data), kernel = "gaussian")
gwpca_res <- gwpca(d1s, vars=colnames(d1s@data), kernel = "gaussian", bw = 50000, scores = TRUE)
pca_res <- princomp(d1s@data, )

gwbss_res_1 <- gwbss(x = elems_ilr, coords = coords, bw = 100000, scatter = "1")
gwbss_res_2 <- gwbss(x = elems_ilr, coords = coords, bw = 100000, scatter = "2")

#gwbss_res_2 <- gwbss2(x =  elems_ilr, coords = coords, bw = 100000, kernel_type = "ball", kernel_parameters = bw,
#                            spatial_type = list(spatial_mean = TRUE, spatial_scores = TRUE))


# scores plots
ic_idx <- 1
plot_map(coords_ll, sbss_res$s[, ic_idx], map = kola_map, quant = TRUE)
plot_map(coords_ll, gwbss_res_1$s[, ic_idx], map = kola_map, quant = TRUE)
plot_map(coords_ll, gwbss_res_2$s[, ic_idx], map = kola_map, quant = TRUE)
a <- do.call(rbind, 
             lapply(1:length(gwpca_res$gwpca.scores), function(idx) gwpca_res$gwpca.scores[[idx]][idx,]))
plot_map(coords_ll, a[, ic_idx], map = kola_map, quant = TRUE)
plot_map(coords_ll, -pca_res$scores[, ic_idx], map = kola_map, quant = TRUE)


plot_map(coords_ll, gwbss_res_1$d[, ic_idx], map = kola_map, quant = FALSE)
plot_map(coords_ll, gwbss_res_2$d[, ic_idx], map = kola_map, quant = FALSE)

ic_idx <- 2
plot_map(coords_ll, sbss_res$s[, ic_idx], map = kola_map, quant = TRUE)
plot_map(coords_ll, gwbss_res_1$s[, ic_idx], map = kola_map, quant = TRUE)
plot_map(coords_ll, gwbss_res_2$s[, ic_idx], map = kola_map, quant = TRUE)


# winnig loadings
ic_idx <- 1
plot_map(coords_ll, apply(gwbss_res_1$loadings,1,function(A){ which.max(A[,ic_idx]) }), map = kola_map, quant = FALSE)
plot_map(coords_ll, apply(gwbss_res_2$loadings,1,function(A){ which.max(A[,ic_idx]) }), map = kola_map, quant = FALSE)
