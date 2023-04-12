# imports
library('SpatialBSS')
library('StatDA')
library('ggmap')
library('sf')
library("GWmodel")
source("utils.R")
source("gwbss.R")
library('abind')
library("gridExtra")

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

# apply gwbss
bw <- 100000

gwbss_res_1 <- gwbss(x = elems_ilr, coords = coords, bw = bw, 
                     spatial_mean = TRUE, 
                     S2_type = "scov",
                     field_order = "gwica")
gwbss_res_2 <- gwbss(x = elems_ilr, coords = coords, bw = bw, 
                     spatial_mean = FALSE, 
                     S2_type = "vario",
                     field_order = "gwica")
gwbss_res_3 <- gwbss(x = elems_ilr, coords = coords, bw = bw, 
                     spatial_mean = TRUE, 
                     S2_type = "sbssw",
                     field_order = "gwica")
gwbss_res_4 <- gwbss(x = elems_ilr, coords = coords, bw = bw, 
                     spatial_mean = TRUE, 
                     S2_type = "sfobi",
                     field_order = "gwica")

# apply sbss
sbss_res <- SpatialBSS:::sbss.default(x = elems_ilr, coords = coords, 
                                      kernel_type = "ring", kernel_parameters = c(0, bw))

# plots
idx <- 2
g_scov <- plot_map(coords_ll, gwbss_res_1$s[, idx], map = kola_map, quant = TRUE, title = "scov")
g_vario <- plot_map(coords_ll, gwbss_res_2$s[, idx], map = kola_map, quant = TRUE, title = "vario")
g_sbssw <- plot_map(coords_ll, gwbss_res_3$s[, idx], map = kola_map, quant = TRUE, title = "sbssw")
g_sfobi <- plot_map(coords_ll, gwbss_res_4$s[, idx], map = kola_map, quant = TRUE, title = "sfobi")
grid.arrange(g_scov, g_vario, g_sbssw, g_sfobi, nrow = 2)




d1s <- SpatialPointsDataFrame(coords, as.data.frame(scale(elems_ilr, center = TRUE, scale = FALSE)))
bw <- GWmodel::bw.gwpca(d1s, vars=colnames(d1s@data), kernel = "gaussian")
gwpca_res <- gwpca(d1s, vars=colnames(d1s@data), kernel = "gaussian", bw = 50000, scores = TRUE)
pca_res <- princomp(d1s@data, )


gwbss_res_cm <- gwbss(x = elems_ilr, coords = coords, bw = 100000, scatter = "2")

gwbss_res_cr <- gwbss2(x = elems_ilr, coords = coords, bw = 100000,
                           spatial_type = list(spatial_mean = TRUE, spatial_scores = TRUE))


# scores plots
pcm1 = plot_map(coords_ll, gwbss_res_cm$s[, 1], map = kola_map, quant = FALSE, title = "cm1")
pcr1 = plot_map(coords_ll, gwbss_res_cr$s[, 1], map = kola_map, quant = FALSE, title = "cr1")

pcm2 = plot_map(coords_ll, gwbss_res_cm$s[, 2], map = kola_map, quant = FALSE, title = "cm2")
pcr2 = plot_map(coords_ll, gwbss_res_cr$s[, 2], map = kola_map, quant = FALSE, title = "cr2")

grid.arrange(pcm1,pcm2,pcr1,pcr2,nrow =2)

ic_idx <- 1
plot_map(coords_ll, sbss_res$s[, ic_idx], map = kola_map, quant = TRUE)
plot_map(coords_ll, gwbss_res_1$s[, ic_idx], map = kola_map, quant = TRUE)
plot_map(coords_ll, gwbss_res_2$s[, ic_idx], map = kola_map, quant = TRUE)
a <- do.call(rbind, 
             lapply(1:length(gwpca_res$gwpca.scores), function(idx) gwpca_res$gwpca.scores[[idx]][idx,]))
plot_map(coords_ll, a[, ic_idx], map = kola_map, quant = TRUE)
plot_map(coords_ll, -pca_res$scores[, ic_idx], map = kola_map, quant = TRUE)


plot_map(coords_ll, gwbss_res_2$d[, 1], map = kola_map, quant = FALSE)
plot_map(coords_ll, gwbss_res_2$d[, 2], map = kola_map, quant = FALSE)
plot_map(coords_ll, gwbss_res_2$d[, 3], map = kola_map, quant = FALSE)
plot_map(coords_ll, gwbss_res_2$d[, 4], map = kola_map, quant = FALSE)

ic_idx <- 2
plot_map(coords_ll, sbss_res$s[, ic_idx], map = kola_map, quant = TRUE)
plot_map(coords_ll, gwbss_res_1$s[, ic_idx], map = kola_map, quant = TRUE)
plot_map(coords_ll, gwbss_res_2$s[, ic_idx], map = kola_map, quant = TRUE)


# winnig loadings
ic_idx <- 1
plot_map(coords_ll, apply(gwbss_res_1$loadings,1,function(A){ which.max(A[,ic_idx]) }), map = kola_map, quant = FALSE)
plot_map(coords_ll, apply(gwbss_res_2$loadings,1,function(A){ which.max(A[,ic_idx]) }), map = kola_map, quant = FALSE)
