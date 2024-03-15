contras_mat <- function(x) {
  c_mat <- matrix(0, nrow = ncol(x), ncol = ncol(x) - 1)
  for (i in 1:ncol(c_mat)) {
    c_mat[1:i, i] <- 1/i
    c_mat[i + 1, i] <- (-1)
    c_mat[, i] <- c_mat[, i] * sqrt(i/(i + 1))
  }
  return(c_mat)
}

clr_func <- function(x) {
  x_clr <- apply(x, 1, function(row) row / prod(row) ^ (1 / length(row)))
  x_clr <- log(x_clr)
  return(t(x_clr))
}



plot_map <- function(coords, variable, map, quant = TRUE, gradient_col = TRUE, title = "") {
  
  # get rid of NAs
  indna  <- is.na(variable)
  coords <- coords[!indna,]
  variable <- variable[!indna]
  
  # make data
  x <- data.frame(x = coords[, 1], y = coords[, 2], variable = variable) 
  #map <- openproj(map)#, projection = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
 

  
  if (quant) {
    q <- quantile(variable, c(0, 0.05, 0.25, 0.75, 0.95, 1))
    leg.round = 2 
    leg.wid = 4
    leg.just = "right"
    leg = rep(NA, length(q) - 1)
    leg[1] = paste("  ", roundpretty(q[1], leg.round), 
                   "-", format(roundpretty(q[1 + 1], leg.round), 
                               width = leg.wid, justify = leg.just))
    for (i in 2:(length(q) - 1)) {
      leg[i] = paste(">", roundpretty(q[i], leg.round), 
                     "-", format(roundpretty(q[i + 1], leg.round), 
                                 width = leg.wid, justify = leg.just))
    }
    
    x$variable <- cut(variable, 
                      breaks = q,
                      include.lowest = TRUE)
    
    shapes <- rev(c(17, 17, 20, 15, 15))
    shape_size <- rev(c(1.7, 1.0, 0.3, 1.2, 1.8))
    shape_size <- 1.6 * rev(c(1.3, 0.6, 0.2, 0.9, 1.5))
    color <- c('royalblue3', 'royalblue1', 'black', 'brown1', 'brown3')
    
    g <-  autoplot.OpenStreetMap(map) + 
      geom_point(data = x, aes(x = x, y = y, color = variable, shape = variable),size = 1.3) + 
      theme_bw() +
      xlab("") +
      ylab("") +
      theme(legend.position = "bottom") +
      scale_shape_manual(values = shapes,
                         labels = leg) + 
      scale_size_manual(values = shape_size,
                        labels = leg) + 
      scale_color_manual(values = color,
                         labels = leg) + 
      guides(shape = guide_legend(nrow = 3, byrow = TRUE)) +
      theme(legend.title = element_blank()) +
      theme(legend.text = element_text(size = 12))+
      theme(axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()) 
  } else {
    g <-  autoplot.OpenStreetMap(map) + 
      geom_point(data = x, aes(x = x, y = y, color = variable, alpha = 1)) + 
      theme_bw() +
      xlab("") +
      ylab("") +
      theme(legend.position = "bottom") +
    theme(legend.title = element_blank()) +
      theme(legend.text = element_text(angle = 45)) +
      theme(legend.text = element_text(size = 11)) +
      theme(axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()) 
    g <- if(gradient_col){ g + scale_color_gradientn(values = scales::rescale(c(-max(abs(variable)),0,max(abs(variable)))),
                                                     limits=c(-max(abs(variable)),max(abs(variable))),
                                                     colours = terrain.colors(10))}else{ g }
  }
  g <- g + ggtitle(title)
  #plot(g)
  return(g)
}


plot_load <- function(coords, val, val_qu, nbr_coords, map, title = ""){
  
# make dataframe  
qs <- quantile(abs(val_qu), seq(0,1,0.25),  na.rm = TRUE)
vcut      <- cut(abs(val_qu), breaks = qs, include.lowest = TRUE, ordered_result = TRUE)
dim(vcut) <- c(nrow(val_qu),ncol(val_qu))
levcut <- levels(vcut)
vcutm <-  t(apply(abs(val_qu), 1, function(u){  cut(u, breaks = qs, include.lowest = TRUE) }))

colnames(vcutm) <- colnames(val)
x   <- data.frame(x = coords[,1], y = coords[,2], val); colnames(x) <- c("x","y", colnames(val))
x   <- reshape2::melt(x, id.vars = c("x","y"))
xqu <- data.frame(x = coords[,1], y = coords[,2], vcutm); colnames(xqu) <- c("x","y", colnames(val))
xqu <- reshape2::melt(xqu, id.vars = c("x","y"))

x <- merge(x, xqu, by = c("x","y","variable"))
colnames(x) <- c("x","y","variable","value","valuequ")
x[,"valuequ"] <-  factor(x[,"valuequ"], levels = levcut, ordered = TRUE)

# make plot 
g <-  ggplot2::autoplot(map)+
      geom_point(data = x, aes(x = x, y = y, shape = factor(value),  color = factor(value), size = factor(valuequ)),
                 stroke=0.7) + 
      #scale_shape_manual(values= seq.int(nbr_coords)) +
      scale_shape_manual(values=rep(c(0:2,5:6,9:10,11:12,14), times=4)) + 
      scale_size_manual(values= seq(1.5,3.5, length.out = 4)) + 
      theme_bw() +
      facet_wrap(~variable, nrow = 2)+
      xlab("") +
      ylab("") +
      theme(legend.position = "bottom") +
      theme(legend.title = element_blank()) +
      theme(legend.text = element_text(angle = 45)) +
      theme(legend.text = element_text(size = 11)) + 
      theme(axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()) +
      guides(shape = guide_legend(nrow = 2, override.aes = list(size = 3)))+
      ggtitle(title)
  g
}


