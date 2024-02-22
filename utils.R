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

plot_map <- function(coords, variable, map, quant = TRUE, title = "") {
  
  x <- data.frame(x = coords[, 1], y = coords[, 2], variable = variable) 
  
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
    
    g <- ggmap::ggmap(map,
                      base_layer = ggplot(aes(x = x, y = y, color = variable, shape = variable, 
                                              size = variable), data = x)) + 
      geom_point(alpha=1) +
      theme_bw() +
      xlab(expression("Longitude"*~degree*W)) +
      ylab(expression("Latitude"*~degree*N)) +
      theme(legend.position = "bottom") +
      scale_shape_manual(values = shapes,
                         labels = leg) + 
      scale_size_manual(values = shape_size,
                        labels = leg) + 
      scale_color_manual(values = color,
                         labels = leg) + 
      guides(shape = guide_legend(nrow = 2, byrow = TRUE)) +
      theme(legend.title = element_blank()) +
      theme(legend.text = element_text(size = 11))
  } else {
    g <- ggmap::ggmap(map,
                      base_layer = ggplot(aes(x = x, y = y, color = variable), data = x)) + 
      geom_point(alpha=1, size = 3) +
      theme_bw() +
      xlab(expression("Longitude"*~degree*W)) +
      ylab(expression("Latitude"*~degree*N)) +
      theme(legend.position = "bottom") +
      #scale_color_gradientn(colours = terrain.colors(7)) +
      #scale_color_gradientn(colours = colorspace::diverge_hcl(7)) + 
      #scale_colo_continuous_diverging() +
      #scale_colour_gradientn(colours =rainbow(7))	+
      #scale_colour_gradientn(low = "grey", high = "brown") + 
      #scale_colour_gradient(low="red", high="blue")
    theme(legend.title = element_blank()) +
      theme(legend.text = element_text(angle = 45)) +
      theme(legend.text = element_text(size = 11))
  }
  g <- g + ggtitle(title)
  plot(g)
  return(g)
}
