LUC_plot_restriction <- function(res, rasterfile, year=NULL, LU=NULL, color = "Greens", label = "Area in ha per pixel"){
  
  plot_df <- terra::as.data.frame(rasterfile, xy = TRUE, na.rm = FALSE)
  
  to.plot <- res$out.res
  
  colnames(plot_df) <- c("x", "y", "ns")
  
  ns = lu.to = times = value = x = y= NULL
  if(is.null(year) & is.null(LU)){
    inputs <- to.plot %>% dplyr::group_by(ns, lu.to, times) %>% dplyr::summarise(value = sum(value),.groups = "keep")
  } else if(!(is.null(year) | is.null(LU))){
    inputs <- to.plot %>% dplyr::group_by(ns, lu.to, times) %>% dplyr::summarise(value = sum(value),.groups = "keep") %>% subset(lu.to==LU & times==year)
  } else if(is.null(year)){
    inputs <- to.plot %>% dplyr::group_by(ns, lu.to, times) %>% dplyr::summarise(value = sum(value),.groups = "keep") %>% subset(lu.to==LU)
  } else {
    inputs <- to.plot %>% dplyr::group_by(ns, lu.to, times) %>% dplyr::summarise(value = sum(value),.groups = "keep") %>% subset(times==year)
  }
  
  plot_df <- merge(plot_df, inputs, by="ns")
  
  plot_obj <- ggplot2::ggplot() +
    geom_tile(data = plot_df, aes(x = x, y = y, fill = value), alpha = 0.8) +
    scale_fill_gradient2(low = "grey", mid = "white", high = "#1b7837",
                         midpoint = 0, limits = c(-1, 1), guide = "none") +
    coord_equal() + ggthemes::theme_map() +
    
    ggplot2::coord_equal() +
    ggthemes::theme_map() +
    ggplot2::theme(legend.position="bottom") +
    ggplot2::theme(legend.key.width=unit(2, "cm")) +
  
  # start a new fill scale just for the two-key legend
    ggnewscale::new_scale_fill() +
    geom_point(
      data = data.frame(key = factor(c("not priority","priority"),
                                     levels = c("not priority","priority"))),
      aes(x = Inf, y = Inf, fill = key), inherit.aes = FALSE, show.legend = TRUE
    ) +
    scale_fill_manual(
      name   = label,
      values = c("not priority" = "grey", "priority" = "#1b7837")
    ) +
    guides(fill = guide_legend(override.aes = list(shape = 22, size = 6))) +
    ggplot2::theme(legend.text = element_text(size = 14))+
    ggplot2::theme(legend.title = element_text(size = 16))
  
  return(list(LUC.plot=plot_obj, plot.df=plot_df))
}
