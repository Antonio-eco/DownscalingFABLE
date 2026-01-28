make_color_scale <- function(limits,
                             label,
                             low = "blue",
                             high = "yellow",
                             na_value = "lightgrey") {
  ggplot2::scale_fill_gradient(
    low = low,
    high = high,
    limits = limits,
    oob = scales::squish,
    na.value = na_value,
    name = label
  )
}


LUC_plot_restoration_percent <- function(res, rasterfile, grid,
                                         year = NULL, LU = NULL, limits,
                                         label = "Area in % of pixel") {
  
  # rasterfile is a SpatRaster with values = ns (ns_int)
  plot_df <- terra::as.data.frame(rasterfile, xy = TRUE, na.rm = FALSE) %>% 
    filter(!is.na(ns_int))
  
  # value column name can vary (e.g., "lyr.1"); standardize to "ns"
  names(plot_df)[3] <- "ns"
  
  to.plot <- res$out.res
  
  # build inputs at the level you want to plot
  inputs <- to.plot %>%
    dplyr::group_by(ns, lu.to, times) %>%
    dplyr::summarise(value = sum(value), .groups = "drop")
  
  # if (!is.null(year)) inputs <- dplyr::filter(inputs, times == year)
  # if (!is.null(LU))   inputs <- dplyr::filter(inputs, lu.to == LU)
  
  # join: many pixels per ns is OK (it paints the polygon)
  plot_df <- dplyr::left_join(plot_df, inputs, by = "ns")
  
  # optional: keep NAs as 0 so background cells don't go grey
  # (choose what you prefer)
  plot_df <- dplyr::mutate(plot_df, value = tidyr::replace_na(value, 0)) %>% 
    mutate(value = ifelse(value == 0, NA, value))
  
  grid <- sf::st_transform(grid, 6933)
  
  plot_obj <- ggplot2::ggplot(plot_df) +
    ggplot2::geom_raster(ggplot2::aes(x = x, y = y, fill = value), alpha = 0.8) +
    make_color_scale(global_limits, label = label) +
    
    # 🔲 GRID BORDERS (GENERAL)
    ggplot2::geom_sf(
      data = grid,
      fill = NA,
      color = "grey",
      linewidth = 0.1
    ) +
    
    ggplot2::coord_sf(expand = FALSE) +
  
    ggthemes::theme_map() +
    ggplot2::theme(
      legend.position = "bottom",
      legend.key.width = grid::unit(1.1, "cm"),
      legend.title = element_text(size = 12),
      panel.background = ggplot2::element_rect(fill = "white")
    )
  
  list(LUC.plot = plot_obj, plot.df = plot_df)
}
