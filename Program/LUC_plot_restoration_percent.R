# -----------------------------------------------------------------------------
# make_color_scale()
# Why:
# - We want all maps to share the same colour meaning.
# - That requires a shared min/max (limits) across maps.
# - oob = squish clips out-of-range values (rare but safe).
# -----------------------------------------------------------------------------
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

# -----------------------------------------------------------------------------
# LUC_plot_restoration_percent()
# Why:
# - The rasterized_layer is a SpatRaster where each pixel stores ns_int.
# - We join those pixels with per-cell values (value) to paint a map.
# - Grid borders are drawn from a vector grid (sf) on top for readability.
#
# Important design choice:
# - Raster is used only as a *canvas* for fast plotting (geom_raster).
# - Borders come from sf polygons; rasters don’t have “borders” as geometry.
# -----------------------------------------------------------------------------

LUC_plot_restoration_percent <- function(res, rasterfile, grid,
                                         year = NULL, LU = NULL, limits,
                                         label = "Area in % of pixel",
                                         crs_equal_area = 6933) {
  
  # 1) Convert raster to a dataframe of pixels: x,y + ns_int (as values)
  plot_df <- terra::as.data.frame(rasterfile, xy = TRUE, na.rm = FALSE)
  
  # The 3rd column name can vary ("lyr.1", etc.). Standardize to "ns".
  names(plot_df)[3] <- "ns"
  
  # Remove pixels that have no id
  plot_df <- plot_df %>% dplyr::filter(!is.na(ns))
  
  # 2) Prepare values to join (one value per ns; here, restoration maps are already per ns)
  to_plot <- res$out.res
  
  # Optional filters (kept for reuse)
  if (!is.null(year)) to_plot <- to_plot %>% dplyr::filter(times == year)
  if (!is.null(LU))   to_plot <- to_plot %>% dplyr::filter(lu.to == LU)
  
  inputs <- to_plot %>%
    dplyr::group_by(ns, lu.to, times) %>%
    dplyr::summarise(
      value = if (all(is.na(value))) NA_real_ else sum(value, na.rm = TRUE),
      .groups = "drop"
    )
  
  # 3) Join pixel canvas with values
  plot_df <- dplyr::left_join(plot_df, inputs, by = "ns")
  
  # 4) Ensure grid is in the same CRS as raster
  grid <- sf::st_transform(grid, crs_equal_area)
  
  # 5) Plot
  plot_obj <- ggplot2::ggplot(plot_df) +
    ggplot2::geom_raster(ggplot2::aes(x = x, y = y, fill = value), alpha = 0.8) +
    make_color_scale(limits = limits, label = label) +
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
      legend.title = ggplot2::element_text(size = 12),
      panel.background = ggplot2::element_rect(fill = "white")
    )
  
  list(LUC.plot = plot_obj, plot.df = plot_df)
}
