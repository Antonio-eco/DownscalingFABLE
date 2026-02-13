### generer gif land use change % per pixel per land cover type

make_ns_canvas_50km <- function(grid_sf, ns_map, crs_equal_area = 6933, res_m = 50000) {
  grid_eq <- sf::st_transform(grid_sf, crs_equal_area)
  
  # attach ns_int to grid polygons
  grid_eq$ns_int <- ns_map$ns_int[match(grid_eq$id_c, ns_map$id_c)]
  
  r_template <- terra::rast(
    terra::ext(grid_eq),
    resolution = res_m,
    crs = paste0("EPSG:", crs_equal_area)
  )
  
  r_ns <- terra::rasterize(
    terra::vect(grid_eq),
    r_template,
    field = "ns_int",
    touches = TRUE
  )
  
  plot_df <- terra::as.data.frame(r_ns, xy = TRUE, na.rm = FALSE)
  names(plot_df)[3] <- "ns"
  plot_df %>% dplyr::filter(!is.na(ns))
}


compute_stepwise_pct_change <- function(to_plot, grid_table, multiplier_to_ha = 1000) {
  # grid_table must have: ns, area (in ha)
  pixel_area <- grid_table %>%
    dplyr::transmute(ns = ns, pixel_area_ha = area)
  
  inputs <- to_plot %>%
    dplyr::group_by(ns, lu.to, times) %>%
    dplyr::summarise(
      value_ha = multiplier_to_ha * sum(value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(times = as.integer(as.character(times))) %>%
    dplyr::arrange(ns, lu.to, times) %>%
    dplyr::group_by(ns, lu.to) %>%
    dplyr::mutate(
      times_from = dplyr::lag(times),
      delta_ha   = value_ha - dplyr::lag(value_ha)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::filter(!is.na(delta_ha)) %>%
    dplyr::left_join(pixel_area, by = "ns") %>%
    dplyr::mutate(value = 100 * delta_ha / pixel_area_ha) %>%
    dplyr::select(ns, lu.to, times_from, times, value)
  
  inputs
}


global_diverging_limits <- function(x) {
  m <- max(abs(x), na.rm = TRUE)
  c(-m, m)
}


plot_one_change_map <- function(canvas_df, inputs, lu, yr, limits,
                                title = NULL,
                                low = "red", mid = "white", high = "blue") {
  
  inputs_sub <- inputs %>%
    dplyr::filter(lu.to == lu, times == yr) %>%
    dplyr::select(ns, value)
  
  plot_df <- canvas_df %>%
    dplyr::left_join(inputs_sub, by = "ns")   # keeps NA where no value
  
  ggplot2::ggplot(plot_df) +
    ggplot2::geom_raster(ggplot2::aes(x = x, y = y, fill = value), alpha = 0.9) +
    ggplot2::scale_fill_gradient2(
      low = low, mid = mid, high = high,
      limits = limits,
      oob = scales::squish,
      na.value = "lightgrey",
      name = "% of pixel (5-year change)"
    ) +
    ggplot2::coord_equal() +
    ggthemes::theme_map() +
    ggplot2::labs(title = title %||% paste0(lu, " – ", yr)) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.key.width = grid::unit(1.1, "cm"),
      panel.background = ggplot2::element_rect(fill = "white")
    )
}


# 0) Get model output table
to_plot <- results_FABLE_Chaturvedi$out.res   # expects ns, lu.to, times, value

# 1) Build canvas once (50 km pixels)
canvas_df <- make_ns_canvas_50km(
  grid_sf = grid,        # <- your sf grid (50km polygons)
  ns_map  = ns_map,
  crs_equal_area = 6933,
  res_m = 50000
)

# 2) Compute all stepwise % pixel changes
# ensure numeric times
inputs <- compute_stepwise_pct_change(
  to_plot = to_plot,
  grid_table = grid_table,          # must have ns + area (ha)
  multiplier_to_ha = 1000           # set to 1 if res$out.res$value already in ha
)

inputs <- inputs %>%
  dplyr::mutate(times = as.integer(as.character(times)))

combos <- inputs %>%
  dplyr::distinct(lu.to, times) %>%
  dplyr::arrange(lu.to, times)

plots <- purrr::pmap(
  list(lu = combos$lu.to, yr = combos$times),
  function(lu, yr) {
    
    p <- plot_one_change_map(
      canvas_df = canvas_df,
      inputs    = inputs,
      lu        = lu,
      yr        = yr,
      limits    = limits,
      title     = paste0("5-year change (% pixel): ", lu, " (", yr - 5, "–", yr, ")")
    )
    
    ggplot2::ggsave(
      filename = here("Output", country, paste0("change_", lu, "_", yr, ".jpeg")),
      plot = p, device = "jpeg",
      width = 6, height = 6, units = "in", dpi = 300
    )
    
    p
  }
)
