### Maps for restoration documentation (FOLUR)
### Country-adaptive, reproducible, and comparable colour scales across maps
###
### What this script does (high level):
### 1) Loads the *downscaling model objects* produced earlier (priors, design matrix, etc.).
### 2) Runs downscalR once (only to get the standard "results" structure used by plotting tools).
### 3) Loads whichever restoration datasets exist for the chosen country (Bastin/Williams/Fesenmyer/Chaturvedi…).
### 4) Converts each restoration dataset to the same (ns, times, lu.from, lu.to, value) schema.
### 5) Computes one common colour scale limit across the final values to ensure maps are comparable.
### 6) Plots and saves all available maps as JPEG.

# ---- Parameters ----
country <- "IND"
stamp   <- "260129"   # identifies the run that produced country_start_areas / X_long / pred_coeff_long
today   <- format(Sys.Date(), "%y%m%d")
tag     <- paste0(stamp, "_", country, "_CurrentTrends_HILDA")

# Threshold below which values are hidden (set to NA to not display)
# (This mimics your previous “<5% should not appear on maps” rule.)
threshold <- 5

# EPSG:6933 is the equal-area CRS used by your 50km grid workflow
crs_equal_area <- 6933

# ---- Libraries ----
library(here)
library(dplyr)
library(tidyr)
library(readr)
library(purrr)
library(ggplot2)
library(geojsonio)
library(stringr)
library(terra)
library(readxl)
library(sf)        # needed because we use geom_sf + st_transform
library(ggthemes)  # used in theme_map
library(scales)    # used for oob = squish in colour scale

if (requireNamespace("conflicted", quietly = TRUE)) {
  conflicted::conflict_prefer("filter",    "dplyr", quiet = TRUE)
  conflicted::conflict_prefer("select",    "dplyr", quiet = TRUE)
  conflicted::conflict_prefer("lag",       "dplyr", quiet = TRUE)
  conflicted::conflict_prefer("rename",    "dplyr", quiet = TRUE)
  conflicted::conflict_prefer("left_join", "dplyr", quiet = TRUE)
}

# ---- Sources ----
source(here("Program", "wrangle_plot.R")) # wrangle_plot.R contains function: format_to_id_map
source(here("Program", "LUC_plot_restoration_percent.R")) 

# ---- Helpers ----
# -----------------------------------------------------------------------------
# safe_rest_df()
# Why:
# - Not all countries have the same restoration datasets.
# - Instead of failing when a file is missing, we return NULL so the script
#   gracefully skips that dataset.
#
# How:
# - path_vec = c("Data", country, "<file.geojson>")
# - loader   = a function that reads the file and returns a data.frame
# -----------------------------------------------------------------------------
safe_rest_df <- function(path_vec, loader) {
  path <- do.call(here, as.list(path_vec))
  if (!file.exists(path)) return(NULL)
  loader(path)
}

# -----------------------------------------------------------------------------
# add_spec()
# Why:
# - We store “what to plot” in a single list called restoration_specs.
# - Each entry contains:
#   df    : the restoration data (id_c + p)
#   label : legend label
#   file  : output filename
# This makes plotting scalable: add one dataset -> add one spec.
# -----------------------------------------------------------------------------
add_spec <- function(specs, name, df, label, file) {
  if (is.null(df) || nrow(df) == 0) return(specs)
  specs[[name]] <- list(df = df, label = label, file = file)
  specs
}

# -----------------------------------------------------------------------------
# make_restoration_result()
# Why:
# - Your plotting functions expect the downscalR-like schema:
#   ns, times, lu.from, lu.to, value
# - Most restoration datasets come as per-grid-cell percentages p with id_c.
# - We also apply your “hide values below threshold” rule here (single source of truth).
#
# Notes:
# - We set values < threshold to NA so they are not painted.
# -----------------------------------------------------------------------------
make_restoration_result <- function(base_results, rest_df, ns_map,
                                    times = 2050,
                                    lu_to = "newforest",
                                    threshold = 5) {
  
  stopifnot(all(c("id_c", "p") %in% names(rest_df)))
  
  out <- rest_df %>%
    dplyr::mutate(id_c = as.character(id_c)) %>%
    dplyr::left_join(
      ns_map %>%
        dplyr::mutate(id_c = as.character(id_c)) %>%
        dplyr::select(id_c, ns_int),
      by = "id_c"
    ) %>%
    dplyr::transmute(
      ns      = ns_int,
      times   = times,
      lu.from = "any",
      lu.to   = lu_to,
      value   = dplyr::if_else(p < threshold, NA_real_, as.numeric(p))
    )
  
  base_results$out.res <- out
  base_results
}

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
# save_plot_jpeg()
# Why:
# - Centralizes export settings (size, dpi) so all outputs are consistent.
# -----------------------------------------------------------------------------
save_plot_jpeg <- function(plot_obj, country, filename, width = 5, height = 6, dpi = 300) {
  ggplot2::ggsave(
    filename = here("Output", country, filename),
    plot     = plot_obj,
    device   = "jpeg",
    width    = width,
    height   = height,
    units    = "in",
    dpi      = dpi
  )
}

# =============================================================================
# Load core objects from the chosen run
# =============================================================================

country_start_areas <- readr::read_rds(here("Output", country, paste0(stamp, "_country_start_areas.rds")))
X_long              <- readr::read_rds(here("Output", country, paste0(stamp, "_X_long.rds")))
pred_coeff_long     <- readr::read_rds(here("Output", country, paste0(stamp, "_pred_coeff_long.rds")))

# IMPORTANT: terra rasters should be read from file formats (GeoTIFF), not RDS.
rasterized_layer <- terra::rast(here("Output", country, paste0(tag, "_rasterized_layer.tif")))

start_map_reproj <- readr::read_rds(here("Output", country, paste0(tag, "_start_map_reproj.rds")))
start_map        <- readr::read_rds(here("Output", country, paste0(tag, "_start_map.rds")))
grid_sf          <- readr::read_rds(here("Output", country, paste0(tag, "_grid_sf.rds")))

ns_map <- readr::read_rds(here("Output", country, paste0(tag, "_ns_map.rds"))) %>%
  dplyr::mutate(id_c = as.character(id_c))

grid <- readr::read_csv(here("Data", "global", "grid50_equal_area.csv")) %>%
  dplyr::filter(iso3 == country) %>%
  dplyr::mutate(id_c = as.character(id_c))

FABLE_CT <- readr::read_rds(here("Output", country, paste0(tag, "_FABLE.rds")))

# =============================================================================
# Run downscale (only to obtain the standard results object structure)
# =============================================================================
# Note: we are not using downscaled outputs as “restoration maps”.
# We only need a downscalR-like container (res$out.res) to reuse plotting tools.

pathways <- list(
  FABLE_CT = FABLE_CT %>% dplyr::filter(times >= 2020)
)

results_by_pathway <- purrr::imap(pathways, function(targets, nm) {
  
  tmp <- downscalr::downscale(
    targets     = targets,
    start.areas = country_start_areas,
    xmat        = X_long,
    betas       = pred_coeff_long
  )
  
  tmp$out.res <- format_to_id_map(tmp, ns_map)
  message("Generated: results_", nm)
  tmp
})

results_FABLE_CT <- results_by_pathway$FABLE_CT

# =============================================================================
# Restoration maps: build country-adaptive “specs”
# =============================================================================
# A “spec” is the minimal info needed to produce one map:
# - df:   id_c + p (% in the cell)
# - label: legend title
# - file: output filename
#
# We only add a spec if its file exists for this country.

restoration_specs <- list()

# ---- Bastin (if file exists) ----
rest_bastin.df <- safe_rest_df(
  c("Data", country, "GlobalTreeRestorationPotential_Bastin.geojson"),
  loader = function(path) {
    geojsonio::geojson_read(path, what = "sp")@data %>%
      dplyr::filter(!is.na(mean)) %>%
      dplyr::transmute(id_c = as.character(id_c), p = mean)
  }
)

restoration_specs <- add_spec(
  restoration_specs,
  name  = "Bastin",
  df    = rest_bastin.df,
  label = "Bastin et al. (2019)\nPercentage of pixel prioritised\nfor forest restoration:",
  file  = paste0(today, "_Bastin_map.jpeg")
)

# ---- Williams (if file exists) ----
rest_williams.df <- safe_rest_df(
  c("Data", country, "PotentialForestRegenerationTropicalForest_Williams.geojson"),
  loader = function(path) {
    geojsonio::geojson_read(path, what = "sp")@data %>%
      dplyr::select(-id) %>%
      dplyr::mutate(id_c = as.character(id_c)) %>%
      dplyr::left_join(grid %>% dplyr::select(id_c, area), by = "id_c") %>%
      dplyr::mutate(p = pmin(100, 100 * sum / area)) %>%
      dplyr::select(id_c, p)
  }
)

restoration_specs <- add_spec(
  restoration_specs,
  name  = "Williams",
  df    = rest_williams.df,
  label = "Williams et al. (2024)\nPercentage of pixel prioritised for\nforest regen. in tropical regions:",
  file  = paste0(today, "_Williams_map.jpeg")
)

# ---- Fesenmyer (if file exists) ----
rest_fesenmyer.df <- safe_rest_df(
  c("Data", country, "ScenarioBasedTotalRestorationPotentialArea_Fesenmyer.geojson"),
  loader = function(path) {
    geojsonio::geojson_read(path, what = "sp")@data %>%
      dplyr::select(id_c, d0) %>%
      dplyr::mutate(id_c = as.character(id_c)) %>%
      dplyr::left_join(grid %>% dplyr::select(id_c, area), by = "id_c") %>%
      dplyr::mutate(p = pmin(100, 100 * d0 / area)) %>%
      dplyr::select(id_c, p)
  }
)

restoration_specs <- add_spec(
  restoration_specs,
  name  = "Fesenmyer",
  df    = rest_fesenmyer.df,
  label = "Fesenmyer et al (2025)\nPercentage of pixel prioritised\nfor forest restoration:",
  file  = paste0(today, "_Fesenmyer_map.jpeg")
)

# ---- Chaturvedi (if file exists) ----

map_chaturvedi <- readxl::read_excel(here("Data/mapping_code.xlsx"), sheet = "Chaturvedi")

rest_chaturvedi.df <- safe_rest_df(
  c("Data", country, "LandscapeRestorationOpportunities.geojson"),
  loader = function(path) {
    geojsonio::geojson_read(path, what = "sp")@data %>%
      dplyr::select(-id) %>%
      dplyr::mutate(id_c = as.character(id_c)) %>%
      tidyr::pivot_longer(-id_c) %>%
      dplyr::mutate(code = stringr::str_remove(name, "X")) %>%
      dplyr::left_join(map_chaturvedi %>% dplyr::mutate(code = as.character(code)),
                       by = "code") %>%
      dplyr::left_join(grid %>% dplyr::select(id_c, area), by = "id_c") %>%
      dplyr::mutate(p = pmin(100, 100 * value / area)) %>%
      dplyr::select(id_c, restoration_type, p)
  }
)

if (!is.null(rest_chaturvedi.df) && nrow(rest_chaturvedi.df) > 0) {
  
  restoration_specs <- add_spec(
    restoration_specs,
    name  = "Chaturvedi_wide",
    df    = rest_chaturvedi.df %>%
      dplyr::filter(restoration_type == "wide_scale_restoration") %>%
      dplyr::select(id_c, p),
    label = "Chaturvedi et al. (2018)\nPercentage of pixel prioritised\nfor wide scale restoration:",
    file  = paste0(today, "_Chaturvedi_wide_map.jpeg")
  )
  
  restoration_specs <- add_spec(
    restoration_specs,
    name  = "Chaturvedi_mosaic",
    df    = rest_chaturvedi.df %>%
      dplyr::filter(restoration_type == "mosaic_scale_restoration") %>%
      dplyr::select(id_c, p),
    label = "Chaturvedi et al. (2018)\nPercentage of pixel prioritised\nfor mosaic scale restoration:",
    file  = paste0(today, "_Chaturvedi_mosaic_map.jpeg")
  )
}

# Fail early if nothing to plot
if (length(restoration_specs) == 0) {
  stop("No restoration datasets found for this country (no matching geojson files in Data/<country>/).")
}


# =============================================================================
# Convert all restoration datasets into results objects + compute shared limits
# =============================================================================

# 1) Build all results objects first (this applies the threshold rule consistently)
rest_results <- purrr::imap(restoration_specs, function(spec, nm) {
  make_restoration_result(
    base_results = results_FABLE_CT,
    rest_df      = spec$df,
    ns_map       = ns_map,
    times        = 2050,
    threshold    = threshold
  )
})

# 2) Compute limits from what is actually plotted (after thresholding)
#    This ensures all maps share the same colour scale.
limits <- range(unlist(purrr::map(rest_results, ~ .x$out.res$value)), na.rm = TRUE)

# =============================================================================
# Plot + save
# =============================================================================

plots <- purrr::imap(restoration_specs, function(spec, nm) {
  
  res_obj <- rest_results[[nm]]
  
  p <- LUC_plot_restoration_percent(
    res        = res_obj,
    rasterfile = rasterized_layer,
    grid       = grid_sf,
    limits     = limits,
    label      = spec$label,
    crs_equal_area = crs_equal_area
  )
  
  save_plot_jpeg(p$LUC.plot, country, spec$file, width = 5, height = 6, dpi = 300)
  message("Saved: ", spec$file)
  
  p
})
