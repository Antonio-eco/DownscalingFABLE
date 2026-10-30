# =============================================================================
# Spatial GHG (biomass) emissions/sequestration by land-cover class,
# CUMULATIVE across all modelled periods (2020-2050)
#
# Inputs:
#   - rasterized_layer.tif  : SpatRaster where each cell's value ("ns_int") is
#                              an integer ID that links the cell to rows in
#                              the long-format table below.
#   - downscaled_LUC.rds    : long table with one row per (cell, year, origin
#                              class "lu.from", destination class "lu.to"),
#                              including a "GHG_biomass" column (positive =
#                              emission, negative = sequestration, in the
#                              table's native units per 5-year period).
#   - br_states.shp / br_biomes.shp : boundaries, any CRS (reprojected here).
#
# Output:
#   - One grid figure with one panel PER CLASS, showing the SUM of
#     GHG_biomass across all periods (2020-2050) per cell, on a diverging
#     green (net sequestration) <-> red (net emission) scale - matching the
#     reference figure's colour convention. State borders (thin, light grey)
#     and biome borders (thick, black) are overlaid, with lon/lat axes.
# =============================================================================

library(terra)
library(dplyr)
library(tidyr)

# ---- 0. User settings ------------------------------------------------------
tif_path    <- "C:/Users/User/Desktop/DownscalingFABLE/Output/BRA/261027_BRA_CT_HILDA_rasterized_layer.tif"
rds_path    <- "C:/Users/User/Desktop/DownscalingFABLE/Output/BRA/261027_BRA_CT_HILDA_downscaled_LUC.rds"
states_path <- "C:/Users/User/Desktop/FABLE/DownscalingFABLE/br_states.shp"
biomes_path <- "C:/Users/User/Desktop/FABLE/DownscalingFABLE/br_biomes.shp"

classes_keep <- c("cropland", "forest", "otherland", "pasture", "urban")

out_png <- "ghg_emissions_cumulative_by_class.png"

legend_title <- "GHG_biomass"   # relabel to match your units, e.g. "Mt CO2"

# ---- 1. Load data -----------------------------------------------------------
luc <- readRDS(rds_path)
r   <- rast(tif_path)

# ---- 2. Aggregate GHG_biomass by cell, destination class & period ----------
# (as specified: replace NA with 0, then sum across all origins (lu.from)
# landing in each lu.to class, per cell & period)
inputs <- luc %>%
  dplyr::mutate(GHG_biomass = tidyr::replace_na(GHG_biomass, 0)) %>%
  dplyr::group_by(ns, lu.to, times) %>%
  dplyr::summarise(GHG_biomass = sum(GHG_biomass), .groups = "drop")

# ---- 3. Link the table's cell ID ("ns") to the raster's cell ID ("ns_int") --
ns_sorted <- sort(unique(luc$ns))
id_tab <- data.frame(ns = ns_sorted, ns_int = seq_along(ns_sorted) + 1L)

check_vals <- values(r)[, 1]
match_rate <- mean(check_vals[!is.na(check_vals)] %in% id_tab$ns_int)
message(sprintf("ns <-> ns_int match rate: %.1f%%", 100 * match_rate))
stopifnot(match_rate > 0.99)

inputs <- inputs %>% mutate(ns_int = id_tab$ns_int[match(ns, id_tab$ns)])

classes <- intersect(classes_keep, unique(inputs$lu.to))

# ---- 4. Sum GHG_biomass across all periods, per cell & class ----------------
cum_tbl <- inputs %>%
  group_by(ns_int, lu.to) %>%
  summarise(GHG_biomass = sum(GHG_biomass, na.rm = TRUE), .groups = "drop")

# ---- 5. Build one cumulative raster per class --------------------------------
change_rasters <- lapply(classes, function(cls) {
  sub <- cum_tbl %>% filter(lu.to == cls)
  rc <- r
  rc[] <- sub$GHG_biomass[match(values(r)[, 1], sub$ns_int)]
  names(rc) <- cls
  rc
})
names(change_rasters) <- classes
change_stack <- rast(change_rasters)

# ---- 5b. Add a "Total (all classes)" layer: per-cell sum across classes ----
change_stack <- c(change_stack, setNames(sum(change_stack), "Total (all classes)"))
plot_layers <- names(change_stack)

# ---- 6. Boundaries; reproject everything to lon/lat for display -------------
states <- vect(states_path)
biomes <- vect(biomes_path)
if (!same.crs(states, "EPSG:4326")) states <- project(states, "EPSG:4326")
if (!same.crs(biomes, "EPSG:4326")) biomes <- project(biomes, "EPSG:4326")

change_stack <- project(change_stack, "EPSG:4326")

add_boundaries <- function() {
  plot(states, border = "lightgray", lwd = 0.3, add = TRUE)
  plot(biomes, border = "black", lwd = 0.8, add = TRUE)
}

# ---- 7. Plot: one cumulative panel per class, in a grid ---------------------
# NOTE: `type = "continuous"` is required - terra auto-switches to a
# categorical legend when a layer has very few distinct values, which would
# otherwise silently ignore the diverging palette/range and misrepresent
# the map.
#
# Each panel uses its OWN symmetric colour scale (based on that layer's max
# |value|), rather than one scale shared across all classes. Cropland's
# emissions are an order of magnitude larger than other classes' small
# sequestration values, so a single shared scale would wash the latter out
# to near-white. If you'd rather compare classes on one common scale, use a
# single `max_abs <- max(abs(values(change_stack)), na.rm = TRUE)` and pass
# that as every panel's `range=` instead of the per-layer `cls_max` below.
pal_symmetric <- colorRampPalette(c("blue", "white", "firebrick"))(100)

# Build a diverging green/white/red palette for an ASYMMETRIC range, with
# white forced exactly on 0 (sized proportionally to each side of zero).
# Used for the "Total" panel, whose positive (emission) values are much
# larger in magnitude than its negative (sequestration) values - a plain
# symmetric scale would compress the negative side to near-invisibility.
make_asym_pal <- function(vmin, vmax, n = 100) {
  n_neg <- max(1, round(n * (0 - vmin) / (vmax - vmin)))
  n_pos <- n - n_neg
  c(
    colorRampPalette(c("blue", "white"))(n_neg),
    colorRampPalette(c("white", "firebrick"))(n_pos + 1)[-1]
  )
}

n_panels  <- length(plot_layers)
ncol_grid <- 3
nrow_grid <- ceiling(n_panels / ncol_grid)

png(out_png, width = 700 * ncol_grid, height = 480 * nrow_grid, res = 150)
par(mfrow = c(nrow_grid, ncol_grid), mar = c(3, 3, 3, 4), oma = c(0, 0, 2, 0))
for (lyr in plot_layers) {
  vals <- values(change_stack[[lyr]])
  if (lyr == "Total (all classes)") {
    # Most of the country has a small, widespread NEGATIVE value (net
    # sequestration) dwarfed by a few large positive (emission) outlier
    # cells. A plain linear colour scale - even an asymmetric one sized to
    # the true min/max - would still render that widespread small signal as
    # essentially white. Instead we colour by a SIGNED-LOG transform
    # (sign(x) * log1p(|x|)), which compresses the extremes and expands
    # resolution near zero, then relabel the legend in the original units.
    rc <- app(change_stack[[lyr]], function(x) sign(x) * log1p(abs(x)))
    vmin_t <- min(values(rc), na.rm = TRUE)
    vmax_t <- max(values(rc), na.rm = TRUE)
    pal_use   <- make_asym_pal(vmin_t, vmax_t)
    range_use <- c(vmin_t, vmax_t)

    orig_min <- min(vals, na.rm = TRUE)
    orig_max <- max(vals, na.rm = TRUE)
    breaks_orig <- sort(unique(c(orig_min, -5, -1, -0.1, 0, 0.1, 1, 5, orig_max)))
    breaks_orig <- breaks_orig[breaks_orig >= orig_min & breaks_orig <= orig_max]
    at_t <- sign(breaks_orig) * log1p(abs(breaks_orig))

    plot(rc, col = pal_use, range = range_use, type = "continuous",
         main = lyr, axes = TRUE,
         plg = list(title = legend_title, cex = 0.8,
                     at = at_t, labels = formatC(breaks_orig, format = "g", digits = 2)))
  } else {
    lyr_max <- max(abs(vals), na.rm = TRUE)
    if (lyr_max <= 0) lyr_max <- 1   # placeholder scale for all-zero layers
    plot(change_stack[[lyr]], col = pal_symmetric, range = c(-lyr_max, lyr_max),
         type = "continuous", main = lyr, axes = TRUE,
         plg = list(title = legend_title, cex = 0.8))
  }
  add_boundaries()
}
mtext(
  "Cumulative GHG_biomass by class, 2020-2050 (blue = sequestration, red = emission)",
  side = 3, outer = TRUE, font = 2, cex = 1.1
)
dev.off()
message("Saved: ", out_png)
