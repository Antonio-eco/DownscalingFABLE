# =============================================================================
# Land-use change by land-cover class, 2020-2050
#
# Inputs:
#   - rasterized_layer.tif  : SpatRaster where each cell's value ("ns_int") is
#                              an integer ID that links the cell to rows in
#                              the long-format table below.
#   - downscaled_LUC.rds    : long table with one row per (cell, year,
#                              origin class "lu.from", destination class
#                              "lu.to"), where "value" is the area (in
#                              1,000 ha) of that origin's land currently
#                              classified as "lu.to" in that year.
#   - br_states.shp         : Brazilian state boundaries (any CRS - reprojected
#                              on load).
#   - br_biomes.shp         : Brazilian biome boundaries (any CRS - reprojected
#                              on load).
#
# Output:
#   - One panel PER CLASS (net change, gain AND loss, 2020 -> 2050), on a
#     diverging red (loss) <-> green (gain) scale, all classes laid out in a
#     single grid figure. Each panel is overlaid with state borders (thin,
#     light grey) and biome borders (thick, black), and has lon/lat axes.
# =============================================================================

library(terra)
library(dplyr)
library(tidyr)

# ---- 0. User settings ------------------------------------------------------
tif_path    <- "C:/Users/User/Desktop/DownscalingFABLE/Output/BRA/261027_BRA_CT_HILDA_rasterized_layer.tif"
rds_path    <- "C:/Users/User/Desktop/DownscalingFABLE/Output/BRA/261027_BRA_CT_HILDA_downscaled_LUC.rds"
states_path <- "C:/Users/User/Desktop/FABLE/DownscalingFABLE/br_states.shp"
biomes_path <- "C:/Users/User/Desktop/FABLE/DownscalingFABLE/br_biomes.shp"

year_start <- "2020"
year_end   <- "2050"

out_png <- "luc_change_2020_2050.png"       # one grid figure, all classes

# ---- 1. Load data -----------------------------------------------------------
luc <- readRDS(rds_path)     # long table: times, ns, lu.to, lu.from, value, ...
r   <- rast(tif_path)        # raster of cell IDs ("ns_int")

# ---- 2. Link the table's cell ID ("ns") to the raster's cell ID ("ns_int") --
# "ns" (character cell code) and "ns_int" (integer stored in the raster) are
# linked through a simple rank order: the raster's integer codes are the
# sorted-rank of the unique "ns" values, offset by 1 (code 1 is reserved,
# e.g. for background/NA). The match-rate check below verifies this holds.
ns_sorted <- sort(unique(luc$ns))
id_tab <- data.frame(ns = ns_sorted, ns_int = seq_along(ns_sorted) + 1L)

check_vals <- values(r)[, 1]
match_rate <- mean(check_vals[!is.na(check_vals)] %in% id_tab$ns_int)
message(sprintf("ns <-> ns_int match rate: %.1f%%", 100 * match_rate))
stopifnot(match_rate > 0.99)

# ---- 3. Total area of each class, per cell & year ---------------------------
# For a given cell/year/destination class, the class's TOTAL area is the sum
# of area flowing in from every possible origin (lu.from).
area_tbl <- luc %>%
  mutate(ns_int = id_tab$ns_int[match(ns, id_tab$ns)]) %>%
  group_by(ns_int, times, lu.to) %>%
  summarise(area = sum(value, na.rm = TRUE), .groups = "drop")

# ---- 4. Net change between year_start and year_end, per cell & class --------
chg_tbl <- area_tbl %>%
  filter(times %in% c(year_start, year_end)) %>%
  pivot_wider(names_from = times, values_from = area,
              names_prefix = "y", values_fill = 0) %>%
  mutate(change = .data[[paste0("y", year_end)]] - .data[[paste0("y", year_start)]])

classes <- sort(unique(chg_tbl$lu.to))
classes <- setdiff(classes, "newforest")   # drop newforest (shown separately / not of interest here)

# ---- 5. Build one change-raster per class -----------------------------------
change_rasters <- lapply(classes, function(cls) {
  sub <- chg_tbl %>% filter(lu.to == cls)
  rc <- r
  rc[] <- sub$change[match(values(r)[, 1], sub$ns_int)]
  names(rc) <- cls
  rc
})
names(change_rasters) <- classes
change_stack <- rast(change_rasters)

# ---- 6. Load boundaries; reproject everything to lon/lat for display --------
# The raster is in an equal-area projected CRS (metres), which is good for
# area calculations but doesn't give readable lat/long axis ticks. For the
# maps we reproject the change rasters and the boundaries to EPSG:4326
# (WGS84 lon/lat), purely for plotting - the area math above was already
# done in the native CRS, so reprojecting afterwards doesn't bias the values.
states <- vect(states_path)
biomes <- vect(biomes_path)
if (!same.crs(states, "EPSG:4326")) states <- project(states, "EPSG:4326")
if (!same.crs(biomes, "EPSG:4326")) biomes <- project(biomes, "EPSG:4326")

change_stack <- project(change_stack, "EPSG:4326")

add_boundaries <- function() {
  plot(states, border = "lightgray", lwd = 0.3, add = TRUE)
  plot(biomes, border = "black", lwd = 0.8, add = TRUE)
}

# ---- 7. Plot: one diverging gain/loss panel per class, in a grid -----------
# Symmetric, continuous colour scale (via `range=`, not `breaks=`, so the
# legend stays a smooth bar) shared across all classes for comparability.
max_abs <- max(abs(values(change_stack)), na.rm = TRUE)
pal_div <- colorRampPalette(c("firebrick", "white", "blue"))(100)

n_classes <- length(classes)
ncol_grid <- 3
nrow_grid <- ceiling(n_classes / ncol_grid)

png(out_png, width = 700 * ncol_grid, height = 480 * nrow_grid, res = 150)
par(mfrow = c(nrow_grid, ncol_grid), mar = c(3, 3, 3, 4), oma = c(0, 0, 2, 0))
for (cls in classes) {
  # axes = TRUE draws the lon/lat tick marks and labels (map is in EPSG:4326)
  plot(change_stack[[cls]], col = pal_div, range = c(-max_abs, max_abs),
       main = cls, axes = TRUE, plg = list(title = "1,000 ha", cex = 0.8))
  add_boundaries()
}
mtext(
  paste0("Land-use change by class, ", year_start, "-", year_end, " (red = loss, blue = gain)"),
  side = 3, outer = TRUE, font = 2, cex = 1.2
)
dev.off()
message("Saved: ", out_png)
