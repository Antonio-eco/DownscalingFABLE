# =============================================================================
# K.  Plot -- otherland 2050, terra, reference-style breaks (1000 ha)
# =============================================================================

library(sf); library(terra); library(RColorBrewer)

geo_raw    <- readxl::read_excel("C:/Users/User/Desktop/FABLE/downscalr/pop_2020.xlsx") %>% select(id_c, .geo)
suffix_off <- 1e6
raw_id     <- geo_raw$id_c
is_suf     <- grepl("_a$", raw_id)
numeric_id <- suppressWarnings(
  ifelse(is_suf, as.numeric(sub("_a$","",raw_id))+suffix_off, as.numeric(raw_id)))

geoms     <- lapply(geo_raw$.geo, function(g) sf::st_read(g,quiet=TRUE)$geometry[[1]])
grid_sf   <- sf::st_sf(id_c=numeric_id, geometry=sf::st_sfc(geoms,crs=4326))
grid_vect <- terra::vect(grid_sf)
tmpl      <- terra::rast(terra::ext(grid_vect), resolution=0.05,
                         crs=terra::crs(grid_vect))
ns_raster <- terra::rasterize(grid_vect, tmpl, field="id_c")

out.res <- readRDS("C:/Users/User/Desktop/fable-scenario-comparison-brazil-main/data/luc/downscaled_LUC_UP51a_CT.rds")

# Stock (diagonal, lu.from == lu.to) = area that stayed in otherland by 2050.
# out.res$value is already in 1000 ha (kha) units, so no conversion is needed
# here -- values are used as-is for the reference plotting scale (kha).
otherland_2020 <- out.res %>%
  filter(lu.from == lu.to, lu.to == "otherland", times == "2020") %>%
  group_by(ns) %>%
  summarise(value_kha = sum(value, na.rm = TRUE), .groups = "drop") %>%
  mutate(id_c = suppressWarnings(
    ifelse(grepl("_a$", ns), as.numeric(sub("_a$","",ns))+suffix_off, as.numeric(ns))))

# value_kha is in 1000 ha; divide by 1000 to report the total in Mha (million ha)
total_Mha <- sum(otherland_2020$value_kha, na.rm = TRUE) / 1000
message(sprintf("otherland 2020 (CT) total = %.2f Mha", total_Mha))

# Reference breaks (1000 ha), scaled proportionally to this data's max
BREAKS_KHA <- c(0, 0.001, 5.9, 16.8, 37.3, 71.1, 130.4, 211.5, 260, 310)
HOTPINKS <- c("#FFFFFF","#fde0ef","#fcc5e0","#fa9fb5","#f768a1",
             "#dd3497","#ae017e","#7a0177","#49006a")

reclass_mat <- as.matrix(otherland_2020[, c("id_c","value_kha")])
r_cont      <- terra::classify(ns_raster, reclass_mat, others = NA)
mx          <- terra::global(r_cont, "max", na.rm = TRUE)[[1]]
breaks      <- BREAKS_KHA / 310 * mx
breaks[length(breaks)] <- mx * 1.0001

rcl_mat <- cbind(c(-Inf, breaks[-c(1,length(breaks))]), breaks[-1], seq_along(breaks[-1]))
r_class <- terra::classify(r_cont, rcl = rcl_mat, right = TRUE)

fmt <- function(x) if (x==0) "0" else if (x<1) formatC(x,format="f",digits=3) else formatC(x,format="f",digits=1)
LEG_LABELS <- paste0(sapply(BREAKS_KHA[-length(BREAKS_KHA)], fmt), " - ",
                     sapply(BREAKS_KHA[-1], fmt))

BRAZIL_EXT <- terra::ext(-74, -33.5, -35, 6)

# Optional shapefiles -- skip cleanly if not present (sandbox has none)
shp_states_path <- "C:/Users/User/Desktop/FABLE/DownscalingFABLE/br_states.shp"
shp_biomes_path <- "C:/Users/User/Desktop/FABLE/DownscalingFABLE/br_biomes.shp"
shp_states <- if (file.exists(shp_states_path)) terra::vect(shp_states_path) else NULL
shp_biomes <- if (file.exists(shp_biomes_path)) terra::vect(shp_biomes_path) else NULL

out_png <- "otherland_2020_mapbiomas_restricted.png"
png(out_png, width = 760, height = 820, res = 105, bg = "white")

terra::plot(
  terra::crop(r_class, BRAZIL_EXT),
  col    = HOTPINKS,
  type   = "classes",
  legend = FALSE,
  axes   = TRUE,
  mar    = c(3, 3, 3, 1),
  main   = sprintf("CT | otherland | 2020 | Total: %.2f Mha", total_Mha),
  cex.main = 0.95
)

if (!is.null(shp_states)) terra::plot(shp_states, border="lightgray", lwd=0.3, add=TRUE)
if (!is.null(shp_biomes)) terra::plot(shp_biomes, border="black",     lwd=0.8, add=TRUE)

usr <- par("usr")
legend(
  x = usr[2] - (usr[2]-usr[1])*0.02, y = usr[3] + (usr[4]-usr[3])*0.02,
  legend = LEG_LABELS, fill = HOTPINKS, border = NA,
  title = "1000 ha", title.adj = 0, bty = "n", cex = 0.72, xjust = 1, yjust = 0
)

dev.off()
message(sprintf("Saved: %s", out_png))

saveRDS(pasture_2020, "pasture_2020_ct_states_kha.rds")


scale_breaks <- c(0, 0.001, 5.9, 16.8, 37.3, 71.1, 130.4, 211.5, 260, 310)

up50 <- readRDS("C:/Users/User/Desktop/fable-scenario-comparison-brazil-main/data/luc/downscaled_LUC_UP50_ct.rds")
up70 <- readRDS("C:/Users/User/Desktop/FABLE/downscalr/downscaled_LUC_usocob_states.rds")
up51 <- readRDS("C:/Users/User/Desktop/DownscalingFABLE/Output/BRA/260924_BRA_NDC_ESACCI_downscaled_LUC.rds")

f50 <- up50[up50$times=="2020" & up50$lu.to=="Forest" & up50$lu.from=="Forest",]
f50_kha <- f50$value * 0.001   # legacy hectare->kha conversion
cat("UP50 (scaled to kha) bin distribution:\n")
print(table(cut(f50_kha, scale_breaks, include.lowest=TRUE)))

f70 <- up70[up70$times=="2020" & up70$lu.to=="Forest" & up70$lu.from=="Forest",]
f70_kha <- f70$value * 0.001   # legacy hectare->kha conversion
cat("UP70 (scaled to kha) bin distribution:\n")
print(table(cut(f70_kha, scale_breaks, include.lowest=TRUE)))

f51 <- up51[up51$times=="2020" & up51$lu.to=="forest",]
cat("\nUP51a (already ~kha) bin distribution:\n")
print(table(cut(f51$value, scale_breaks, include.lowest=TRUE)))

