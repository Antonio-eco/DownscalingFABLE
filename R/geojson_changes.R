library(sf)

lc2015 <- st_read("C:/Users/User/Desktop/DownscalingFABLE/Data/BRA/LandCoverESACCI2015.geojson")

hilda2015 <- st_read("C:/Users/User/Desktop/DownscalingFABLE/hilda_landcover_2015_brazil.geojson")
hilda_luc <- st_read("C:/Users/User/Desktop/DownscalingFABLE/hilda_landcover_change_2015_2020_brazil.geojson")

plot(hilda2015["X42"], main = "Forest (evergreen broad leaf) 2015")

library(ggplot2)

ggplot(hilda2015) +
  geom_sf(aes(fill = X55)) +   # substitute whichever column you want, e.g. X22, X33, X55...
  scale_fill_viridis_c(name = "Area (ha)") +   # or use your REDS/GREENS palette from earlier
  theme_minimal() +
  labs(title = "HILDA+ 2015 — Forest (deciduous, broad leaf)")

# Inspect it
head(lc2015)
names(lc2015)
class(lc2015)   # "sf" "data.frame" — behaves like a normal dataframe, plus a geometry column
nrow(lc2015)    # 3652 cells, matching what I see in the raw file

lcfable <- st_read("C:/Users/User/Desktop/FABLE/downscalr/grade_long_fable.geojson")
head(lcfable)
unique(lcfable$tema)
class(lcfable)   # "sf" "data.frame" — behaves like a normal dataframe, plus a geometry column
nrow(lcfable)    # 3652 cells, matching what I see in the raw file




#Land Cover Classification according to FABLE-C Brazil
tema_to_lu <- c(
  # --- Forest-type classes ---
  "forest"          = "X42",  
  "floodabforest"   = "X40",  
  "mangrove"        = "X40",  
  "savanna"         = "X42",  
  "woodsandbank"    = "X40",
  "wetland"         = "X42",
  "grassland"       = "X40",
  "othernonforest"  = "X40",
  "herbsandbank"    = "X42",
  
  
  # --- Other Land-type classes ---
  "rockyoutcrop"    = "X66",
  "beachsandbank"   = "X66",
  "hypersalinetf"   = "X66",
  "othernonveg"     = "X66",
  "mining"          = "X55",
  "photovoltaicpp"  = "X55",
  "forestplantatio" = "X55",
  
  # --- Pasture-type classes ---
  "pasture"         = "X33",
  
  # --- Cropland-type classes ---
  "soybean"         = "X22",
  "sugarcane"       = "X22",
  "rice"            = "X22",
  "othertempcrop"   = "X22",
  "cotton"          = "X22",
  "palmoil"         = "X23",
  "coffee"          = "X23",
  "citrus"          = "X23",
  "otherpercrop"    = "X24",
  "mosaicrural"     = "X24",
  
  # --- Urban-type classes ---
  "urban"          = "X11"
)

LU_CLASSES  <- c("X11", "X22", "X23", "X24", "X33", "X40", "X42", "X44", "X55",
                 "X66", "X77")

# Fix valor: comma decimal separator -> period, coerce to numeric
# valor arrives as character with "," as decimal (e.g. "30496,82327")
# gsub replaces ALL commas; as.numeric then parses correctly
lcfable <- lcfable %>%
  mutate(valor = as.numeric(gsub(",", ".", valor)))

n_na <- sum(is.na(lcfable$valor))
if (n_na > 0) warning(n_na, " NA values after valor conversion — check raw data")
message("valor conversion OK | range: [",
        round(min(lcfable$valor, na.rm = TRUE), 2), ", ",
        round(max(lcfable$valor, na.rm = TRUE), 2), "]")

# Sanity cap: valor cannot exceed the cell's total area (area column).
# Some temas (notably pivocentral) have erroneous values orders of magnitude
# larger than the cell area (e.g. 338 million ha for a 248k ha cell).
# Clamp these to the cell area so they don't corrupt start.areas or xmat.
n_over <- sum(lcfable$valor > lcfable$area, na.rm = TRUE)
if (n_over > 0) {
  message("Clamping ", n_over, " rows where valor > cell area (data errors in GeoJSON)")
  lcfable <- lcfable %>%
    mutate(valor = pmin(valor, area))
}

exclude <- c("nodata", "primary", "secondaryveg", "primaryloss",
"secondarygain", "secondaryloss", "carbsoilmedia", "desmatamento", "residuo",
"VegNativaFlorestal", "alto", "intermediario", "degradado", "carb_t_rec",
"protecaointegral", "terraindigena", "usosustentavel", "areasprotegidas",
"outrossist", "antrop", "aquaculture", "inundacao", "pivocentral")

# Filter to base year and mapped temas only
stock_long_2015 <- lcfable %>%
  filter(ano == 2015, tema %in% names(tema_to_lu)) %>%
  mutate(
    lu    = tema_to_lu[tema],
    valor = as.numeric(valor)
  ) %>%
  group_by(id_c, lu) %>%
  summarise(area_ha = sum(valor, na.rm = TRUE), .groups = "drop")

# Pivot to wide: one column per LU class (mirrors hilda_br_2015 structure)
mapbiomas_2015 <- stock_long_2015 %>%
  tidyr::pivot_wider(names_from = lu, values_from = area_ha, values_fill = 0) %>%
  # Ensure all 5 LU columns exist even if a class has no area anywhere
  { for (cl in LU_CLASSES) if (!cl %in% names(.)) .[[cl]] <- 0; . }

message(sprintf(
  "MapBiomas stock (%d): %d cells | classes: %s",
  2015, nrow(mapbiomas_2015),
  paste(LU_CLASSES, sapply(LU_CLASSES, function(cl) {
    sprintf("%.0f ha", sum(mapbiomas_2015[[cl]], na.rm = TRUE))
  }), sep = "=", collapse = ", ")
))


mapbiomas_2020 <- mapbiomas_2020 %>%
  mutate(across(c(3:13), ~ .x / 1000))

mapbiomas_2015 <- mapbiomas_2015 %>%
  mutate(across(c(3:13), ~ .x / 1000))

st_write(mapbiomas_2020, "LandCoverMapbiomas2020.geojson", driver = "GeoJSON")
st_write(mapbiomas_2015, "LandCoverMapbiomas2015.geojson", driver = "GeoJSON")

check <- st_read("LandCoverMapbiomas2015.geojson")
check


mapbiomas_2015 %>%
  filter(id_c == 84502)


