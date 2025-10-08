# library(sf)
# library(here)
# 
# population <- st_read(here("Data", "global", "Population_2020_Zonal_Statistics.geojson"),  quiet = FALSE)
# 
# # Keep only Polygon and MultiPolygon
# population_poly <- population[st_geometry_type(population) %in% c("POLYGON", "MULTIPOLYGON"), ]
# 
# # Keep only Points
# population_points <- population[st_geometry_type(population) == "POINT", ]
# 
# # Extract polygons from GeometryCollection
# population_geom_collections <- population[st_geometry_type(population) == "GEOMETRYCOLLECTION", ]
# population_fixed_geom <- st_collection_extract(population_geom_collections, "POLYGON")
# 
# # Combine all geometries into one sf object
# population_clean <- rbind(population_poly, population_points, population_fixed_geom)

#-----------------------------------------------------------------
# from chatgpt
#-----------------------------------------------------------------

# Packages
library(sf)
library(dplyr)
library(purrr)
library(fs)
library(here)

# ---- 1. Input & output folders ----
in_dir  <- here("Data", "global")                 # your current folder
out_dir <- path(in_dir, "cleaned")                # new subfolder for cleaned files
dir_create(out_dir)                               # create if missing

# All GeoJSON files in the folder
in_files <- dir_ls(in_dir, regexp = "\\.geojson$", type = "file")

# ---- Cleaning function for one file ----
clean_one <- function(path) {
  message("Reading: ", path)
  x <- st_read(path, quiet = TRUE, stringsAsFactors = FALSE)
  
  # Drop Z/M if present
  x <- st_zm(x, drop = TRUE, what = "ZM")
  
  # Add filename (handy for joins/debug)
  x <- mutate(x, .file = path_file(path))
  
  # Geometry types per feature
  gtypes <- st_geometry_type(x, by_geometry = TRUE)
  
  # Keep polygon-like, make valid (avoids downstream issues)
  poly_like <- x[gtypes %in% c("POLYGON", "MULTIPOLYGON"), , drop = FALSE]
  if (nrow(poly_like)) poly_like <- st_make_valid(poly_like)
  
  # Keep point-like (POINT and MULTIPOINT)
  pts <- x[gtypes %in% c("POINT", "MULTIPOINT"), , drop = FALSE]
  
  # Extract POLYGON & MULTIPOLYGON from GeometryCollections
  gc <- x[gtypes == "GEOMETRYCOLLECTION", , drop = FALSE]
  gc_poly <- if (nrow(gc)) {
    suppressWarnings({
      # Extract polygons
      poly_part <- st_collection_extract(gc, "POLYGON")
      # Extract linestrings or points if you want them separately, but not needed here
      
      # Convert polygons to multipolygons for consistency
      poly_part <- st_cast(poly_part, "MULTIPOLYGON", warn = FALSE)
    })
    poly_part
  } else {
    gc[0, ]
  }
  
  
  # Combine; result is sfc_GEOMETRY (mixed) which is fine for many workflows
  out <- bind_rows(poly_like, pts, gc_poly)
  
  # Drop empties
  out <- out[!st_is_empty(out), , drop = FALSE]
  
  # Optional: standardize attribute-geometry relationship (silences some warnings)
  st_set_agr(out, "constant")
  
  # Output path in "cleaned" subfolder
  out_path <- path(out_dir, path_file(path))
  message("Writing cleaned file → ", out_path)
  
  st_write(out, out_path, delete_dsn = TRUE, quiet = TRUE)
  
  return(out_path)
}

# ---- 4. Run on all files ----
out_paths <- map_chr(in_files, clean_one)

message("Done. Cleaned ", length(out_paths), " files written to: ", out_dir)
