library(sf)
population <- st_read(here("Data", "global", "Population_2020_Zonal_Statistics.geojson"),  quiet = FALSE)

# Keep only Polygon and MultiPolygon
population_poly <- population[st_geometry_type(population) %in% c("POLYGON", "MULTIPOLYGON"), ]

# Keep only Points
population_points <- population[st_geometry_type(population) == "POINT", ]

# Extract polygons from GeometryCollection
population_geom_collections <- population[st_geometry_type(population) == "GEOMETRYCOLLECTION", ]
population_fixed_geom <- st_collection_extract(population_geom_collections, "POLYGON")

# Combine all geometries into one sf object
population_clean <- rbind(population_poly, population_points, population_fixed_geom)


