datafiles <- rbind.data.frame(
  c("population" , "Population2020.geojson"),
  c("protectedareas" , "ProtectedAreas.geojson"),
  c("travel", "TravelTime.geojson"),
  c("landcover","LandCoverHILDA2015.geojson"),
  c("landcover","LandCoverCopernicus2019.geojson"),
  c("landcoverchange", "LandCoverChangeHILDA2015_2019.geojson"),
  c("forestmanagement", "ForestManagement.geojson"),
  c("altitude", "Altitude.geojson"),
  c("slope", "Slope.geojson"),
  c("gaez", "GAEZCropDistribution2015.geojson"),
  c("livestock", "GLW4WorldGriddedLivestock2020.geojson"),
  c("ab_be_biomass", "AbovegroundBelowgroundBiomassCarbonDensity.geojson")
) %>% 
  setNames(c("rdataset", "filename"))


LoadSpatialData <- function(country, datafiles){
  
  spatialdata_list <- list()
  
  for(cur in seq_len(nrow(datafiles))){
    
    temp <- geojson_read(
      here("Data", country, datafiles[cur, 2]),
      what = "sp"
    )
    
    temp <- temp@data %>%
      select(-id)
    
    spatialdata_list[[ datafiles[cur, 1] ]] <- temp
  }
  
  return(spatialdata_list)
}