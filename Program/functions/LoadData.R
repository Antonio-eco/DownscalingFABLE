
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


## All mapping data is stored in one excel file names "mapping_code"
LoadMappingData <- function(filename, datafiles){
  
  mappingdata_list <- list()
  
  for(cur in seq_len(nrow(datafiles))){
    
    temp <-  readxl::read_excel(here("Data", paste0(filename, ".xlsx")),
                                sheet = datafiles[cur, 2])
    
    mappingdata_list[[ datafiles[cur, 1] ]] <- temp
  }
  
  return(mappingdata_list)
}
