#install.packages(c("geojsonio"))

############################## ISSUES TO BE ADRESSED ##################################

#Different final area between copernicus and FABLE targets -> needs scaling

#######################################################################################

library(downscalr)

library(geojsonio)
library(here)
library(tidyr)
library(dplyr)
library(ggplot2)

#do I still use it???
# library(sp)
library(raster)


# Loading data ------------------------------------------------------------


population <- geojson_read(here("Data", "Population2020GRC.geojson"),  what = "sp")
population.df <- population@data

ProtectedAreas <- geojson_read(here("Data", "ProtectedAreasGRC.geojson"),  what = "sp")
ProtectedAreas.df <- ProtectedAreas@data

Travel <- geojson_read(here("Data", "TravelTimeGRC.geojson"),  what = "sp")
Travel.df <- Travel@data

LandCover <- geojson_read(here("Data", "LandCoverCopernicus2019GRC.geojson"),  what = "sp")
LandCover.df <- LandCover@data

LandCoverChange <- geojson_read(here("Data", "LandCoverChangeHILDA2018_2019GRC.geojson"),  what = "sp")
LandCoverChange.df <- LandCoverChange@data

ForestManagement <- geojson_read(here("Data", "ForestManagementGRC.geojson"),  what = "sp")
ForestManagement.df <- ForestManagement@data

#raster to use for plots
ext <- extent(Travel)
res <- 0.5
raster_layer <- raster(ext, res = res)
rasterized_layer <- rasterize(Travel, raster_layer, field = "fid")



##### WILL NEED TO REVIEW THE MAPPING!!
##### Was done quickly

map_HILDA <- readxl::read_excel("Data/mapping_code.xlsx", 
                        sheet = "HILDA")

map_Copernicus <- readxl::read_excel("Data/mapping_code.xlsx", 
                                sheet = "Copernicus")

# Changing structure ------------------------------------------------------

greece_luc <- LandCoverChange.df %>% 
  pivot_longer(-c(id, fid), values_to = "AreaPerCode") %>% 
  mutate(AreaPerCode = as.numeric(AreaPerCode)) %>% 
  #group by cell
  group_by(fid) %>% 
  #compute the total area in that cell
  mutate(TotalArea = sum(AreaPerCode, na.rm = T)) %>% 
  #compute the share of each land cover per cell
  mutate(value = AreaPerCode/TotalArea) %>% 
  mutate(name = stringr::str_remove(name, "X")) %>% 
  #associate land cover code with right name
  left_join(map_HILDA %>%  mutate(code = as.character(code)), by = c("name" = "code")) %>% 
  #set the right variable names to use in DownscalR package
  rename(ns = fid,
         lu.from = from,
         lu.to = to) %>% 
  mutate(Ts = 2018) %>% 
  dplyr::select(ns, lu.from, lu.to, Ts, value) %>% 
  mutate(value = ifelse(is.na(value), 0, value))

#Right now only has land cover data but should contain all the explanatory variables to compute the priors
Xmat <- LandCover.df %>% 
  pivot_longer(-c(id, fid), values_to = "AreaPerCover") %>%
  mutate(AreaPerCover = as.numeric(AreaPerCover))  %>% 
  mutate(name = stringr::str_remove(name, "X")) %>% 
  #associate the right land cover names
  left_join(map_Copernicus %>%  mutate(code = as.character(code)), by = c("name" = "code")) %>% 
  group_by(fid, LandCover) %>% 
  #compute the area per cover class in each cell
  summarise(value = sum(AreaPerCover, na.rm = T)) %>% 
  dplyr::select(fid, LandCover, value) %>% 
  #pivot wider to add the intercept value
  pivot_wider(names_from = "LandCover", values_from = "value") %>% 
  mutate(Intercept = 1) %>% 
  #pivot back longer to fit the DownscalR package structure
  pivot_longer(-fid, values_to = "value", names_to = "ks") %>% 
  rename(ns = fid)

lu_levels <- LandCover.df %>% 
  pivot_longer(-c(id, fid), values_to = "AreaPerCover") %>%
  mutate(AreaPerCover = as.numeric(AreaPerCover))  %>% 
  mutate(name = stringr::str_remove(name, "X")) %>% 
  #associate the right land cover names
  left_join(map_Copernicus %>%  mutate(code = as.character(code)), by = c("name" = "code"))%>% 
  filter(!(LandCover %in% c("ocean", "water"))) %>% 
  group_by(fid, LandCover) %>% 
  #Compute the area per cover class per cell
  summarise(value = sum(AreaPerCover, na.rm = T)) %>% 
  rename(ns = fid,
         lu.from = LandCover) %>% 
  dplyr::select(ns, lu.from, value)

greece_df <- list(xmat = tibble(Xmat), lu_levels = tibble(lu_levels))

# Following vignette code -------------------------------------------------

example_LU_from <- "cropland"

Yraw <- greece_luc 

Y <- dplyr::filter(Yraw,lu.from == example_LU_from & Ts == 2018) %>% #changes 2000 to 2010
  tidyr::pivot_wider(names_from = lu.to) %>%
  tibble::column_to_rownames(var = "ns") %>%
  dplyr::select(-c(lu.from, Ts))

X <- greece_df$xmat %>% tidyr::pivot_wider(names_from = "ks") %>%
  dplyr::arrange(match(ns,Y$ns)) %>%
  tibble::column_to_rownames(var = "ns") 

## Estimate MNL model

#Compute MNL model with standard settings
baseline <- which(colnames(Y) == example_LU_from)

# Decreased niter = 100 and nburn = 50 for faster computation 
results_MNL <- mnlogit(as.matrix(X), as.matrix(Y), baseline, niter = 100, nburn = 50, A0 = 10^4)
# Strongly recommended if you have more time: 
# results_MNL <- mnlogit(as.matrix(X), as.matrix(Y), baseline, niter = 2000, nburn = 1000, A0 = 10^4)

## Downscale module

### Preliminaries

# Set example year of downscaling target. e.g. 2010
example_time <- 2015

# Prepare inputs for downscale function
## get target data and filter for year and lu.from according to example, i.e. 2010 & cropland
GRC_FABLE <- readxl::read_excel("Data/GRC_FABLE.xlsx") %>% #In calculator, sheet 4_calc_land, 
  #table calc_landmatrix, columns LandCoverInit	YearStart	YearEnd	ToForest	ToOtherLand	ToCropland	ToPasture	ToUrban	ToNewForest
  dplyr::select(-YearStart) %>% 
  pivot_longer(-c(YearEnd, LandCoverInit), values_to = "value", names_to = "lu.to") %>% 
  mutate(lu.to = stringr::str_remove(lu.to, "To")) %>% 
  mutate(lu.from = tolower(LandCoverInit),
         lu.to = tolower(lu.to),
         times = YearEnd) %>% 
  dplyr::select(times, lu.from, lu.to, value) %>% 
  #do not keep land cover that remains the same
  filter(!(lu.from == lu.to)) %>% 
  #do not keep any land cover turning into forest because not possible
  filter(!(lu.from != "forest" & lu.to == "forest"))

#focus on the LU chosen before
# Targeted land use change from the chose Land cover to other land cover in the year selected
greece_targets_crop_2020 <- dplyr::filter(GRC_FABLE, lu.from == example_LU_from &
                                         times == example_time)

## pivot explanatory input data into long format
X_long <-  X %>% tibble::rownames_to_column(var="ns") %>%
  tidyr::pivot_longer(!ns, names_to="ks", values_to = "value")

## prior coefficients: compute mean of posterior draws & pivot to long format  
pred_coeff <- as.data.frame(apply(results_MNL$postb[,-baseline,],c(1,2),mean))

pred_coeff_long <- pred_coeff %>% tibble::rownames_to_column(var="ks") %>%
  tidyr::pivot_longer(!ks, names_to = "lu.to", values_to="value") %>% 
  dplyr::arrange(lu.to) %>% 
  tibble::add_column(lu.from=example_LU_from, .after="ks")

## starting area should provide information on area of each (higher resolution) grid cell, i.e ns
greece_start_areas_crop <- dplyr::filter(greece_df$lu_levels, lu.from == example_LU_from)

## Downscaling computation and results

results_DS <- downscale(targets = greece_targets_crop_2020, 
                        start.areas =  greece_start_areas_crop, 
                        xmat = X_long, 
                        betas = pred_coeff_long)

downscaled_LUC <- results_DS$out.res


# call plot function & plot downscaled LUC projections from Cropland for the year 2010
LUC_greece_2010_Crop_plot <- LUC_plot(results_DS, rasterized_layer)
LUC_greece_2010_Crop_plot$LUC.plot

for(cur_LC in colnames(Y)){
  for(cur_Year in seq(2020, 2050, 5)){
    
  }
}
