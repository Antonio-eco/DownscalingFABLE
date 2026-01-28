### Maps for restoration documentation FOLUR


library(geojsonio)
library(here)
library(raster)
library(dplyr)
library(tidyr)
library(downscalr)
library(stringr)
#library(kableExtra)
library(knitr)
library(readr)
library(forcats)
library(purrr)
library(writexl)
library(ggplot2)

if (requireNamespace("conflicted", quietly = TRUE)) {
  conflicted::conflict_prefer("filter",    "dplyr", quiet = TRUE)
  conflicted::conflict_prefer("select",    "dplyr", quiet = TRUE)
  conflicted::conflict_prefer("lag",       "dplyr", quiet = TRUE)
  conflicted::conflict_prefer("rename",    "dplyr", quiet = TRUE)
  conflicted::conflict_prefer("left_join", "dplyr", quiet = TRUE)
}

source(here("Program", "LUC_plot_cum_allLUC.R"))
source(here("Program", "LUC_plot_compare_pathways.R"))
source(here("Program", "LUC_plot_restriction.R"))
source(here("Program", "wrangle_plot.R"))
source(here("Program", "fct_IKI_India.R"))
source(here("Program", "LUC_plot_restoration_percent.R"))

date <- "251024"
country <- "IND"
stamp <-  date



## created on the R markdown
country_start_areas  <- read_rds(here("Output", country, paste0(stamp, "_country_start_areas.rds")))
X_long               <- read_rds(here("Output", country, paste0(stamp, "_X_long.rds")))
pred_coeff_long      <- read_rds(here("Output", country, paste0(stamp, "_pred_coeff_long.rds")))
rasterized_layer     <- terra::rast(here("Output", country, paste0("260128", "_IND_CT_BxG_HILDA", "_rasterized_layer.tif")))
start_map_reproj     <- read_rds(here("Output", country, paste0("251218", "_IND_CT_Bastin_HILDA", "_start_map_reproj.rds")))
start_map            <- read_rds(here("Output", country, paste0("251219", "_IND_CT_BxG_HILDA", "_start_map.rds")))
grid_sf              <- read_rds(here("Output", country, paste0("260128", "_IND_CT_BxG_HILDA", "_grid_sf.rds")))

ns_map               <- read_rds(here("Output", country, paste0("251218", "_IND_CT_Bastin_HILDA", "_ns_map.rds")))
grid <- read_csv(here("Data", "global", "grid50_equal_area.csv")) %>% 
  filter(iso3 == country)  

FABLE_CT            <- read_rds(here("Output", country, paste0("251020", "_IND_CurrentTrends_HILDA_FABLE.rds")))

l_pathways <- c("FABLE_CT")

for(cur.path in l_pathways){
  
  temp <- get(cur.path) %>% filter(times >= 2020)
  temp_results <- downscale(
    targets      = temp,
    start.areas  = country_start_areas,
    xmat         = X_long,
    betas        = pred_coeff_long
  )
  
  temp_results$out.res <- format_to_id_map(temp_results, ns_map)
  
  assign(paste0("results_", cur.path), temp_results)
  print(paste0("Generated:  ", "results_", cur.path))
}


## Restoration maps ------

#Global tree cover (existing and potential) (%) (remove existing? which map is used for existing?)
rest_bastin.df <- geojson_read(here("Data", country, "GlobalTreeRestorationPotential_Bastin.geojson"),  what = "sp")@data  %>% 
  filter(!is.na(mean)) %>% 
  left_join(grid %>% select(id_c, area)) %>% 
  mutate(p = mean) %>% 
  select(-id) 

# #kilograms of CO2 equivalent per year (kgCO2e yr-1)
# #transformed into t CO2/ 1000 ha
# rest_griscom <- geojson_read(here("Data", country, "ReforestationPotential_Griscom.geojson"),  what = "sp")
# rest_griscom.df <- rest_griscom@data  %>% 
#   select(-id) %>% 
#   left_join(grid %>% select(id_c, area)) %>% 
#   mutate(p = (sum / 1000)/(area/1000)) #p in t CO2/ 1000 ha

#mapping
map_chaturvedi <- readxl::read_excel(here("Data/mapping_code.xlsx"), 
                                     sheet = "Chaturvedi")

#Landscape restoration opportunities in India (absolute hectares)
rest_chaturvedi.df <- geojson_read(here("Data", 
                                        country, 
                                        "LandscapeRestorationOpportunities.geojson"),  
                                   what = "sp")@data %>% 
  select(-id) %>% 
  pivot_longer(-c(id_c)) %>% 
  mutate(code = stringr::str_remove(name, "X")) %>%
  left_join(map_chaturvedi %>% mutate(code = as.character(code))) %>% 
  left_join(grid %>% select(id_c, area)) %>% 
  #use the area in the cell to compute % per cell
  mutate(p = pmin(100, 100 * value/area))


# Scenario based total  restoration potential area (hectares)
rest_fesenmyer.df <- geojson_read(here("Data", 
                                       country, 
                                       "ScenarioBasedTotalRestorationPotentialArea_Fesenmyer.geojson"),  
                                  what = "sp")@data %>% 
  select(id_c, d0) %>% 
  left_join(grid %>% select(id_c, area)) %>% 
  #use the area in the cell to compute % per cell
  mutate(p = pmin(100, 100 * d0/area))

#Area with potential for forest regeneration per pixel
#Values are given in hectares
rest_williams.df <- geojson_read(here("Data", 
                                      country, 
                                      "PotentialForestRegenerationTropicalForest_Williams.geojson"),  
                                 what = "sp")@data  %>% 
  select(-id) %>% 
  left_join(grid %>% select(id_c, area)) %>% 
  #use the area in the cell to compute % per cell
  mutate(p = pmin(100, 100 * sum/area))


########## plot maps

## Bastin

results_DS_bastin <- results_FABLE_CT
results_DS_bastin$out.res <- rest_bastin.df %>% 
  left_join(ns_map %>% 
              select(-ns)) %>% 
  rename(ns = ns_int,
         value = p) %>% 
  mutate(lu.from = "any",
         lu.to = "newforest",
         times = 2050) %>% 
  mutate(value = ifelse(value < 5, 0, value)) %>% #greater than 5% to appear
select(ns, times, lu.from, lu.to, value)



## Williams

results_DS_williams <- results_FABLE_CT
results_DS_williams$out.res <- rest_williams.df %>% 
  left_join(ns_map %>% 
              select(-ns)) %>% 
  rename(ns = ns_int,
         value = p) %>% 
  mutate(lu.from = "any",
         lu.to = "newforest",
         times = 2050) %>% 
  mutate(value = ifelse(value < 5, 0, value)) %>% #greater than 5% to appear
  select(ns, times, lu.from, lu.to, value)




## Fesenmyer

results_DS_fesenmyer <- results_FABLE_CT
results_DS_fesenmyer$out.res <- rest_fesenmyer.df %>% 
  left_join(ns_map %>% 
              select(-ns)) %>% 
  rename(ns = ns_int,
         value = p) %>% 
  mutate(lu.from = "any",
         lu.to = "newforest",
         times = 2050) %>% 
  mutate(value = ifelse(value < 5, 0, value)) %>% #greater than 5% to appear
  select(ns, times, lu.from, lu.to, value)



## Chaturvedi - wide scale restoration


results_DS_chaturvedi_wide <- results_FABLE_CT
results_DS_chaturvedi_wide$out.res <- rest_chaturvedi.df %>% 
  select(-value) %>% 
  filter(restoration_type == "wide_scale_restoration") %>% 
  left_join(ns_map %>% 
              select(-ns)) %>% 
  rename(ns = ns_int,
         value = p) %>% 
  mutate(lu.from = "any",
         lu.to = "newforest",
         times = 2050) %>% 
  mutate(value = ifelse(value < 5, 0, value)) %>% #greater than 5% to appear
  select(ns, times, lu.from, lu.to, value)



## Chaturvedi - mosaic scale restoration


results_DS_chaturvedi_mosaic <- results_FABLE_CT
results_DS_chaturvedi_mosaic$out.res <- rest_chaturvedi.df %>% 
  select(-value) %>% 
  filter(restoration_type == "mosaic_scale_restoration") %>% 
  left_join(ns_map %>% 
              select(-ns)) %>% 
  rename(ns = ns_int,
         value = p) %>% 
  mutate(lu.from = "any",
         lu.to = "newforest",
         times = 2050) %>% 
  mutate(value = ifelse(value < 5, 0, value)) %>% #greater than 5% to appear
  select(ns, times, lu.from, lu.to, value)

##limits

limits <- range(
  results_DS_bastin$out.res$value,
  results_DS_williams$out.res$value,
  results_DS_fesenmyer$out.res$value,
  results_DS_chaturvedi_mosaic$out.res$value,
  results_DS_chaturvedi_wide$out.res$value,
  na.rm = TRUE
)


##### plot

plot_bastin_restricted <- LUC_plot_restoration_percent(results_DS_bastin, 
                                                       rasterized_layer, 
                                                       grid_sf,
                                                       limits,
                                                       label = "Bastin et al. (2019)\nPercentage of pixel prioritised\nfor forest restoration:")

ggplot2::ggsave(
  filename = here("Output", country, paste0(format(Sys.Date(), "%y%m%d"), "_Bastin_map.jpeg")),
  plot     = plot_bastin_restricted$LUC.plot,
  device   = "jpeg",
  width    = 5,
  height   = 6,
  units    = "in",
  dpi      = 300
)

plot_williams_restricted <- LUC_plot_restoration_percent(results_DS_williams, 
                                                         rasterized_layer, 
                                                         grid_sf,
                                                         limits,
                                                         label = "Williams et al. (2024)\nPercentage of pixel prioritised for\nforest regen. in tropical regions:")

ggplot2::ggsave(
  filename = here("Output", country, paste0(format(Sys.Date(), "%y%m%d"), "_Williams_map.jpeg")),
  plot     = plot_williams_restricted$LUC.plot  ,
  device   = "jpeg",
  width    = 5,
  height   = 6,
  units    = "in",
  dpi      = 300
)

plot_fesenmyer_restricted <- LUC_plot_restoration_percent(results_DS_fesenmyer, 
                                                          rasterized_layer, 
                                                          grid_sf,
                                                          limits,
                                                          label = "Fesenmyer et al (2025)\nPercentage of pixel prioritised\nfor forest restoration:")



ggplot2::ggsave(
  filename = here("Output", country, paste0(format(Sys.Date(), "%y%m%d"), "_Fesenmyer_map.jpeg")),
  plot     = plot_fesenmyer_restricted$LUC.plot ,
  device   = "jpeg",
  width    = 5,
  height   = 6,
  units    = "in",
  dpi      = 300
)

plot_chaturvedi_wide_restricted <- LUC_plot_restoration_percent(results_DS_chaturvedi_wide, 
                                                                rasterized_layer, 
                                                                grid_sf,
                                                                limits,
                                                                label = "Chaturvedi et al. (2018)\nPercentage of pixel prioritised\nfor wide scale restoration:")


ggplot2::ggsave(
  filename = here("Output", country, paste0(format(Sys.Date(), "%y%m%d"), "_Chaturvedi_wide_map.jpeg")),
  plot     = plot_chaturvedi_wide_restricted$LUC.plot ,
  device   = "jpeg",
  width    = 5,
  height   = 6,
  units    = "in",
  dpi      = 300
)

plot_chaturvedi_mosaic_restricted <- LUC_plot_restoration_percent(results_DS_chaturvedi_mosaic, 
                                                                  rasterized_layer, 
                                                                  grid_sf,
                                                                  limits,
                                                                  label = "Chaturvedi et al. (2018)\nPercentage of pixel prioritised\nfor mosaic scale restoration:")


ggplot2::ggsave(
  filename = here("Output", country, paste0(format(Sys.Date(), "%y%m%d"), "_Chaturvedi_mosaic_map.jpeg")),
  plot     = plot_chaturvedi_mosaic_restricted$LUC.plot ,
  device   = "jpeg",
  width    = 5,
  height   = 6,
  units    = "in",
  dpi      = 300
)