## IKI-LTS results on India
## Clara Douzal
## First created 

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

date <- "251024"
country <- "IND"
stamp <-  date



## created on the R markdown
country_start_areas  <- read_rds(here("Output", country, paste0(stamp, "_country_start_areas.rds")))
X_long               <- read_rds(here("Output", country, paste0(stamp, "_X_long.rds")))
pred_coeff_long      <- read_rds(here("Output", country, paste0(stamp, "_pred_coeff_long.rds")))
rasterized_layer     <- read_rds(here("Output", country, paste0("251218", "_IND_CT_Bastin_HILDA", "_rasterized_layer.rds")))
start_map_reproj     <- read_rds(here("Output", country, paste0("251218", "_IND_CT_Bastin_HILDA", "_start_map_reproj.rds")))
start_map            <- read_rds(here("Output", country, paste0("251219", "_IND_CT_BxG_HILDA", "_start_map.rds")))

ns_map               <- read_rds(here("Output", country, paste0("251218", "_IND_CT_Bastin_HILDA", "_ns_map.rds")))
grid <- read_csv(here("Data", "global", "grid50_equal_area.csv")) %>% 
  filter(iso3 == country)  

FABLE_CT            <- read_rds(here("Output", country, paste0("251020", "_IND_CurrentTrends_HILDA_FABLE.rds")))
FABLE_CT_Bastin     <- read_rds(here("Output", country, paste0("251217", "_IND_CT_Bastin_HILDA_FABLE.rds")))
FABLE_CT_BxG        <- read_rds(here("Output", country, paste0("251219", "_IND_CT_BxG_HILDA_FABLE.rds")))
#FABLE_CT_Fesenmyer  <- read_rds(here("Output", country, paste0("251217", "_IND_CT_Fesenmyer_HILDA_FABLE.rds")))

l_pathways <- c("FABLE_CT", "FABLE_CT_BxG")

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


## Restrictions------

#Global tree cover (existing and potential) (%) (remove existing? which map is used for existing?)
rest_bastin <- geojson_read(here("Data", country, "GlobalTreeRestorationPotential_Bastin.geojson"),  what = "sp")
rest_bastin.df <- rest_bastin@data  %>% 
  filter(!is.na(mean)) %>% 
  left_join(grid %>% select(id_c, area)) %>% 
  mutate(rest_area = pmin(area, (mean/100)*area)) %>% 
  select(-id) 

#kilograms of CO2 equivalent per year (kgCO2e yr-1)
#transformed into t CO2/ 1000 ha
rest_griscom <- geojson_read(here("Data", country, "ReforestationPotential_Griscom.geojson"),  what = "sp")
rest_griscom.df <- rest_griscom@data  %>% 
  select(-id) %>% 
  left_join(grid %>% select(id_c, area)) %>% 
  mutate(p = (sum / 1000)/(area/1000)) #p in t CO2/ 1000 ha
 

## Bastin restriction -----

bastin_treecover <- remove_existingForest(rest_bastin.df, ns_map, start_map_reproj)

q75 <- quantile(bastin_treecover[["potential_area"]], prob=c(.75), type=1)

colSums(bastin_treecover %>% select(potential_area) %>% filter(potential_area >30.77))

bastin_restriction <- bastin_treecover %>% 
  filter(potential_area >q75) %>% 
  select(ns)

bastin_rest_worthy     <- bastin_restriction
bastin_not_rest_worthy <- setdiff(ns_map["ns"], bastin_rest_worthy)

## plot to see where the cells are
bastin_notworthy <- ns_map %>% 
  filter(ns %in% bastin_rest_worthy[[1]]) %>% 
  mutate(value = 1) %>% 
  rbind(ns_map %>% 
          filter(ns %in% bastin_not_rest_worthy[[1]])
        %>% mutate(value = -1))%>%
  mutate(lu.from = "any",
         lu.to = "newforest",
         times = 2025) %>% 
  select(-ns) %>% 
  rename(ns = ns_int) %>% 
  select(ns, times, lu.from, lu.to, value)

results_DS_bastin <- results_FABLE_CT
results_DS_bastin$out.res <- bastin_notworthy

plot_bastin_restricted <- LUC_plot_restriction(results_DS_bastin, rasterized_layer, 
                                               color = "PRGn", label = "Priority areas for restoration:")

plot_bastin_restricted$LUC.plot  


#### close ------
# tiff(
#   here("Output", country, paste0(tag, "_Bastin_restriction_map.tiff")),
#   units = "in", height = 6, width = 10, res = 300
# )
# plot_bastin_restricted$LUC.plot
# dev.off()


## Griscom -----

#Use Griscom's global restoration map in t CO2/ 1000 ha
q75 <- quantile(rest_griscom.df[["p"]], prob=c(.75), type=1)

griscom_not_restoration_worthy <- rest_griscom.df %>% 
  filter(p < q75) %>% 
  merge(ns_map) %>% 
  select(ns)

griscom_rest_worthy <- setdiff(ns_map["ns"], griscom_not_restoration_worthy)

griscom_notworthy <- ns_map %>% 
  filter(ns %in% griscom_rest_worthy[[1]]) %>% 
  mutate(value = 1) %>% 
  rbind(ns_map %>% 
          filter(ns %in% griscom_not_restoration_worthy[[1]])
        %>% mutate(value = -1))%>%
  mutate(lu.from = "any",
         lu.to = "newforest",
         times = 2025) %>% 
  select(-ns) %>% 
  rename(ns = ns_int) %>% 
  select(ns, times, lu.from, lu.to, value)

results_DS_griscom <- results_FABLE_CT
results_DS_griscom$out.res <- griscom_notworthy

plot_griscom_restricted <- LUC_plot_restriction(results_DS_griscom, rasterized_layer, 
                                        color = "PRGn", label = "Priority areas for restoration:")

plot_griscom_restricted$LUC.plot

## overlap bastin and griscom

overlap <- griscom_notworthy %>% 
  rename(g_value = value) %>% 
  left_join(
    bastin_notworthy %>% 
      rename(b_value = value)
  ) %>% 
  mutate(value = g_value + b_value) %>% 
  mutate(value = ifelse(value == 2, 1, -1))

results_DS_overlap <- results_FABLE_CT
results_DS_overlap$out.res <- overlap

plot_overlap_restricted <- LUC_plot_restriction(results_DS_overlap, rasterized_layer, 
                                                color = "PRGn", label = "Priority areas for restoration:")

plot_overlap_restricted$LUC.plot


overlap_area <- overlap  %>% 
  left_join(bastin_treecover%>%
              select(ns_int, potential_area) %>% 
              rename(ns = ns_int)) %>% 
  filter(value == 1)

colSums(overlap_area %>% select(potential_area))
# potential_area 
# 4263.46 


# downscale with BxG target and constraint -----------------------------------------------

#Use less restrictive Griscom's global restoration map in t CO2/ 1000 ha
# If I keep the same restriction on Griscom I cannot downscale my pathway because it is too restrictive.
q75 <- quantile(rest_griscom.df[["p"]], prob=c(.4), type=1)

griscom_not_restoration_worthy <- rest_griscom.df %>% 
  filter(p < q75) %>% 
  merge(ns_map) %>% 
  select(ns)

griscom_rest_worthy <- setdiff(ns_map["ns"], griscom_not_restoration_worthy)

griscom_notworthy <- ns_map %>% 
  filter(ns %in% griscom_rest_worthy[[1]]) %>% 
  mutate(value = 1) %>% 
  rbind(ns_map %>% 
          filter(ns %in% griscom_not_restoration_worthy[[1]])
        %>% mutate(value = -1))%>%
  mutate(lu.from = "any",
         lu.to = "newforest",
         times = 2025) %>% 
  select(-ns) %>% 
  rename(ns = ns_int) %>% 
  select(ns, times, lu.from, lu.to, value)

overlap <- griscom_notworthy %>% 
  rename(g_value = value) %>% 
  left_join(
    bastin_notworthy %>% 
      rename(b_value = value)
  ) %>% 
  mutate(value = g_value + b_value) %>% 
  mutate(value = ifelse(value == 2, 1, -1))

BxG_restworthy <- setdiff(ns_map["ns_int"], overlap %>% filter(value == -1) %>% select(ns) %>% rename(ns_int = ns)) %>% 
  left_join(ns_map) %>% 
  select(ns)


restriction_BxG <- data.frame(
  expand.grid(ns               = unique(start_map$ns),
              lu.from          = c("cropland", "forest", "otherland", "pasture", "urban"),
              lu.to            = c("cropland", "newforest", "otherland", "pasture", "urban", "forest"),
              stringsAsFactors = FALSE) %>%
    filter(!(lu.from == "forest" & lu.to == "newforest")) %>%
    filter(!(lu.from != "forest" & lu.to == "forest")) %>%
    filter(!(lu.from == lu.to)) %>%
    #allowed
    mutate(value = 0L) %>%
    #not allowed
    mutate(value = ifelse((ns %in% (setdiff(ns_map["ns"], BxG_restworthy))$ns & lu.to == "newforest"), 1L, value))
)

## Bastin potential using FABLE output and restriction

results_DS_rest <- downscale(
  targets      = FABLE_CT_BxG %>% filter(times >= 2020),
  start.areas  = country_start_areas,
  xmat         = X_long,
  betas        = pred_coeff_long,
  restrictions = restriction_BxG
)

results_DS_rest$out.res <- format_to_id_map(results_DS_rest, ns_map)

plot_BxG_restoration <- LUC_plot(results_DS_rest, rasterized_layer,
                                    year = 2025, LU = "newforest",
                                    color = "Blues", label = "Restoration tree cover potential (kha):")
# tiff(
#   here("Output", country, paste0(tag, "_BxG_potential_map.tiff")),
#   units = "in", height = 6, width = 10, res = 300
# )

plot_BxG_restoration$LUC.plot +
  ggplot2::theme(legend.text = element_text(size = 14))+
  ggplot2::theme(legend.title = element_text(size = 16))

dev.off()

reforestation_area_BxG <- data.frame(results_DS_rest$out.res %>% group_by(times, lu.to, lu.from) %>% summarise(value = sum(value)) %>% filter(lu.to == "newforest" & lu.from != "newforest"))

colSums(reforestation_area_BxG %>% select(value))
# There is not enough otherland left in 2050 to meet the afforestation target
# this is not a probleme in the FABLE -C bu one here because we use a different we restrict 
# value 
# 3687.091 

##refined indicator -------------------------------------------------------

## Refine indicator using Griscom

Seq_per_area <- rest_griscom.df %>%  #tCO2e/1000ha 
  left_join(ns_map) %>% 
  select(ns_int, p) %>% #p in t CO2/ 1000 ha
  rename(ns = ns_int)

vect_years <-  data.frame(
  times = as.character(c(2025, 2030, 2035, 2040, 2045, 2050)),
  factor = c(25, 20, 15, 10, 5, 1)
)

BxG_seq_cell_per_5year <- results_DS_rest$out.res %>% #in hectares
  filter(lu.to == "newforest") %>%
  left_join(Seq_per_area) %>% #in t CO2e / 1000 ha
  mutate(seq_per_cell = p/50*5 * value)  #in t CO2e

## Would need to restructure the format to not have directly full potential but something similar to the FABLE-C

# BxG_seq_year <- BxG_seq_cell %>% 
#   filter(lu.from != "newforest") %>% 
#   left_join(vect_years) %>% 
#   mutate(seq_per_cell = seq_per_cell/50*5) %>% 
#   group_by(times) %>% 
#   summarise(seq_per_year = sum(seq_per_cell))

BxG_seq_potential <- results_DS_rest$out.res %>% #in hectares
  filter(lu.to == "newforest") %>%
  group_by( ns) %>% 
  summarise(value = sum(value, na.rm = T)) %>% 
  left_join(Seq_per_area) %>%
  mutate(seq_pot = value * p /value /50 /1000)

seq_rates_stats <- as.data.frame(t(summary(BxG_seq_potential)$p))

summary(BxG_seq_potential$p)
summary(BxG_seq_potential$seq_pot)


