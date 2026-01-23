
library(here)
library(scales)
library(ggnewscale)

date <- "251020"
country <- "IND"
stamp <-  date


# results_FABLE_CT <- read_rds(here("Output", country, paste0(date, "_", country, "_CurrentTrends_HILDA_downscaled_LUC.rds")))
# results_FABLE_CT_Bastin <- read_rds(here("Output", country, paste0(date, "_", country, "_CT_Bastin_HILDA_downscaled_LUC.rds")))
# results_FABLE_CT_Potapov <- read_rds(here("Output", country, paste0(date, "_", country, "_CT_Potapov_HILDA_downscaled_LUC.rds")))
# results_FABLE_NC_Bastin <- read_rds(here("Output", country, paste0(date, "_", country, "_NC_Bastin_HILDA_downscaled_LUC.rds")))
# results_FABLE_NC_Potapov <- read_rds(here("Output", country, paste0(date, "_", country, "_NC_Potapov_HILDA_downscaled_LUC.rds")))

country_start_areas <- read_rds(here("Output", country, paste0(stamp, "_country_start_areas.rds")))
X_long <- read_rds(here("Output", country, paste0(stamp, "_X_long.rds")))
pred_coeff_long <- read_rds(here("Output", country, paste0(stamp, "_pred_coeff_long.rds")))

FABLE_CT <- read_rds(here("Output", country, paste0(stamp, "_IND_CurrentTrends_HILDA_FABLE.rds")))
FABLE_CT_Bastin <- read_rds(here("Output", country, paste0(stamp, "_IND_CT_Bastin_HILDA_FABLE.rds")))
FABLE_CT_Potapov <- read_rds(here("Output", country, paste0(stamp, "_IND_CT_Potapov_HILDA_FABLE.rds")))
FABLE_NC_Bastin <- read_rds(here("Output", country, paste0(stamp, "_IND_NC_Bastin_HILDA_FABLE.rds")))
FABLE_NC_Potapov <- read_rds(here("Output", country, paste0(stamp, "_IND_NC_Potapov_HILDA_FABLE.rds")))

l_pathways <- c("FABLE_CT", "FABLE_CT_Bastin", "FABLE_CT_Potapov", "FABLE_NC_Bastin", "FABLE_NC_Potapov")


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

### Restrictions

#Use Griscom's global restoration map in t CO2/ 1000 ha
q75 <- quantile(rest_griscom.df[["p"]], prob=c(.60), type=1)

#cells that are not restoration worthy
not_restoration_worthy <- rest_griscom.df %>% 
  filter(p < q75) %>% 
  merge(ns_map) %>% 
  select(ns)
#cells that are restoration worthy
rest_worthy <- setdiff(ns_map["ns"], not_restoration_worthy)

## plot to see where the cells are
test_notworthy <- ns_map %>% 
  filter(ns %in% rest_worthy[[1]]) %>% 
  mutate(value = 1) %>% 
  rbind(ns_map %>% 
          filter(ns %in% not_restoration_worthy[[1]])
        %>% mutate(value = -1))%>%
  mutate(lu.from = "any",
         lu.to = "newforest",
         times = 2025) %>% 
  select(-ns) %>% 
  rename(ns = ns_int) %>% 
  select(ns, times, lu.from, lu.to, value)

results_DS_test <- results_FABLE_CT
results_DS_test$out.res <- test_notworthy

plot_restricted <- LUC_plot_restriction(results_DS_test, rasterized_layer, 
                                        color = "PRGn", label = "Priority areas for restoration:")

plot_restricted$LUC.plot

ggplot2::ggsave(
  filename = here("Output", country, paste0(tag, "restriction_map_griscom.tiff")), 
  plot = plot_restricted$LUC.plot, 
  units = "in", 
  height = 6, width = 10, dpi = 300)

## /!\ /!\ /!\ /!\ /!\  
# restriction table needs to include exactly the same (lu.from, lu.to) transition pairs as in targets 
# (excluding diagonals where lu.from == lu.to). Using any subset/superset leads to errors or unexpected behavior.
## /!\ /!\ /!\ /!\ /!\ 

restriction <- data.frame(
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
    mutate(value = ifelse((ns %in% not_restoration_worthy[[1]] & lu.to == "newforest"), 1L, value)) 
)


area_target_restoration <- data.frame()

for(cur.path in l_pathways){
  
  temp <- get(cur.path) %>% 
    filter(lu.to   == "newforest") %>% 
    filter(lu.from != "urban") %>% 
    pivot_wider(names_from = lu.from, values_from = value) 
  
  area_target_restoration <- rbind.data.frame(
    area_target_restoration,
    cbind(pathway = cur.path, t(colSums(temp %>% select(-c(times, lu.to)))))
  )
}

cells_notworthy <- ns_map %>% 
  filter(ns %in% unique(c(not_restoration_worthy[[1]]))) %>% 
  select(ns_int) %>% 
  rename(ns = ns_int) %>% 
  mutate(worthy = 0)

cells_worthy <- ns_map %>% 
  filter(!(ns %in% unique(c(not_restoration_worthy[[1]])))) %>% 
  select(ns_int) %>% 
  rename(ns = ns_int) %>% 
  mutate(worthy = 1)

area_restriction <- start_map_reproj %>% 
  left_join(ns_map) %>% 
  select(ns_int, lu.from, value) %>% 
  rename(ns = ns_int) %>% 
  left_join(rbind(cells_notworthy, cells_worthy)) %>% 
  filter(worthy == 1) %>% 
  filter(lu.from != "forest") %>% 
  filter(lu.from != "urban") %>% 
  pivot_wider(names_from = lu.from, values_from = value)

enough_land <- rbind.data.frame(
  area_target_restoration, 
  cbind(pathway = "restriction", t(colSums(area_restriction %>% select(cropland, otherland, pasture)))))


for(cur.path in l_pathways){
  
  temp <- get(cur.path) %>% filter(times >= 2020)
  temp_results <- downscale(
    targets      = temp,
    start.areas  = country_start_areas,
    xmat         = X_long,
    betas        = pred_coeff_long,
    restrictions = restriction   
  )
  
  temp_results$out.res <- format_to_id_map(temp_results, ns_map)
  
  assign(paste0("results_", cur.path, "_restriction1"), temp_results)
  print(paste0("Generated:  ", "results_", cur.path, "_restriction1"))
}


# Aggregate allocation vs targets

# helper to normalize times consistently
norm_year <- function(x) as.integer(as.character(x))

alloc_agg_fct <- function(res) {
  
  temp <- res$out.res
  
  alloc <- temp %>% 
    mutate(times = norm_year(times)) %>%
    group_by(times, lu.from, lu.to) %>%
    summarise(alloc = sum(value, na.rm = TRUE), .groups = "drop")
  
  return(alloc)
}

consistency_fct <- function(alloc, targ){
  
  consistency <- full_join(alloc, targ, by = c("times","lu.from","lu.to")) %>%
    mutate(across(c(alloc, target), ~ replace_na(.x, 0)),
           pct_err = 100 * (alloc - target) / pmax(abs(target), 1e-9)) %>%
    filter(lu.from != lu.to)
  
  # if you filter by min/max year anywhere, make sure you're numeric:
  min_year <- min(consistency$times, na.rm = TRUE)
  print(
    knitr::kable(consistency %>% filter(times >= min_year)) 
  )
  
  
}


for(cur.path in l_pathways){
  
  
  temp_rest <- get(paste0("results_", cur.path, "_restriction1"))
  
  alloc_agg_rest <- alloc_agg_fct(temp_rest)
  
  targ_agg <- get(cur.path) %>%
    mutate(times = norm_year(times)) %>%           # <- force same type as alloc_agg
    filter(times >= 2020) %>%
    group_by(times, lu.from, lu.to) %>%
    summarise(target = sum(value, na.rm = TRUE), .groups = "drop")
  
  # now the join will match on identical types
  print(cur.path)
  summary(consistency_fct(alloc_agg_rest, targ_agg))
}

### plot of all pathways x restriction maps

results_several_pathways <- list(results_FABLE_CT, 
                                 results_FABLE_CT_Bastin, 
                                 results_FABLE_CT_Potapov,
                                 results_FABLE_NC_Bastin,
                                 results_FABLE_NC_Potapov)
  

LUC_plot_pathways <- LUC_plot_compare_pathways(results_several_pathways, rasterized_layer, grid, l_pathways)
LUC_plot_pathways$LUC.plot

results_several_pathways_restriction <- list(results_FABLE_CT_restriction1, 
                                 results_FABLE_CT_Bastin_restriction1, 
                                 results_FABLE_CT_Potapov_restriction1,
                                 results_FABLE_NC_Bastin_restriction1,
                                 results_FABLE_NC_Potapov_restriction1)

LUC_plot_pathways <- LUC_plot_compare_pathways_percent(results_several_pathways_restriction, rasterized_layer, grid, ns_map, l_pathways)
LUC_plot_pathways$LUC.plot

ggplot2::ggsave(
  filename = here("Output", country, paste0(tag, "2020_2050_per_pathway_cumulative_final_Land-use_change-map.tiff")), 
  plot = LUC_plot_pathways$LUC.plot, 
  units = "in", 
  height = 6, width = 10, dpi = 300)


################################################################################
### Stricter restriction #######################################################
################################################################################

#Use Griscom's global restoration map in t CO2/ 1000 ha
q60_griscom <- quantile(rest_griscom.df[["p"]], prob=c(.45), type=1)

not_restoration_worthy_griscom <- rest_griscom.df %>% 
  filter(p < q60_griscom) %>% 
  merge(ns_map) %>% 
  select(ns)

#Use walker's unrealized potenial for increased storage of carbon on land
q60_walker <- quantile(rest_walker.df[["sum"]], prob=c(.45), type=1)

not_restoration_worthy_walker <- rest_walker.df %>% 
  filter(sum < q60_walker) %>% 
  merge(ns_map) %>% 
  select(ns)

#cells that are restoration worthy
rest_worthy <- setdiff(setdiff(ns_map["ns"], not_restoration_worthy_griscom),
                       not_restoration_worthy_walker)


## plot to see where the cells are
test_notworthy <- ns_map %>% 
  filter(ns %in% rest_worthy[[1]]) %>% 
  mutate(value = 1) %>% 
  rbind(ns_map %>% 
          filter(!(ns %in% rest_worthy[[1]]))
        %>% mutate(value = -1))%>%
  mutate(lu.from = "any",
         lu.to = "newforest",
         times = 2025) %>% 
  select(-ns) %>% 
  rename(ns = ns_int) %>% 
  select(ns, times, lu.from, lu.to, value)

results_DS_test <- results_FABLE_CT
results_DS_test$out.res <- test_notworthy

plot_restricted <- LUC_plot(results_DS_test, rasterized_layer)

plot_restricted$LUC.plot

## /!\ /!\ /!\ /!\ /!\  
# restriction table needs to include exactly the same (lu.from, lu.to) transition pairs as in targets 
# (excluding diagonals where lu.from == lu.to). Using any subset/superset leads to errors or unexpected behavior.
## /!\ /!\ /!\ /!\ /!\ 

restriction_full <- data.frame(
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
    mutate(value = ifelse((ns %in% (setdiff(ns_map["ns"], rest_worthy))$ns & lu.to == "newforest"), 1L, value)) 
)


area_target_restoration <- data.frame()

for(cur.path in l_pathways){
  
  temp <- get(cur.path) %>% 
    filter(lu.to   == "newforest") %>% 
    filter(lu.from != "urban") %>% 
    pivot_wider(names_from = lu.from, values_from = value) 
  
  area_target_restoration <- rbind.data.frame(
    area_target_restoration,
    cbind(pathway = cur.path, t(colSums(temp %>% select(-c(times, lu.to)))))
  )
}

cells_notworthy <- ns_map %>% 
  filter(!(ns %in% rest_worthy[[1]])) %>% 
  select(ns_int) %>% 
  rename(ns = ns_int) %>% 
  mutate(worthy = 0)

cells_worthy <- ns_map %>% 
  filter(ns %in% rest_worthy[[1]]) %>% 
  select(ns_int) %>% 
  rename(ns = ns_int) %>% 
  mutate(worthy = 1)

area_restriction <- start_map_reproj %>% 
  left_join(ns_map) %>% 
  select(ns_int, lu.from, value) %>% 
  rename(ns = ns_int) %>% 
  left_join(rbind(cells_notworthy, cells_worthy)) %>% 
  filter(worthy == 1) %>% 
  filter(lu.from != "forest") %>% 
  filter(lu.from != "urban") %>% 
  pivot_wider(names_from = lu.from, values_from = value)

enough_land_2 <- rbind.data.frame(
  area_target_restoration, 
  cbind(pathway = "restriction", t(colSums(area_restriction %>% select(cropland, otherland, pasture)))))


for(cur.path in l_pathways){
  
  temp <- get(cur.path) %>% filter(times >= 2020)
  temp_results <- downscale(
    targets      = temp,
    start.areas  = country_start_areas,
    xmat         = X_long,
    betas        = pred_coeff_long,
    restrictions = restriction_full   
  )
  
  temp_results$out.res <- format_to_id_map(temp_results, ns_map)
  
  assign(paste0("results_", cur.path, "_restriction_full"), temp_results)
  print(paste0("Generated:  ", "results_", cur.path, "_restriction_full"))
}

for(cur.path in l_pathways){
  
  
  temp_rest <- get(paste0("results_", cur.path, "_restriction_full"))
  
  alloc_agg_rest <- alloc_agg_fct(temp_rest)
  
  targ_agg <- get(cur.path) %>%
    mutate(times = norm_year(times)) %>%           # <- force same type as alloc_agg
    filter(times >= 2020) %>%
    group_by(times, lu.from, lu.to) %>%
    summarise(target = sum(value, na.rm = TRUE), .groups = "drop")
  
  # now the join will match on identical types
  print(cur.path)
  summary(consistency_fct(alloc_agg_rest, targ_agg))
}

### plot of all pathways x restriction maps


results_several_pathways_restriction_full <- list(results_FABLE_CT_restriction_full, 
                                             results_FABLE_CT_Bastin_restriction_full, 
                                             results_FABLE_CT_Potapov_restriction_full,
                                             results_FABLE_NC_Bastin_restriction_full,
                                             results_FABLE_NC_Potapov_restriction_full)

LUC_plot_pathways_full <- LUC_plot_compare_pathways_percent(results_several_pathways_restriction_full, 
                                                       rasterized_layer, grid, ns_map, l_pathways)
LUC_plot_pathways_full$LUC.plot

# ggplot2::ggsave(
#   filename = here("Output", country, paste0(tag, "2020_2050_per_pathway_cumulative_final_Land-use_change-map.tiff")), 
#   plot = LUC_plot_pathways$LUC.plot, 
#   units = "in", 
#   height = 6, width = 10, dpi = 300)

################################################################################
#################### Chaturvedi ##################################
################################################################################

target_chaturvedi <- rest_chaturvedi.df %>% 
  filter(restoration_type == "wide_scale_restoration") %>% 
  mutate(restoration_area = p*area/100) %>% 
  left_join(ns_map) %>% 
  left_join(start_map_reproj %>% filter(lu.from == "cropland") %>% mutate(cropland = value) %>% select(-c(lu.from, value))) %>%
  left_join(start_map_reproj %>% filter(lu.from == "urban") %>% mutate(urban = value) %>% select(-c(lu.from, value))) %>%
  left_join(start_map_reproj %>% filter(lu.from == "newforest") %>% mutate(newforest = value) %>% select(-c(lu.from, value))) %>%
  left_join(start_map_reproj %>% filter(lu.from == "pasture") %>% mutate(pasture = value) %>% select(-c(lu.from, value))) %>%
  left_join(start_map_reproj %>% filter(lu.from == "forest") %>% mutate(forest = value) %>% select(-c(lu.from, value))) %>% 
  left_join(start_map_reproj %>% filter(lu.from == "otherland") %>% mutate(otherland = value) %>% select(-c(lu.from, value)))

chaturvedi_restriction <- target_chaturvedi %>% 
  filter(restoration_area > 0)

colSums(chaturvedi_restriction %>% select(pasture))   /sum(colSums(chaturvedi_restriction %>% select(pasture, otherland, cropland)))
colSums(chaturvedi_restriction %>% select(otherland)) /sum(colSums(chaturvedi_restriction %>% select(pasture, otherland, cropland)))
colSums(chaturvedi_restriction %>% select(cropland)) /sum(colSums(chaturvedi_restriction %>% select(pasture, otherland, cropland)))

  
colSums(target_chaturvedi %>% select(restoration_area))

bastin_treecover <- rest_bastin.df %>% 
  left_join(ns_map) %>% 
  left_join(start_map_reproj %>% filter(lu.from == "cropland") %>% mutate(cropland = value) %>% select(-c(lu.from, value))) %>%
  left_join(start_map_reproj %>% filter(lu.from == "urban") %>% mutate(urban = value) %>% select(-c(lu.from, value))) %>%
  left_join(start_map_reproj %>% filter(lu.from == "newforest") %>% mutate(newforest = value) %>% select(-c(lu.from, value))) %>%
  left_join(start_map_reproj %>% filter(lu.from == "pasture") %>% mutate(pasture = value) %>% select(-c(lu.from, value))) %>%
  left_join(start_map_reproj %>% filter(lu.from == "forest") %>% mutate(forest = value) %>% select(-c(lu.from, value))) %>% 
  left_join(start_map_reproj %>% filter(lu.from == "otherland") %>% mutate(otherland = value) %>% select(-c(lu.from, value))) %>% 
  mutate(potential_area = pmax(0, rest_area/1000 - forest))   #value being the current forest cover %>% 
  
bastin_restriction <- bastin_treecover %>% 
  filter(potential_area >0)

colSums(bastin_treecover %>% select(potential_area))
colSums(bastin_treecover %>% filter(lu.from == "forest") %>% select(rest_area))

colSums(bastin_restriction %>% select(pasture)) /sum(colSums(bastin_restriction %>% select(pasture, otherland)))
colSums(bastin_restriction %>% select(otherland)) /sum(colSums(bastin_restriction %>% select(pasture, otherland)))

################################################################################
#################### Crossing with other maps ##################################
################################################################################


## cross otherland with biomass