
library(here)
library(scales)
library(ggnewscale)

date <- "251020"
country <- "IND"
stamp <-  date


country_start_areas <- read_rds(here("Output", country, paste0(stamp, "_country_start_areas.rds")))
X_long <- read_rds(here("Output", country, paste0(stamp, "_X_long.rds")))
pred_coeff_long <- read_rds(here("Output", country, paste0(stamp, "_pred_coeff_long.rds")))


FABLE_NC            <- read_rds(here("Output", country, paste0("251021", "_IND_NationalCommitments_HILDA_FABLE.rds")))
FABLE_NC_Bastin     <- read_rds(here("Output", country, paste0(stamp, "_IND_NC_Bastin_HILDA_FABLE.rds")))
FABLE_NC_Chaturvedi <- read_rds(here("Output", country, paste0("251021", "_IND_NC_Chaturvedi_HILDA_FABLE.rds")))

l_pathways <- c("FABLE_NC", "FABLE_NC_Bastin", "FABLE_NC_Chaturvedi")

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

### Map NC with Bonn challenge with no restriction -----

cum_diff_NC <- results_FABLE_NC
cum_diff_NC$out.res <- cum_difference(cum_diff_NC)

LUC_plot_2020_2050 <- LUC_plot(cum_diff_NC, rasterized_layer, label = "Cumulative LUC over 2020-2050")
plot_to_save <- LUC_plot_2020_2050$LUC.plot +
  ggplot2::scale_fill_gradient2(low = "red",
                                mid = "white",
                                high = "blue")+
  ggplot2::theme(legend.key.width = ggplot2::unit(0.8, "cm"),
                 panel.background = ggplot2::element_rect(fill = "lightgrey"))
plot_to_save      

ggplot2::ggsave(
  filename = here("Output", country, paste0(stamp, "2020_2050_cumulative_final_Land-use_change-map_NC.tiff")), 
  plot = plot_to_save, 
  units = "in", 
  height = 6, width = 10, dpi = 300)
rm(plot_to_save)



detailed_diff <- results_FABLE_NC
detailed_diff$out.res <- pairwise_net_difference(results_FABLE_NC, ns_map) 

Transition_plot_20_50 <- LUC_plot_cum_allLUC(
  detailed_diff, 
  rasterized_layer,
  label = "Transitions in % of pixel"
)

Transition_plot_20_50$LUC.plot + 
  theme(strip.text = element_text(size = 14),
        legend.title = element_text(size = 14))
# View(Transition_plot_20_50)
# View(Transition_plot_20_50$plot.df)
# View(detailed_diff$out.res)

ggplot2::ggsave(
  filename = here("Output", country, paste0(stamp, "2020_2050_cumulative_pair_wise_Land-use_change-map_NC_Bonn.tiff")), 
  plot = Transition_plot_20_50$LUC.plot + 
    theme(strip.text = element_text(size = 14),
          legend.title = element_text(size = 14)), 
  units = "in", 
  height = 6, width = 10, dpi = 300)
rm(plot_to_save)

## Restrictions------
 
## Bastin restriction -----

bastin_treecover <- rest_bastin.df %>% 
  left_join(ns_map) %>% 
  left_join(start_map_reproj %>% filter(lu.from == "cropland") 
            %>% mutate(cropland = value) %>% select(-c(lu.from, value))) %>%
  left_join(start_map_reproj %>% filter(lu.from == "urban") 
            %>% mutate(urban = value) %>% select(-c(lu.from, value))) %>%
  left_join(start_map_reproj %>% filter(lu.from == "newforest") 
            %>% mutate(newforest = value) %>% select(-c(lu.from, value))) %>%
  left_join(start_map_reproj %>% filter(lu.from == "pasture") 
            %>% mutate(pasture = value) %>% select(-c(lu.from, value))) %>%
  left_join(start_map_reproj %>% filter(lu.from == "forest") 
            %>% mutate(forest = value) %>% select(-c(lu.from, value))) %>% 
  left_join(start_map_reproj %>% filter(lu.from == "otherland") 
            %>% mutate(otherland = value) %>% select(-c(lu.from, value))) %>% 
  mutate(rest_area = rest_area/1000) %>% 
  mutate(potential_area_1 = rest_area - forest) %>% 
  mutate(potential_area = pmax(0,  rest_area - forest))   #value being the current forest cover %>% 

colSums(bastin_treecover %>% select(potential_area) %>% filter(potential_area >35))

bastin_restriction <- bastin_treecover %>% 
  filter(potential_area >35) %>% 
  select(ns)

rest_worthy     <- bastin_restriction
not_rest_worthy <- setdiff(ns_map["ns"], rest_worthy)

## plot to see where the cells are
test_notworthy <- ns_map %>% 
  filter(ns %in% rest_worthy[[1]]) %>% 
  mutate(value = 1) %>% 
  rbind(ns_map %>% 
          filter(ns %in% not_rest_worthy[[1]])
        %>% mutate(value = -1))%>%
  mutate(lu.from = "any",
         lu.to = "newforest",
         times = 2025) %>% 
  select(-ns) %>% 
  rename(ns = ns_int) %>% 
  select(ns, times, lu.from, lu.to, value)

results_DS_plot <- results_FABLE_NC
results_DS_plot$out.res <- test_notworthy

plot_bastin_restricted <- LUC_plot_restriction(results_DS_plot, rasterized_layer, 
                                        color = "PRGn", label = "Priority areas for restoration:")

plot_bastin_restricted$LUC.plot  

tiff(
  here("Output", country, paste0(tag, "_Bastin_restriction_map.tiff")),
  units = "in", height = 6, width = 10, res = 300
)
plot_bastin_restricted$LUC.plot
dev.off()

restriction_Bastin <- data.frame(
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

## Bastin potential ----

results_DS_plot <- results_FABLE_NC

results_DS_plot$out.res <- bastin_treecover %>% 
  select(-ns) %>% 
  rename(ns = ns_int) %>% 
  mutate(times = 2025,
         lu.from = "any",
         lu.to = "newforest",
         value = potential_area) %>% 
  select(ns, times, lu.from, lu.to, value)
  

plot_bastin_restoration <- LUC_plot(results_DS_plot, rasterized_layer, 
                                    year = 2025, LU = "newforest",
                                               color = "Blues", label = "Restoration tree cover potential (kha):") 
tiff(
  here("Output", country, paste0(tag, "_Bastin_potential_map.tiff")),
  units = "in", height = 6, width = 10, res = 300
)

plot_bastin_restoration$LUC.plot +
  ggplot2::theme(legend.text = element_text(size = 14))+
  ggplot2::theme(legend.title = element_text(size = 16))

dev.off()


## Chaturvedi restriction -----


target_chaturvedi <- rest_chaturvedi.df %>% 
  filter(restoration_type == "wide_scale_restoration") %>% 
  mutate(restoration_area = (p/100)*area/1000) %>% 
  left_join(ns_map)

colSums(target_chaturvedi %>% select(restoration_area))

chaturvedi_restriction <- target_chaturvedi %>% 
  filter(restoration_area > 0) %>% 
  select(ns)

rest_worthy     <- chaturvedi_restriction
not_rest_worthy <- setdiff(ns_map["ns"], rest_worthy)

## plot to see where the cells are
test_notworthy <- ns_map %>% 
  filter(ns %in% rest_worthy[[1]]) %>% 
  mutate(value = 1) %>% 
  rbind(ns_map %>% 
          filter(ns %in% not_rest_worthy[[1]])
        %>% mutate(value = -1))%>%
  mutate(lu.from = "any",
         lu.to = "newforest",
         times = 2025) %>% 
  select(-ns) %>% 
  rename(ns = ns_int) %>% 
  select(ns, times, lu.from, lu.to, value)

results_DS_plot <- results_FABLE_NC
results_DS_plot$out.res <- test_notworthy

plot_chaturvedi_restricted <- LUC_plot_restriction(results_DS_plot, rasterized_layer, 
                                               color = "PRGn", label = "Priority areas for restoration:")

plot_chaturvedi_restricted$LUC.plot  

tiff(
  here("Output", country, paste0(tag, "_Chaturvedi_restriction_map.tiff")),
  units = "in", height = 6, width = 10, res = 300
)
plot_chaturvedi_restricted$LUC.plot
dev.off()

restriction_Chaturvedi <- data.frame(
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

## Chaturvedi potential ----

results_DS_plot <- results_FABLE_NC
results_DS_plot$out.res <- ns_map %>% 
  left_join(rest_chaturvedi.df %>% filter(restoration_type == "wide_scale_restoration")) %>% 
  select(-ns) %>% 
  rename(ns = ns_int) %>% 
  mutate(times = 2025,
         lu.from = "any",
         lu.to = "newforest",
         value = area) %>% 
  select(ns, times, lu.from, lu.to, value)

results_DS_plot$out.res <- target_chaturvedi %>% 
  select(-ns) %>% 
  rename(ns = ns_int) %>% 
  mutate(times = 2025,
         lu.from = "any",
         lu.to = "newforest",
         value = restoration_area) %>% 
  select(ns, times, lu.from, lu.to, value)

plot_chaturvedi_restoration <- LUC_plot(results_DS_plot, rasterized_layer, 
                                    year = 2025, LU = "newforest",
                                    color = "Blues", label = "Restoration potential (kha):") 
tiff(
  here("Output", country, paste0(tag, "_Chaturvedi_potential_map.tiff")),
  units = "in", height = 6, width = 10, res = 300
)

plot_chaturvedi_restoration$LUC.plot +
  ggplot2::theme(legend.text = element_text(size = 14))+
  ggplot2::theme(legend.title = element_text(size = 16))

dev.off()

### Map NC Bastin and NC Bastin restriction

results_FABLE_NC_Bastin_restriction <- downscale(
  targets      = temp,
  start.areas  = country_start_areas,
  xmat         = X_long,
  betas        = pred_coeff_long,
  restrictions = restriction_Bastin
)

results_FABLE_NC_Bastin_restriction$out.res <- format_to_id_map(results_FABLE_NC_Bastin_restriction, ns_map)

results_several_pathways <- list(results_FABLE_NC, results_FABLE_NC_Bastin, results_FABLE_NC_Bastin_restriction)
pathways_names <- c(  "Bonn Challenge", "Bastin", "Bastin & restrictions")

LUC_plot_pathways_bastin <- LUC_plot_compare_pathways(results_several_pathways, rasterized_layer, grid, pathways_names, label = "Area in 1000 ha per pixel")
LUC_plot_pathways_bastin$LUC.plot

ggplot2::ggsave(
  filename = here("Output", country, paste0(tag, "2020_2050_per_pathway_cumulative_final_LUC-map_Bastin.tiff")), 
  plot = LUC_plot_pathways_bastin$LUC.plot+ 
    theme(strip.text = element_text(size = 12),
          legend.title = element_text(size = 12)), , 
  units = "in", 
  height = 5.5, width = 10, dpi = 300)

### Map NC Chaturvedi and NC Chaturvedi restriction

results_FABLE_NC_Chaturvedi_restriction <- downscale(
  targets      = temp,
  start.areas  = country_start_areas,
  xmat         = X_long,
  betas        = pred_coeff_long,
  restrictions = restriction_Chaturvedi
)

results_FABLE_NC_Chaturvedi_restriction$out.res <- format_to_id_map(results_FABLE_NC_Chaturvedi_restriction, ns_map)

results_several_pathways <- list(results_FABLE_NC, results_FABLE_NC_Chaturvedi, results_FABLE_NC_Chaturvedi_restriction)
pathways_names <- c(  "Bonn Challenge", "Chaturvedi", "Chaturvedi & restrictions")

LUC_plot_pathways_Chaturvedi <- LUC_plot_compare_pathways(results_several_pathways, rasterized_layer, grid, pathways_names, label = "Area in 1000 ha per pixel")
LUC_plot_pathways_Chaturvedi$LUC.plot

ggplot2::ggsave(
  filename = here("Output", country, paste0(tag, "2020_2050_per_pathway_cumulative_final_LUC-map_Chaturvedi.tiff")), 
  plot = LUC_plot_pathways_Chaturvedi$LUC.plot+ 
    theme(strip.text = element_text(size = 12),
          legend.title = element_text(size = 12)), , 
  units = "in", 
  height = 5.5, width = 10, dpi = 300)


### map bastin pariwise ----

detailed_diff_rest <- results_FABLE_NC_Bastin_restriction
detailed_diff_rest$out.res <- pairwise_net_difference(results_FABLE_NC_Bastin_restriction, ns_map)

Transition_plot_20_50_rest <- LUC_plot_cum_allLUC(
  detailed_diff_rest, 
  rasterized_layer,
  label = "Transitions in % of pixel"
)

Transition_plot_20_50_rest$LUC.plot

ggplot2::ggsave(
  filename = here("Output", country, paste0(tag, "2020_2050_restricted_cumulative_pair_wise_Land-use_change-map.tiff")), 
  plot = Transition_plot_20_50_rest$LUC.plot+ 
    theme(strip.text = element_text(size = 12),
          legend.title = element_text(size = 12)), 
  units = "in", 
  height = 6, width = 10, dpi = 300)
rm(plot_to_save)