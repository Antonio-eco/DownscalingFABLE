res_restriction <- results_DS_rest_full

cells_notworthy <- ns_map %>% 
  filter(ns %in% unique(c(not_restoration_worthy_griscom[[1]],
                          not_restoration_worthy_walker[[1]]))) %>% 
  select(ns_int) %>% 
  rename(ns = ns_int) %>% 
  mutate(worthy = 0)

cells_worthy <- ns_map %>% 
  filter(!(ns %in% unique(c(not_restoration_worthy_griscom[[1]],
                          not_restoration_worthy_walker[[1]])))) %>% 
  select(ns_int) %>% 
  rename(ns = ns_int) %>% 
  mutate(worthy = 1)


special_grid <- grid %>%  
  select(id_c, area) %>% 
  left_join(ns_map) %>% 
  select(ns_int, area) %>% 
  rename(ns = ns_int)

res_restriction$out.res <- res_restriction$out.res %>% 
  left_join(rbind(cells_notworthy, cells_worthy)) %>% 
  left_join(special_grid)

downscaled_est <- res_restriction$out.res

test <- downscaled_est %>% 
  filter(lu.to == "newforest") %>% 
  mutate( restoreButShouldnot = ifelse(value >0 & worthy == 0, 1, 0))

area_rest <- start_map_reproj %>% 
  left_join(ns_map) %>% 
  select(ns_int, lu.from, value) %>% 
  rename(ns = ns_int) %>% 
  left_join(rbind(cells_notworthy, cells_worthy)) %>% 
  filter(worthy == 1) %>% 
  filter(lu.from != "forest") %>% 
  filter(lu.from != "urban")

colSums(area_rest["value"])


area_rest2 <- start_map_reproj %>% 
  left_join(ns_map) %>% 
  select(ns_int, lu.from, value) %>% 
  rename(ns = ns_int) %>% 
  left_join(rbind(cells_notworthy, cells_worthy)) %>% 
  group_by(ns) %>% 
  summarise(value = sum(value))%>% 
  left_join(special_grid) %>% 
  mutate(area = area /1000)

colSums(area_rest)


colSums(area_rest2)

area_rest3 <- start_map_reproj %>% 
  left_join(ns_map) %>% 
  select(ns_int, lu.from, value) %>% 
  rename(ns = ns_int) %>% 
  left_join(rbind(cells_notworthy, cells_worthy)) %>% 
  filter(worthy == 1) %>% 
  filter(lu.from != "forest") %>% 
  filter(lu.from != "urban") %>% 
  pivot_wider(names_from = lu.from, values_from = value)

colSums(area_rest3)

area_target_restoration <- FABLE %>% 
  filter(lu.to == "newforest") %>% 
  pivot_wider(names_from = lu.from, values_from = value)

colSums(area_target_restoration %>% select(-c(times, lu.to)))

start_map

area_rest3 <- start_map %>% 
  left_join(ns_map) %>% 
  select(ns_int, lu.from, value) %>% 
  rename(ns = ns_int) %>% 
  left_join(rbind(cells_notworthy, cells_worthy)) %>% 
  group_by(ns) %>% 
  summarise(value = sum(value))%>% 
  left_join(special_grid) %>% 
  mutate(area = area /1000)

colSums(area_rest3)

temp_restriction <- restriction %>% mutate(transition = paste0(lu.from, "_", lu.to))
temp_targets <- FABLE %>% filter(times >= 2020) %>% mutate(transition = paste0(lu.from, "_", lu.to))
setdiff(temp_restriction$transition, temp_targets$transition)
character(0)
setdiff(temp_targets$transition, temp_restriction$transition)
