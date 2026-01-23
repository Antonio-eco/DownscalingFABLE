##function IKI


remove_existingForest <- function(df, ns_map, start_map_reproj) {#rest_bastin.df, ns_map, start_map_reproj
  
  bastin_treecover <- df %>% 
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
  
  return(bastin_treecover)
}
