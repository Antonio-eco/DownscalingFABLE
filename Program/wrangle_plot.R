#wrangle for plot

#Make the results use the same integer codes for plotting as standardized raster

format_to_id_map <- function(res, ns_map){
  
  temp <- res$out.res
  
  formated <- temp %>%
    dplyr::mutate(ns = as.character(ns)) %>%
    dplyr::left_join(ns_map %>% dplyr::select(ns, ns_int), by = "ns") %>%
    dplyr::mutate(ns = ns_int) %>%
    dplyr::select(-ns_int)
  
  return(formated)
  
}


## Create 'gain' and 'loss' datasets out of the downscale results
loss_and_gains <- function(res){
  
  temp <- res$out.res
  
  to.plot <- temp %>%
    filter(lu.from != lu.to) %>%
    # Group by `lu.to`, `ns`, and `times` for gains
    group_by(lu.to, ns, times) %>%
    summarise(gain = sum(value, na.rm = TRUE), .groups = "drop") %>%
    mutate(Type = "gain", lu = lu.to) %>%
    
    # Bind the losses dataset by grouping by `lu.from`
    bind_rows(
      temp %>%
        filter(lu.from != lu.to) %>%
        group_by(lu.from, ns, times) %>%
        summarise(loss = -sum(value, na.rm = TRUE), .groups = "drop") %>%
        mutate(Type = "loss", lu = lu.from)
    ) %>%
    
    # Summarize the final dataset grouping by land-use (`lu`), time, and pixel (`ns`)
    group_by(lu, times, ns) %>%
    summarise(value = sum(gain, loss, na.rm = TRUE), .groups = "drop") %>%
    rename(lu.to = lu)
 
return(to.plot)  
   
}

## Compute cumulative difference between 2020 and 2050 at le final land cover level

cum_difference <-  function(res){
  
  temp <- res$out.res
  
  to.plot <- temp %>% 
    #keep only 2020 and 2050
    filter(times %in% c(2020, 2050)) %>% 
    group_by(lu.to, ns, times) %>% 
    #aggregate for each year, each cells and each final land use the total land coverage at the end of the period
    #we do not keep the information on which land cover is being substituted
    summarise(value = sum(value)) %>% 
    pivot_wider(names_from = times, values_from = value) %>% 
    #compute difference in land cover in each cells 
    mutate(value= `2050` - `2020`) %>% 
    mutate(times = 2050)
  
  return(to.plot)
  
}


## compute pair-wise net transitions over 2020-2050
pairwise_net_difference <- function(res, ns_map){
  
  temp <- res$out.res
  
  to.plot <- temp %>%
    #do not keep 
    filter(lu.from != lu.to) %>%
    #Positive = net flow from the alphabetically-first class (a) to the second (b); 
    #negative = net the other way.
    #will be used to compute the net flow from a to b
    mutate(
      a = pmin(lu.from, lu.to),
      b = pmax(lu.from, lu.to),
      sign = if_else(lu.from == a & lu.to == b, 1, -1)
    ) %>% 
    #adding the missing b to a flow 
    rbind(
      temp %>%
        filter(lu.from != lu.to) %>%
        mutate(
          b = pmin(lu.from, lu.to),
          a = pmax(lu.from, lu.to),
          sign = if_else(lu.from == a & lu.to == b, 1, -1)
        )
    ) %>%
    #compute the net flow for each pair
    group_by(ns, a, b) %>%
    summarise(net = sum(sign * value, na.rm = TRUE), .groups = "drop") %>%
    arrange( ns, a, b) %>% 
    rename(lu.from = a, 
           lu.to   = b) %>% 
    left_join(ns_map %>% select(-ns), by = c("ns" = "ns_int")) %>% 
    left_join(grid %>% select(id_c, area)) %>% 
    mutate(value = 100 * net/(area/1000)) %>% 
    select(ns, lu.from, lu.to, value) %>% 
    filter(lu.from != "urban") %>% 
    filter(lu.to != "forest") %>% 
    filter(lu.from != "newforest") 
  
  return(to.plot)
}


