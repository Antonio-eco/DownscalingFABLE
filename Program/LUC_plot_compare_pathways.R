# function to plot for several pathways, land use transitions in proportion of the total pixel from 2020 to 2050

LUC_plot_compare_pathways <- function(list_res, rasterfile, grid, names, Path=NULL, LU=NULL, color = "Greens", label = "Area in ha per pixel"){
  
  plot_df <- terra::as.data.frame(rasterfile, xy = TRUE, na.rm = FALSE) %>%
    #to avoid that a cell from my grid is attributed to several x and y pairs
    dplyr::distinct(layer, .keep_all = TRUE)
  
  # plot_df <- grid %>% 
  #   select(x, y, id_c) %>% 
  #   merge(ns_map %>% select(-ns) %>% mutate(ns_int = as.numeric(ns_int))) %>% 
  #   rename(ns = ns_int) %>% 
  #   select(-id_c)
  
  to.plot <- data.frame()
  for(cur in 1:length(list_res)){
    
    temp <- list_res[[cur]]$out.res %>% 
      group_by(lu.to, lu.from, ns, times) %>% 
      summarise(value = sum(value)) %>% 
      filter(times %in% c(2025, 2050)) %>% 
      pivot_wider(names_from = times, values_from = value) %>% 
      mutate(value= `2050` - `2025`) %>% 
      mutate(times = 2050) %>% 
      mutate(pathway = names[cur])
    
    to.plot <- rbind.data.frame(to.plot, temp )
  }
  
  
  colnames(plot_df) <- c("x", "y", "ns")
  
  ns = lu.to = pathway = value = x = y= NULL
  if(is.null(Path) & is.null(LU)){
    inputs <- to.plot %>% dplyr::group_by(ns, lu.to, pathway) %>% dplyr::summarise(value = sum(value),.groups = "keep")
  } else if(!(is.null(Path) | is.null(LU))){
    inputs <- to.plot %>% dplyr::group_by(ns, lu.to, pathway) %>% dplyr::summarise(value = sum(value),.groups = "keep") %>% subset(lu.to==LU & pathway==Path)
  } else if(is.null(Path)){
    inputs <- to.plot %>% dplyr::group_by(ns, lu.to, pathway) %>% dplyr::summarise(value = sum(value),.groups = "keep") %>% subset(lu.to==LU)
  } else {
    inputs <- to.plot %>% dplyr::group_by(ns, lu.to, pathway) %>% dplyr::summarise(value = sum(value),.groups = "keep") %>% subset(pathway==Path)
  }
  
  plot_df <- merge(plot_df, inputs, by="ns")
  
 
  
  plot_obj <- ggplot2::ggplot() +
    ggplot2::geom_tile(data=plot_df, aes(x=x, y=y, fill=value, group=lu.to), alpha=0.8) +
    #ggplot2::scale_fill_distiller(palette = color , name=label, direction = 1)+
    ggplot2::scale_fill_gradient2(low = "red",
                                  mid = "white",
                                  high = "blue",
                                  name=label)+
    ggplot2::coord_equal() +
    ggthemes::theme_map() +
    ggplot2::theme(legend.position="bottom") +
    ggplot2::theme(legend.key.width=unit(2, "cm"))+
    ggplot2::theme(panel.background = element_rect(fill = "lightgrey"))+
    if(is.null(Path) & is.null(LU)){
      ggplot2::facet_grid(pathway~lu.to, switch = "y")
    } else if(!(is.null(Path) | is.null(LU))){
      ;
    } else if(is.null(Path)){
      ggplot2::facet_wrap(~pathway)
    } else {
      ggplot2::facet_wrap(~lu.to)
    }
  
  return(list(LUC.plot=plot_obj, plot.df=plot_df))
}
