# function to plot for one pathway, land use transitions in proportion of the total pixel from 2020 to 2050

LUC_plot_cum_allLUC <- function(res, rasterfile, LU.from=NULL, LU=NULL, color = "Greens", label = "Area in ha per pixel"){
  
  plot_df <- terra::as.data.frame(rasterfile, xy = TRUE, na.rm = FALSE)%>%
    #to avoid that a cell from my grid is attributed to several x and y pairs
    dplyr::distinct(layer, .keep_all = TRUE)
  
  to.plot <- res$out.res
  
  colnames(plot_df) <- c("x", "y", "ns")
  
  ns = lu.to = lu.from = value = x = y= NULL
  if(is.null(LU.from) & is.null(LU)){
    inputs <- to.plot %>% dplyr::group_by(ns, lu.to, lu.from)
  } else if(!(is.null(LU.from) | is.null(LU))){
    inputs <- to.plot %>% dplyr::group_by(ns, lu.to, lu.from) %>% dplyr::summarise(value = sum(value),.groups = "keep") %>% subset(lu.to==LU & lu.from==LU.from)
  } else if(is.null(LU.from)){
    inputs <- to.plot %>% dplyr::group_by(ns, lu.to, lu.from) %>% dplyr::summarise(value = sum(value),.groups = "keep") %>% subset(lu.to==LU)
  } else {
    inputs <- to.plot %>% dplyr::group_by(ns, lu.to, lu.from) %>% dplyr::summarise(value = sum(value),.groups = "keep") %>% subset(lu.from==LU.from)
  }
  
  plot_df <- merge(plot_df, inputs, by="ns") %>%
    mutate(lu.from = factor(lu.from, 
                            levels = c("cropland", "otherland", "pasture", "forest")),
           lu.to   = factor(lu.to,   
                            levels = c("cropland", "otherland", "pasture", "newforest", "urban")))
  
  row_order <- c()
  
  df2 <- df 
  
  
  plot_obj <- ggplot2::ggplot() +
    ggplot2::geom_tile(data=plot_df, aes(x=x, y=y, fill=value, group=lu.to), alpha=0.8) +
    ggplot2::scale_fill_distiller(palette = color , name=label, direction = 1)+
    ggplot2::scale_fill_gradient2(low = "red",
                                  mid = "white",
                                  high = "blue",
                                  name=label)+
    ggplot2::coord_equal() +
    ggthemes::theme_map() +
    ggplot2::theme(legend.position="bottom") +
    ggplot2::theme(legend.key.width=unit(2, "cm"))+
    ggplot2::theme(panel.background = element_rect(fill = "lightgrey"))+
    if(is.null(LU.from) & is.null(LU)){
      ggplot2::facet_grid(lu.from~lu.to, switch = "y")
    } else if(!(is.null(LU.from) | is.null(LU))){
      ;
    } else if(is.null(LU.from)){
      ggplot2::facet_wrap(~lu.from)
    } else {
      ggplot2::facet_wrap(~lu.to)
    }
  
  return(list(LUC.plot=plot_obj, plot.df=plot_df))
}
