## This is the only script you need to run to generate your downscaling data and maps
## Author: Clara Douzal SDSN
##Craeted 25/10/2023




# Libraries ---------------------------------------------------------------

library(purrr)
library(here)

# Render Function ---------------------------------------------------------


render_report <- function(country, pathway) {
  
  template <- here("Program", "Downscaling_report.Rmd")
  outpath <- here("Output", country)
  if(!dir.exists(outpath)){
    dir.create(outpath)
  }
  out_file <- paste0(outpath, "/", gsub("-", "",Sys.Date()),"_" ,paste(country, pathway, "Downscaling", sep = "_"))
  
  parameters <- list(variable1 = country,
                     variable2 = pathway)
  
  rmarkdown::render(template,
                    output_file = out_file,
                    params = parameters)
  invisible(TRUE)
}


# Running over multiple parameters ---------------------------------------- 


params_list <- list(
  list("CHN"),
  list("CurrentTrends")) 

pmap(params_list, render_report)