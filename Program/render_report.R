# Clara Douzal SDSN

# Master file for generating automatically individual report for downscaling


# latest version of downscalR ---------------------------------------------

# devtools::install_github("tkrisztin/downscalr", ref="HEAD", repos = "http://cran.us.r-project.org")


# Libraries ---------------------------------------------------------------

library(purrr)
library(here)

# Render Function ---------------------------------------------------------


render_report <- function(country, pathway) {
  
  # country = "ARG"
  # pathway = "CurrentTrends"
  
  template <- here("Program", "Downscaling_report.Rmd")
  outpath <- here("Output", country)
  if(!dir.exists(outpath)){
    dir.create(outpath)
  }
  out_file <- paste0(outpath, "/", format(Sys.Date(),format = "%y%m%d"), "_" , pathway, "_", country, "_Dowscaling_report")
  
  parameters <- list(variable1 = country,
                     variable2 = pathway)
  
  rmarkdown::render(template,
                    output_file = out_file,
                    params = parameters)
  invisible(TRUE)
}


# Running over multiple parameters ---------------------------------------- 

params_list_ct <- list(list(
  "IND"),
  list("CurrentTrends")) 

pmap(params_list_ct, render_report)

