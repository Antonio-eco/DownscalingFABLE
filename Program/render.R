###has not been used, written for later

library(here)
library(rmarkdown)

country <- "IND"
pathway <- "GlobalSustainability"

out_dir  <- here("Output", country)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

stamp    <- format(Sys.time(), "%y%m%d_%H%M%S")
out_file <- sprintf("%s_%s_%s_Downscaling_report", stamp, country, pathway)

rmarkdown::render(here("Program", "Downscaling_report.Rmd"),
       params = list(
         country   = country,
         pathway   = pathway,
         mnl_niter = 100,
         mnl_nburn = 50,
         seed      = 1212
       ),
       output_file = out_file,
       output_dir  = out_dir,
       clean = TRUE,
       envir = new.env(parent = globalenv())
)
