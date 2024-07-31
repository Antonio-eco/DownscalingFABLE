#install.packages("tidyr","tibble","dplyr")
library(tidyr)
library(tibble)
library(dplyr)

# just install the downscalr package
##install.packages("devtools")
#devtools::install_github("tkrisztin/downscalr", ref="HEAD", repos = "http://cran.us.r-project.org")

# install the downscalr package with vignette 
# (note: this may take a couple minutes since some computing is required to render the vignette). 
#install.packages(c("devtools", "knitr", "rmarkdown"))
#devtools::install_github("tkrisztin/downscalr", ref="HEAD", build_vignettes = TRUE, repos = "http://cran.us.r-project.org")


#devtools::install_github("tkrisztin/downscalr",ref="HEAD")
library(downscalr)

##Preliminaries

# Set example land use change origin to Cropland
example_LU_from <- "Cropland"

# load data
data(argentina_df,argentina_luc)

# Prepare data for econometric Multinomial logistic (MNL) regression

#' Data of land-use changes in Argentina between 2000 and 2010, as well as 2010 and 2020.
#' 
#' @format A data frame with 277632 rows and 5 variables:
#' \describe{
#'   \item{ns}{Pixel ids for Argentina}
#'   \item{lu.from}{Origin land-use class of converion}
#'   \item{lu.to}{Destination land-use class to which lu.from is converted}
#'   \item{Ts}{Time-period; either 2000 (changes 2000 to 2010) or 2010 (changes 2010 to 2020)}
#'   \item{value}{Percent of lu.from in pixel being converted}
Yraw <- argentina_luc

Xraw <- argentina_df$xmat

Y <- dplyr::filter(Yraw,lu.from == example_LU_from & Ts == 2000) %>% #changes 2000 to 2010
  tidyr::pivot_wider(names_from = lu.to) %>%
  tibble::column_to_rownames(var = "ns") %>%
  dplyr::select(-c(lu.from, Ts))

X <- Xraw %>% tidyr::pivot_wider(names_from = "ks") %>%
  dplyr::arrange(match(ns,Y$ns)) %>%
  tibble::column_to_rownames(var = "ns") #%>% 
  # mutate(SumLandCoverShare = Cropland + Forest + OtherLand + Pasture + Plantations + Urban) %>% 
  # mutate(SumLandCoverShare = ifelse(SumLandCoverShare == 0, 1, SumLandCoverShare)) %>% 
  # mutate(Cropland = Cropland/SumLandCoverShare,
  #        Forest = Forest/SumLandCoverShare, 
  #        OtherLand = OtherLand/SumLandCoverShare, 
  #        Pasture = Pasture/SumLandCoverShare, 
  #        Plantations = Plantations/SumLandCoverShare, 
  #        Urban = Urban/SumLandCoverShare)

# Xcheck <- X %>% 
#   mutate(SumLandCoverShare = Cropland + Forest + OtherLand + Pasture + Plantations + Urban) %>% 
#   select(Cropland, Forest, OtherLand, Pasture, Plantations, Urban, SumLandCoverShare) %>% 
#   arrange(SumLandCoverShare)
# Xcheck$ns <- row.names(Xcheck)
# Xcheck$lu.from <- "test"
# Xcheck$lu.to <- "test"
# Xcheck$Ts <- 2000
# Xcheck$value <- Xcheck$SumLandCoverShare



## Estimate MNL model

#Compute MNL model with standard settings
baseline <- which(colnames(Y) == example_LU_from)

# Decreased niter = 100 and nburn = 50 for faster computation 
results_MNL <- mnlogit(as.matrix(X), as.matrix(Y), baseline, niter = 100, nburn = 50, A0 = 10^4)
# Strongly recommended if you have more time: 
# results_MNL <- mnlogit(as.matrix(X), as.matrix(Y), baseline, niter = 2000, nburn = 1000, A0 = 10^4)

## Downscale module

### Preliminaries

# Set example year of downscaling target. e.g. 2010
example_time <- 2050

# Prepare inputs for downscale function
## get target data and filter for year and lu.from according to example, i.e. 2010 & cropland
data(argentina_FABLE)

arg_targets_crop_2010 <- dplyr::filter(argentina_FABLE, lu.from == example_LU_from &
                                         times == example_time)

## pivot explanatory input data into long format
X_long <-  X %>% tibble::rownames_to_column(var="ns") %>%
  tidyr::pivot_longer(!ns, names_to="ks", values_to = "value")

## prior coefficients: compute mean of posterior draws & pivot to long format  
pred_coeff <- as.data.frame(apply(results_MNL$postb[,-baseline,],c(1,2),mean))

pred_coeff_long <- pred_coeff %>% tibble::rownames_to_column(var="ks") %>%
  tidyr::pivot_longer(!ks, names_to = "lu.to", values_to="value") %>% 
  dplyr::arrange(lu.to) %>% 
  tibble::add_column(lu.from=example_LU_from, .after="ks")

## starting area should provide information on area of each (higher resolution) grid cell, i.e ns
arg_start_areas_crop <- dplyr::filter(argentina_df$lu_levels, lu.from == example_LU_from)

## Downscaling computation and results

results_DS <- downscale(targets = arg_targets_crop_2010, 
                        start.areas =  arg_start_areas_crop, xmat = X_long, betas = pred_coeff_long)

downscaled_LUC <- results_DS$out.res

## Visualise results

#load raster data for argentina
#data(argentina_raster)
argentina_raster = terra::rast(system.file("extdata", "argentina_raster.grd", package = "downscalr"))

# call plot function & plot downscaled LUC projections from Cropland for the year 2010
LUC_argentina_2010_Crop_plot <- LUC_plot(results_DS, argentina_raster)
LUC_argentina_2010_Crop_plot$LUC.plot

results_DS$out.res <- Xcheck %>% select(colnames(results_DS$out.res))

LUC_plot(Xcheck %>% select(colnames(results_DS)), argentina_raster)
