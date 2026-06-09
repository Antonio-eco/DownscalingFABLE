# DownscalingFABLE/R/run_country.R

#renv::restore()

remotes::install_github("FABLE-consortium/FABLEDownscalR", dependencies = TRUE)

# (1) Load your package (devtools::load_all is done by .Rprofile)
library(FABLEDownscalR)
library(dplyr)
library(here)


# (2) Read config (you can implement fdr_read_yaml in the package)
cfg <- fdr_read_config(here::here("config", "IND.yml"))
cfg$stamp <- format(Sys.Date(), "%y%m%d")
set.seed(cfg$seed)

# (3) Load raw inputs (geojson + mapping + grid + FABLE)
inputs <- fdr_load_inputs(
  data_root        = cfg$data_root,
  country          = cfg$country,
  start_map_source = cfg$start_map_source,
  pathway          = cfg$pathway
)

# (4) Build land-cover change calibration table (DownscalR format)
luc <- lc_build_country_luc(
  LandCoverChange_df = inputs$spatial$landcoverchange,
  map_HILDA_LUC      = inputs$mapping$map_LUC,
  Ts                 = 2015 # year and expected_lu could be changed accordingly to land use map
)

# (5) Harmonise start map to match FABLE baseline totals
harm <- fdr_harmonize_start_map(
  LandCover_df  = inputs$spatial$landcover,
  map_LC        = inputs$mapping$map_LC,
  LC_targets    = inputs$LC_targets
)

# (6) Build ns_map + rasterized ID layer (resolution controlled by YAML)
id <- fdr_build_id_maps(
  grid_sp         = inputs$grid_sp,      # e.g. Travel as sp/sf geometry for cells
  ns_map          = harm$ns_map,
  pixel_res_m     = cfg$pixel_res_m
)

# (7) Build priors (X matrix etc.) + drop cells with NA covariates
priors <- fdr_build_priors(
  inputs         = inputs,
  start_map      = harm$start_map_reproj,
  good_ns_only   = TRUE
)

# (8) Fit MNL + downscale
results <- fdr_run_downscaling( # not clear if restrictions are applied as the call in results <- downscalr::downscale is missing restrictions = NULL
  targets      = fdr_wrangle_fable_targets(inputs$FABLE_targets, min_year = 2020), # why is min_year set at 2020 if FABLE data goes back to 2000
  country_luc  = luc$country_luc,
  priors       = priors,
  mnl_niter    = cfg$mnl_niter,
  mnl_nburn    = cfg$mnl_nburn
)

# (9) Save outputs (consistent naming)
fdr_save_outputs(
  country  = cfg$country,
  tag      = fdr_make_tag(cfg),
  outputs  = list(
    start_map_reproj   = harm$start_map_reproj,
    ns_map             = harm$ns_map,
    rasterized_layer   = id$rasterized_layer,
    grid_sf            = id$grid_sf,
    X_long             = results$X_long,
    betas              = results$pred_coeff_long,
    country_start_areas= results$country_start_areas,
    downscaled_LUC     = results$downscaled_LUC
  )
)

message("✅ Done: ", cfg$country, " (", cfg$pathway, ")")
# 
# p <- fdr_plot_downscaled_maps(
#   out_res = results$out.res,                 # or res$downscaled_LUC if you keep that
#   rasterized_layer = id$rasterized_layer,
#   ns_map = id$ns_map              
# )
# 
# print(p)

p_main_LU <- fdr_plot_downscaled_LU_one(
  out_res          = results$out.res, 
  rasterized_layer = id$rasterized_layer,
  ns_map           = id$ns_map)

ggplot2::ggsave(
  filename = here("Output", cfg$country, paste0("MainLU_", cfg$pathway, ".tiff")), 
  plot = p_main_LU, 
  units = "in", 
  height = 10, width = 10, dpi = 300)

p_LU <- fdr_plot_downscaled_LU(
  out_res          = results$out.res, 
  rasterized_layer = id$rasterized_layer,
  ns_map           = id$ns_map
)

ggplot2::ggsave(
  filename = here("Output", cfg$country, paste0("LU_", cfg$pathway, ".tiff")), 
  plot = p_LU, 
  units = "in", 
  height = 10, width = 10, dpi = 300)

p_LUC <- fdr_plot_downscaled_LUC(
  out_res          = results$out.res, 
  rasterized_layer = id$rasterized_layer,
  ns_map           = id$ns_map
)

ggplot2::ggsave(
  filename = here("Output", cfg$country, paste0("LUC_", cfg$pathway, ".tiff")), 
  plot = p_LUC, 
  units = "in", 
  height = 10, width = 10, dpi = 300)
