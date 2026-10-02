# DownscalingFABLE/R/run_country.R

# ==============================================================================
# downscalr_patch.R
#
# Runtime monkey-patch for a bug in downscalr::solve_biascorr.mnl() (R/solve_biascorr.R).
#
# Bug: 'error_restrictions = restr.mat[,error_ind]' drops to a 1D vector when
# exactly ONE (lu.from,lu.to) target triggers the bias-correction fallback
# (i.e. sum(error_ind) == 1). The subsequent 'error_restrictions[,ccc]' then
# fails with: 'Error in error_restrictions[, ccc] : incorrect number of
# dimensions'. Fix: add drop = FALSE so it stays a (1-column) matrix.
#
# This block redefines solve_biascorr.mnl with that one-line fix and installs
# it into the downscalr namespace via assignInNamespace(). Must run AFTER
# library(downscalr). No package reinstall or separate file needed.
#
# Upstream fix: change line ~207 of R/solve_biascorr.R in tkrisztin/downscalr
# from 'restr.mat[,error_ind]' to 'restr.mat[,error_ind,drop=FALSE]'.
# ==============================================================================

solve_biascorr.mnl = function(targets,areas,xmat,betas,priors = NULL,restrictions=NULL,
                              options = downscale_control()) {
  lu.from <- unique(targets$lu.from)
  lu.to <- unique(targets$lu.to)
  ks = unique(betas$ks)
  
  out.solver <- list()
  curr.lu.from <- lu.from[1]
  full.out.res = NULL
  for(curr.lu.from in lu.from){
    err.txt = paste0(curr.lu.from," ",options$err.txt)
    
    # Extract targets
    curr.targets = dplyr::filter(targets,lu.from == curr.lu.from)$value
    names(curr.targets) <- targets$lu.to[targets$lu.from == curr.lu.from]
    curr.lu.to = names(curr.targets)
    
    # Extract betas
    curr.betas = dplyr::filter(betas,lu.from == curr.lu.from & lu.to %in% curr.lu.to) %>%
      tidyr::pivot_wider(names_from = "lu.to",values_from = "value",id_cols = "ks") %>%
      tibble::column_to_rownames(var = "ks")
    curr.betas = as.matrix(curr.betas)
    
    # Extract xmat
    ## IMPORTANT CHECK ORDER OF VARIABLES FIXED
    curr.xmat = dplyr::filter(xmat,ks %in% row.names(curr.betas)) %>%
      tidyr::pivot_wider(names_from = "ks",values_from = "value",id_cols = "ns")  %>%
      tibble::column_to_rownames(var = "ns") %>%
      dplyr::select(rownames(curr.betas))
    curr.xmat = as.matrix(curr.xmat)
    
    # Extract areas
    curr.areas = dplyr::filter(areas,lu.from == curr.lu.from)$value
    names(curr.areas) <- areas$ns[areas$lu.from == curr.lu.from]
    # BUGFIX: MW, order of curr.areas wrong, need to re-arrange based on xmat
    if (nrow(curr.xmat) > 0) {
      curr.areas = curr.areas[match(rownames(curr.xmat),names(curr.areas))]
    }
    
    # Extract priors
    ## IMPORTANT CHECK ORDER OF VARIABLES SIMILARLY TO XMAT
    if (!is.null(priors) && any(priors$lu.from == curr.lu.from)) {
      curr.priors = dplyr::filter(priors,lu.from == curr.lu.from & lu.to %in% curr.lu.to) %>%
        tidyr::pivot_wider(names_from = lu.to,values_from = "value",id_cols = "ns") %>%
        tibble::column_to_rownames(var = "ns")
      curr.prior_weights = dplyr::filter(priors,lu.from == curr.lu.from & lu.to %in% curr.lu.to) %>%
        tidyr::pivot_wider(names_from = lu.to,values_from = "weight",id_cols = "ns") %>%
        tibble::column_to_rownames(var = "ns")
      name_match = match(names(curr.areas),row.names(curr.priors))
      curr.priors = curr.priors[name_match,,drop = FALSE]
      curr.prior_weights = curr.prior_weights[name_match,,drop = FALSE]
      # check if betas have been provided for priors already
      mixed_priors =c()
      nonmixed_priors = colnames(curr.priors)
      if (any(colnames(curr.betas) %in% colnames(curr.priors))) {
        mixed_priors = colnames(curr.betas)[which(colnames(curr.betas) %in% colnames(curr.priors))]
        nonmixed_priors = nonmixed_priors[!nonmixed_priors %in% mixed_priors]
        #warning(paste0(err.txt,
        #               "Priors provided for lu.from/lu.to combinations for which betas exist.\n These will be overwritten."))
        #curr.betas = curr.betas[,-which(colnames(curr.betas) %in% colnames(curr.priors)),drop = FALSE]
      }
      curr.priors = as.matrix(curr.priors)
    } else {curr.priors = NULL}
    
    # Extract restrictions
    if (!is.null(restrictions) && any(restrictions$lu.from == curr.lu.from)) {
      curr.restrictions = dplyr::filter(restrictions,lu.from == curr.lu.from) %>%
        tidyr::pivot_wider(names_from = lu.to,values_from = "value",id_cols = "ns") %>%
        tibble::column_to_rownames(var = "ns")
      curr.restrictions = curr.restrictions[match(names(curr.areas),row.names(curr.restrictions)),,drop = FALSE]
      curr.restrictions = as.matrix(curr.restrictions)
    } else {curr.restrictions = NULL}
    
    p = length(curr.targets)
    n = length(curr.areas)
    p1 = ncol(curr.betas)
    k = nrow(curr.betas)
    
    if (p1 > 0 ) {
      if (ncol(curr.xmat)!=k || nrow(curr.xmat)!=n) {
        stop(paste0(err.txt,"Dimensions of xmat, areas and betas do not match."))
      }
    }
    if (!is.null(curr.priors)) {
      #p2 = ncol(curr.priors)
      p2 = length(nonmixed_priors)
      p2_mixed = length(mixed_priors)
      if (any(curr.priors<0)) {stop(paste0(err.txt,"Priors must be strictly non-negative."))}
    } else {p2 = 0;p2_mixed = 0}
    
    # check restrictions for consistency
    if (!is.null(curr.restrictions) & any(colnames(curr.restrictions) %in% names(curr.targets))) {
      restr.mat = matrix(0,n,p); colnames(restr.mat) = names(curr.targets)
      restr.mat[,colnames(curr.restrictions)] = curr.restrictions
    } else {restr.mat = NULL}
    
    # out.res contains downscaled estimates; priors.mu econometric & other priors for estimation
    out.res = priors.mu = matrix(0,n,p)
    colnames(out.res) = colnames(priors.mu) = names(curr.targets)
    # match econometric priors
    priors.mu[,colnames(curr.betas)] = curr.xmat %*% curr.betas
    # make sure the priors are numerically well behaved
    priors.mu[priors.mu > options$MAX_EXP] = options$MAX_EXP
    priors.mu[priors.mu < -options$MAX_EXP] = -options$MAX_EXP
    priors.mu[,colnames(curr.betas)] = exp(priors.mu[,colnames(curr.betas)])
    # match other priors (if they exist)
    if (p2 > 0) {
      priors.mu[,nonmixed_priors] = curr.priors[,nonmixed_priors]
    }
    if (p2_mixed > 0) {
      w1 = curr.prior_weights[,mixed_priors,drop = FALSE]
      #re-scale exogeneous prior to priors.mu
      eco.priors_sum = apply(priors.mu[,mixed_priors,drop = FALSE],c(2),sum)
      exo.priors = curr.priors[,mixed_priors,drop = FALSE]
      exo.priors_sum = apply(exo.priors, 2, sum)
      exo.priors = t(
        (t(exo.priors) / exo.priors_sum) * eco.priors_sum   )
      priors.mu[,mixed_priors] = as.matrix((1-w1)*priors.mu[,mixed_priors] + w1*exo.priors)
    }
    # remove targets that are all zero
    not.zero = (curr.targets != 0)
    if (all(curr.targets == 0)) {
      
      #catch case if all targets are equal zero
      out.solver[[curr.lu.from]] = NULL
    } else {
      
      #cut out zero targets from targets and priors
      if (any(curr.targets == 0) && !all(curr.targets == 0)) {
        curr.targets = curr.targets[not.zero]
        priors.mu = priors.mu[,not.zero,drop = FALSE]
        if (!is.null(curr.restrictions)) {restr.mat = restr.mat[,not.zero,drop = FALSE]}
      }
      
      #proceed with bias correction
      x0 = curr.targets / sum(curr.targets + 1)
      opts <- list(algorithm = options$algorithm,
                   xtol_rel = options$xtol_rel,
                   xtol_abs = options$xtol_abs,
                   maxeval = options$maxeval
      )
      # check if the optimiser uses gradient or not
      if (grepl("_LD_",opts$algorithm)) {
        eval_grad_f = grad_sqr_diff.mnl
      } else {
        eval_grad_f = NULL
      }
      res.x = nloptr::nloptr(x0 = x0,
                             eval_f = sqr_diff.mnl,
                             eval_grad_f = eval_grad_f,
                             lb = rep(exp(-options$MAX_EXP),length(x0)),
                             ub = rep(exp(options$MAX_EXP),length(x0)),
                             opts=opts,
                             mu = priors.mu,areas = curr.areas,targets = curr.targets,
                             restrictions = restr.mat,cutoff = options$cutoff)
      res.x$par = res.x$solution
      
      out.mu = mu.mnl(res.x$solution[1:length(curr.targets)],priors.mu,curr.areas,restr.mat,options$cutoff)
      
      # Check where out.mu deviates more from the target than max_diff
      error_ind = (curr.targets - colSums(out.mu))^2 > options$max_diff
      # If any diff is larger, do individual logit models boosted by grid search
      if (any(error_ind)) {
        error_targets = curr.targets[error_ind]
        error_restrictions = restr.mat[,error_ind,drop=FALSE]
        
        # Calculate residual areas -  res_areas
        res_areas = curr.areas
        if (any(!error_ind)) {
          res_areas = res_areas - rowSums(out.mu[,!error_ind,drop=FALSE])
        }
        
        # Loop over remaining targets
        for (ccc in 1:length(error_targets)) {
          curr_error_target = error_targets[ccc]
          curr_error_restrictions = error_restrictions[,ccc]
          curr_error_mu = priors.mu[,names(curr_error_target),drop=FALSE]
          
          # Do an iterated grid search to find correct scaling coefficient for priors
          curr_scaling = iterated_grid_search(min_param = -options$MAX_EXP,
                                              max_param = options$MAX_EXP,
                                              func = sqr_diff.mnl,
                                              max_iterations = 10,
                                              precision_threshold = 1e-3,
                                              exp_transform = TRUE,
                                              mu = curr_error_mu,areas = res_areas,
                                              targets = curr_error_target,
                                              restrictions = curr_error_restrictions,
                                              cutoff = options$cutoff)
          
          # Re-scale prior
          curr_error_mu = curr_error_mu * curr_scaling$best_param
          
          # Optimize with scaled priors
          res.x = nloptr::nloptr(x0 = 1,
                                 eval_f = sqr_diff.mnl,
                                 eval_grad_f = eval_grad_f,
                                 lb = exp(-options$MAX_EXP),
                                 ub = exp(options$MAX_EXP),
                                 opts=opts,
                                 mu = curr_error_mu,areas = res_areas,targets = curr_error_target,
                                 restrictions = curr_error_restrictions,cutoff = options$cutoff)
          
          # Calculate areas of current target with mu.mnl
          curr_error_out.mu =
            mu.mnl(res.x$solution,
                   curr_error_mu,res_areas,curr_error_restrictions,
                   options$cutoff)
          
          # Add calculated areas to out.mu
          out.mu[,names(curr_error_target)] = curr_error_out.mu
          
          # Substract areas from res.areas
          res_areas = res_areas - curr_error_out.mu
          
          # Add note to res.x
          res.x$message = "INDIVIDUAL LOGIT BOOSTED BY GRID SEARCH: Standard optimization failed to converge"
        }
      }
      
      if (all(not.zero)) {out.res = out.mu
      } else {out.res[,not.zero] = out.mu}
      out.solver[[curr.lu.from]] = res.x
    }
    
    # add residual own flows in output
    out.res2 = data.frame(ns = names(curr.areas),
                          curr.areas - rowSums(out.res),out.res)
    colnames(out.res2)[2] = paste0(curr.lu.from)
    # pivot into long format
    res.agg <- out.res2 %>%
      pivot_longer(cols = -c("ns"),names_to = "lu.to") %>%
      bind_cols(lu.from = curr.lu.from)
    
    # aggregate results over dataframes
    if(curr.lu.from==lu.from[1]){
      full.out.res <- res.agg
    } else {
      full.out.res = bind_rows(full.out.res,res.agg)
    }
  }
  return(list(out.res = full.out.res, out.solver = out.solver))
}

environment(solve_biascorr.mnl) <- asNamespace("downscalr")
assignInNamespace("solve_biascorr.mnl", solve_biascorr.mnl, ns = "downscalr")
message("Patched downscalr::solve_biascorr.mnl (error_restrictions drop=FALSE fix)")
# ==================================================================
# USER SETTINGS
# ==================================================================

# Select the configuration file to use for this run.
# The file must be located in the config/ folder.
config_file <- "BRA.yml"


# ------------------------------------------------------------------
# 1. Load packages
# ------------------------------------------------------------------

# Restore the project environment once after cloning the repository:
#install.packages("renv")
#renv::restore()

# In normal use, load the installed package.
# In developer mode, devtools::load_all() will already have loaded it.
if (!"package:FABLEDownscalR" %in% search()) {
  library(FABLEDownscalR)
}

library(dplyr)
library(here)
library(ggnewscale)
library(ggpattern)
library(rnaturalearth)
library(rnaturalearthdata)
library(countrycode)
library(tidyr)
library(stringr)

source(here::here("R", "run_downscale_restricted.R"))  # adjust path as needed
# 1. Restore the environment (once)
renv::restore()   # then restart R

Carbon_LUC <- readRDS(here::here("R", "grid50_carbon_stocks.rds")
                      
# ------------------------------------------------------------------
# 2. Read configuration
# ------------------------------------------------------------------

cfg <- fdr_read_config(
  here::here("config", config_file)
)

# Automatically create the date stamp used in output filenames
cfg$stamp <- format(Sys.Date(), "%y%m%d")

# Ensure reproducible MNL estimation
set.seed(cfg$seed)


# ------------------------------------------------------------------
# 3. Prepare country information
# ------------------------------------------------------------------

country <- countrycode::countrycode(
  cfg$country,
  origin = "iso3c",
  destination = "country.name"
)

border_sf <- rnaturalearth::ne_countries(
  country = country,
  scale = "medium",
  returnclass = "sf"
)

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
  map_LUC            = inputs$mapping$map_LUC,
  Ts                 = 2015
)

transitions_table <- inputs$mapping$map_LUC

transitions_table %>%
  dplyr::filter(from == "otherland", to == "otherland")

# (5) a  Harmonise start map to match FABLE baseline totals
harm_starting <- fdr_harmonize_start_map(
  LandCoverStarting_df  = inputs$spatial$landcoverstarting,
  LandCoverInitial_df   = inputs$spatial$landcoverinitial,
  type          = "starting",
  map_LC        = inputs$mapping$map_LC,
  LC_targets    = inputs$LC_targets
)

# (5) b Harmonise initial map to match FABLE baseline totals
harm_initial <- fdr_harmonize_start_map(
  LandCoverStarting_df  = inputs$spatial$landcoverstarting,
  LandCoverInitial_df   = inputs$spatial$landcoverinitial,
  type = "initial",
  map_LC        = inputs$mapping$map_LC,
  LC_targets    = inputs$LC_targets
)

# (6) Build ns_map + rasterized ID layer (resolution controlled by YAML)
id <- fdr_build_id_maps(
  grid_sp         = inputs$grid_sp,      # e.g. Travel as sp/sf geometry for cells
  ns_map          = harm_starting$ns_map,
  pixel_res_m     = cfg$pixel_res_m
)

# (7) Build priors (X matrix etc.) + drop cells with NA covariates
priors <- fdr_build_priors(
  inputs         = inputs,
  start_map      = harm_initial$start_map_reproj,
  good_ns_only   = TRUE
)

# (8) Build restrictions (from Areas_Protegidas_FABLE.xls, sheet "FABLE")

# Three protected-area categories (all in ha):
#   Proteção Integral  — strict protection (no extractive use)
#   Terra Indígena     — indigenous territory
#   Uso Sustentável    — sustainable use (limited extractive allowed)
# Any cell with total PA area > 0 is restricted: Cropland, Pasture, and Urban
# expansion into it is forbidden (max.change = 0 / value = 1).

message("Reading Areas_Protegidas_FABLE.xls | sheet: FABLE ...")
pa_raw <- readxl::read_xls(here::here("R", "Areas_Protegidas_FABLE.xls"), sheet = "FABLE")

# Select id_c + all (ha) columns
pa_ha_cols <- names(pa_raw)[grepl("\\(ha\\)", names(pa_raw))]
message("PA (ha) columns: ", paste(pa_ha_cols, collapse = ", "))

pa_df <- pa_raw %>%
  select(id_c...1, all_of(pa_ha_cols)) %>%
  mutate(
    id_c     = as.character(id_c...1),
    pa_total = rowSums(across(all_of(pa_ha_cols)), na.rm = TRUE)
  )

protected_cells <- pa_df %>% filter(pa_total > 0) %>% pull(id_c)
message("Protected cells (PA overlap): ", length(protected_cells), " / ", nrow(pa_df))

all_restricted_cells <- protected_cells
message("Total restricted cells (PA only): ", length(all_restricted_cells))

restrictions <- expand.grid(lu.from = "forest", lu.to = c("pasture", "cropland", "otherland", "urban"),
                               stringsAsFactors = FALSE) %>%
  filter(lu.from != lu.to) %>%
  tidyr::expand_grid(ns = all_restricted_cells) %>%
  mutate(value = 1) %>%
  select(ns, lu.from, lu.to, value)

restrictions <- as.data.frame(restrictions)



# (9) Fit MNL + downscale
results <- fdr_run_downscaling_restricted(
  targets      = fdr_wrangle_fable_targets(inputs$FABLE_targets, min_year = 2020),
  country_luc  = luc$country_luc,
  priors       = priors,
  mnl_niter    = cfg$mnl_niter,
  mnl_nburn    = cfg$mnl_nburn,
  EF_LUC       = inputs$EF_LUC,
  Carbon_LUC   = Carbon_LUC,
  restrictions = restrictions
)

# (10) Save outputs (consistent naming)
fdr_save_outputs(
  country     = cfg$country,
  tag         = fdr_make_tag(cfg),
  output_root = cfg$output_root,
  outputs     = list(
    start_map_reproj   = harm_starting$start_map_reproj,
    ns_map             = harm_starting$ns_map,
    rasterized_layer   = id$rasterized_layer,
    grid_sf            = id$grid_sf,
    X_long             = results$X_long,
    betas              = results$pred_coeff_long,
    country_start_areas= results$country_start_areas,
    downscaled_LUC     = results$downscaled_LUC,
    luc_hist           = luc$country_luc
  )
)
