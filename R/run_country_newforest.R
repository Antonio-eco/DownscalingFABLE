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

Carbon_LUC <- readRDS("C:/Users/User/Desktop/FABLE_Spatial_LUC_EF-main/output/run_20260928_1644/grid50_carbon_stocks.rds")
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



# (8)  Build restrictions  (from Areas_Protegidas_FABLE.xls, sheet "FABLE")

# Three protected-area categories (all in ha):
#   Proteção Integral  — strict protection (no extractive use)
#   Terra Indígena     — indigenous territory
#   Uso Sustentável    — sustainable use (limited extractive allowed)
# Any cell with total PA area > 0 is restricted: Cropland, Pasture, and Urban
# expansion into it is forbidden (max.change = 0 / value = 1).

message("Reading Areas_Protegidas_FABLE.xls | sheet: FABLE ...")
pa_raw <- readxl::read_xls("C:/Users/User/Desktop/FABLE/downscalr/Areas_Protegidas_FABLE.xls", sheet = "FABLE")

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

all_lu_classes <- c("pasture", "cropland", "otherland", "urban")  # match your actual LU class names

restrictions <- expand.grid(
  lu.from = "forest",
  lu.to   = c("pasture", "cropland"),
  stringsAsFactors = FALSE
) %>%
  filter(lu.from != lu.to) %>%
  tidyr::expand_grid(ns = all_restricted_cells) %>%
  mutate(value = 1) %>%
  select(ns, lu.from, lu.to, value)

restrictions <- as.data.frame(restrictions)

# NOTE: NewForest is NOT given an explicit forest/urban -> newforest
# restriction here. It isn't needed -- brazil_FABLE's newforest_targets below
# only ever has lu.from %in% c("otherland","cropland","pasture") (per the NDC
# sheet), so lu.from == "forest"/"urban" never has a "newforest" column in
# curr.targets to begin with inside solve_biascorr.mnl, making that
# transition structurally impossible already. Adding a restriction row for it
# anyway actively BREAKS solve_biascorr.mnl: its restr.mat is built with
# colnames = names(curr.targets) (i.e. only classes that ARE targets for that
# lu.from), but the restriction-matching check only verifies
# any(colnames(curr.restrictions) %in% names(curr.targets)) instead of all(),
# so it still tries restr.mat[, "newforest"] <- ... for forest/urban, which
# don't have that column, producing "subscript out of bounds".

# ── Harmonize toward manually-specified, per-class target totals ───────────
# Each class in HARMONIZE_TARGETS is rescaled independently: every cell's
# value for that class is multiplied by (target / current class total), so
# the class hits the target total exactly while every cell keeps the same
# SHARE of that class's total it had before (no reshuffling across cells,
# no change to other classes). Units: raw ha, matching priors$lu_levels.
# Leave a class out of the vector to leave it untouched.
HARMONIZE_TARGETS <- c(
  "forest"   = 577126,
  "pasture"  = 159279.2196,
  "cropland" = 76951,
  "otherland" = 33408,
  "urban" = 3916
)

if (length(HARMONIZE_TARGETS) > 0) {
  current_totals <- priors$lu_levels %>%
    filter(lu.from %in% names(HARMONIZE_TARGETS)) %>%
    group_by(lu.from) %>%
    summarise(current = sum(value), .groups = "drop")
  
  scale_tbl <- current_totals %>%
    mutate(
      target = HARMONIZE_TARGETS[lu.from],
      factor = target / current
    )
  
  message("\nHarmonizing priors$lu_levels to manual per-class targets:")
  print(as.data.frame(scale_tbl %>% select(lu.from, current, target, factor)),
        row.names = FALSE)
  
  priors$lu_levels <- priors$lu_levels %>%
    left_join(scale_tbl %>% select(lu.from, factor), by = "lu.from") %>%
    mutate(value = ifelse(!is.na(factor), value * factor, value)) %>%
    select(-factor)
  
  message("priors$lu_levels: total AFTER harmonization = ",
          round(sum(priors$lu_levels$value) / 1e3, 2), " M ha")
}

brazil_FABLE <- readxl::read_excel("Data/BRA/FABLE.xlsx", sheet = "NDC") %>%
  rename(lu.from = LandCoverInit, times = YearEnd) %>%
  select(-YearStart) %>%
  pivot_longer(starts_with("To"), names_to = "lu.to", values_to = "value") %>%
  mutate(
    lu.to = sub("^To", "", lu.to),
  ) %>%
  select(lu.from, times, lu.to, value) %>%
  mutate(lu.from = recode(lu.from,
                     "Forest" = "forest",
                     "Cropland" = "cropland",
                     "Pasture" = "pasture",
                     "OtherLand" = "otherland",
                     "Urban" = "urban"),
    lu.to = recode(lu.to,
                   "Forest" = "forest",
                   "Cropland" = "cropland",
                   "Pasture" = "pasture",
                   "OtherLand" = "otherland",
                   "Urban" = "urban")) %>%
  filter(lu.from != "NewForest", lu.to != "NewForest", times >= 2020)

message("brazil_FABLE: ", nrow(brazil_FABLE), " rows | periods: ",
        paste(sort(unique(brazil_FABLE$times)), collapse = ", "))

# ── NewForest targets, read directly from the "NDC" sheet ───────────────────
# The "CT" sheet (above) has no meaningful NewForest schedule, so NewForest is
# dropped from it entirely. The "NDC" sheet's ToNewForest column already IS
# the complete, correct schedule -- no extrapolation needed: 6 Mha flow in
# 2020-2025 (OtherLand+Cropland+Pasture = 638.58+1297.15+4064.27 = 6,000
# [1000 ha]), another 6 Mha flow in 2025-2030 (same split, same total -- this
# is cumulative 12 Mha reached BY 2030, not a typo/stale value), then flat at
# ZERO further flow for every period after 2030 through 2050 (confirmed
# directly from the sheet: all six LandCoverInit classes show ToNewForest = 0
# for 2030-2035 onward). So the cumulative NewForest stock is 0 (2020) -> 6
# Mha (2025) -> 12 Mha (2030) -> stays at 12 Mha through 2050.
#
# Source classes are {OtherLand, Cropland, Pasture} -> NewForest (land that is
# ITSELF afforested), per the sheet's literal LandCoverInit values -- this
# still differs from the verbal "forest->pasture / forest->cropland /
# pasture->otherland redirected" framing from before, but is corroborated by
# the GHG post-processing code elsewhere in this script, which independently
# hardcodes ef_biomass for exactly lu.to=="newforest" & lu.from %in%
# c("cropland","pasture","otherland"). Flag back if the redirect framing was
# actually intended instead.

newforest_targets <- readxl::read_excel("Data/BRA/FABLE.xlsx", sheet = "NDC") %>%
  dplyr::filter(LandCoverInit %in% c("OtherLand", "Cropland", "Pasture"),
                YearEnd >= 2020) %>%
  dplyr::transmute(
    lu.from = recode(LandCoverInit,
                      "OtherLand" = "otherland",
                      "Cropland"  = "cropland",
                      "Pasture"   = "pasture"),
    times   = YearEnd,
    lu.to   = "newforest",
    value   = ToNewForest
  )

newforest_cumulative_by_period <- newforest_targets %>%
  dplyr::group_by(times) %>%
  dplyr::summarise(flow_Mha = sum(value) / 1000, .groups = "drop") %>%
  dplyr::arrange(times) %>%
  dplyr::mutate(cumulative_Mha = cumsum(flow_Mha))

message("NewForest schedule, read from NDC sheet:")
print(as.data.frame(newforest_cumulative_by_period), row.names = FALSE)

stopifnot(
  "NewForest cumulative total by 2030 should be 12 Mha" =
    abs(newforest_cumulative_by_period$cumulative_Mha[
      newforest_cumulative_by_period$times == 2030] - 12) < 0.01,
  "NewForest flow after 2030 should be zero (stays flat at 12 Mha)" =
    all(newforest_cumulative_by_period$flow_Mha[
      newforest_cumulative_by_period$times > 2030] == 0)
)

brazil_FABLE <- dplyr::bind_rows(brazil_FABLE, newforest_targets) %>%
  dplyr::arrange(lu.from, lu.to, times)

# (9) Fit MNL (per-covariate priors) + flat-prior blend + downscale ----------
# Replaces the call to fdr_run_downscaling_restricted()'s internal MNL step
# (which used a single scalar A0 ridge-grid and which.max(colSums(Y)) as the
# baseline class) with: (a) an explicit per-covariate prior-variance vector
# A0 passed straight into downscalr::mnlogit(), using each class's OWN column
# as baseline, and (b) a uniform "flat prior" blended in at the solve_biascorr
# stage via the `priors`/`weight` columns (see the patched solve_biascorr.mnl
# at the top of this script). downscalr::downscale() is then called directly
# so both of these feed the spatial allocation step, instead of only the
# econometric betas.

LU_CLASSES <- c("forest", "pasture", "cropland", "otherland", "urban", "newforest")

brazil_luc <- luc$country_luc %>% dplyr::mutate(ns = as.character(ns))

# ── NewForest covariates: high_vigor / medium_vigor / low_vigor ─────────────
# Built from the "vigor_pastagem" theme in the geoparquet (classe: alto /
# intermediario / degradado = high / medium / low pasture vigor). Uses the
# most recent available year (2020) as the static covariate, same convention
# as every other time-invariant covariate in xmat. `valor` is stored as text
# with a comma decimal separator (Brazilian format), so it's parsed before
# use. Values are converted to within-cell SHARES (alto/(alto+interm+degr))
# rather than raw areas, so they're on a comparable [0,1] scale to other
# covariates and aren't dropped as near-constant/mis-scaled during cleaning.
vigor_raw <- arrow::read_parquet(
  "R/Long_Fable_ago.geoparquet",
  col_select = c("id_c", "ano", "tema", "classe", "valor")
)

vigor_wide <- vigor_raw %>%
  dplyr::filter(tema == "vigor_pastagem", ano == "2020") %>%
  dplyr::mutate(
    ns    = as.character(id_c),
    valor = as.numeric(gsub(",", ".", valor))
  ) %>%
  dplyr::select(ns, classe, valor) %>%
  tidyr::pivot_wider(names_from = classe, values_from = valor, values_fill = 0) %>%
  dplyr::mutate(
    vigor_total  = alto + intermediario + degradado,
    high_vigor   = ifelse(vigor_total > 0, alto / vigor_total, 0),
    medium_vigor = ifelse(vigor_total > 0, intermediario / vigor_total, 0),
    low_vigor    = ifelse(vigor_total > 0, degradado / vigor_total, 0)
  ) %>%
  dplyr::select(ns, high_vigor, medium_vigor, low_vigor)

message("vigor_pastagem (2020): ", nrow(vigor_wide), " cells with vigor data ",
        "(cells outside this set get high/medium/low_vigor = 0)")

vigor_long <- vigor_wide %>%
  tidyr::pivot_longer(-ns, names_to = "ks", values_to = "value")

# Append to priors$xmat BEFORE the cleaning/standardisation step below, so
# high_vigor/medium_vigor/low_vigor go through the same sd-filter/QR
# full-rank treatment as every other covariate. Cells present in priors$xmat
# but absent from vigor_wide implicitly get 0 for all three columns once
# pivot_wider(values_fill = 0) runs just below.
priors$xmat <- dplyr::bind_rows(priors$xmat, vigor_long)

# ── Build + clean the covariate matrix (same cleaning fdr_run_downscaling_
#    restricted() did internally: drop near-constant cols, standardise,
#    drop collinear cols via QR) ──────────────────────────────────────────
fdr_make_full_rank_X <- function(X, tol_sd = 1e-10, tol_qr = 1e-6, scale_X = TRUE) {
  stopifnot(is.matrix(X))
  sds <- apply(X, 2, sd, na.rm = TRUE)
  keep_sd <- which(is.finite(sds) & sds > tol_sd)
  X <- X[, keep_sd, drop = FALSE]
  if (ncol(X) < 2) stop("X has <2 usable columns after sd filtering.")
  if (scale_X) {
    mu  <- colMeans(X, na.rm = TRUE)
    sdv <- apply(X, 2, sd, na.rm = TRUE)
    sdv[!is.finite(sdv) | sdv < 1e-12] <- 1
    X <- sweep(sweep(X, 2, mu, "-"), 2, sdv, "/")
  }
  q <- qr(X, tol = tol_qr, LAPACK = TRUE)
  keep_qr <- q$pivot[seq_len(q$rank)]
  X <- X[, keep_qr, drop = FALSE]
  list(X = X, keep_cols = colnames(X))
}

X_long_raw <- priors$xmat

xmat_full <- X_long_raw %>%
  tidyr::pivot_wider(names_from = "ks", values_from = "value", values_fill = 0) %>%
  dplyr::arrange(ns) %>%
  tibble::column_to_rownames("ns") %>%
  as.matrix()

clean    <- fdr_make_full_rank_X(xmat_full, tol_sd = 1e-10, tol_qr = 1e-6, scale_X = TRUE)
xmat_std <- clean$X                 # standardised, full-rank covariate matrix (rows = ns)
keep_ks  <- clean$keep_cols

# Keep long xmat consistent with the cleaned/standardised matrix -- this is
# what gets passed into downscalr::downscale() as `xmat`, so its `value`
# column must match xmat_std, not the raw (unstandardised) priors$xmat.
X_long <- tibble::as_tibble(xmat_std, rownames = "ns") %>%
  tidyr::pivot_longer(-ns, names_to = "ks", values_to = "value")

# ── Per-covariate prior variance (A0) ────────────────────────────────────────
# mnlogit() builds beta_prior_var <- diag(k) * A0 (see downscalr::mnlogit).
# A0 as a scalar gives every covariate the same diffuse prior (1e4, near-flat).
# Passing a length-k VECTOR gives each covariate its own prior variance on the
# diagonal. Near-separation in a handful of extreme-covariate training cells
# combined with a flat A0=1e4 lets posterior betas explode (observed:
# travel_time->Urban beta=9.19, log_livestock->Urban beta=-7.79,
# travel_time->Pasture beta=8.27 -- all from cells with extreme standardised
# travel_time, ~3+ SD out). These extreme betas make the same handful of
# cells look overwhelmingly attractive across MULTIPLE destination classes,
# and because xmat is static across all 7 chained periods in downscale(),
# the same distortion re-applies every period, cascading land through those
# cells (Forest->Cropland one period, that Cropland->Urban the next, etc.)
# rather than spreading it plausibly. Tightening these covariates' priors
# shrinks the posterior toward zero unless the data strongly supports it.

A0_DEFAULT <- 1e4
A0_PER_COVARIATE <- c(
  log_MeanTravelTime = 100,   # was "travel_time" -- never matched, silently fell
                              # back to A0_DEFAULT=1e4, causing the near-
                              # separation Inf/NaN betas below
  log_pop            = 100,
  alllivestock        = 100,  # was "log_livestock" -- same issue; no single
  buffalo             = 100,  # aggregated livestock column exists, so the
  cattle              = 100,  # tightened prior is applied to each species
  chicken             = 100,
  goats               = 100,
  pigs                = 100,
  sheep               = 100,
  # Forest-management type shares: sparse, near-zero for most cells and
  # extreme for a handful -- the same near-separation shape that caused the
  # travel_time/livestock blow-up, confirmed here with Agroforestry producing
  # NA/Inf betas for every lu.from under the untightened A0_DEFAULT=1e4.
  Agroforestry        = 100,
  ForestPlanted       = 100,
  ForestRegNoHuman    = 100,
  ForestRegWithHuman  = 100,
  OilPalmPlantation   = 100,
  SRP                 = 100
  # No "log_crop" equivalent in this dataset's xmat_full -- see covariate_
  # weights note above re: yield_level_mean.
)
build_A0_vector <- function(covariate_names) {
  a0 <- rep(A0_DEFAULT, length(covariate_names))
  names(a0) <- covariate_names
  override <- intersect(names(A0_PER_COVARIATE), covariate_names)
  a0[override] <- A0_PER_COVARIATE[override]
  a0
}

active_cells <- brazil_luc %>%
  filter(lu.from != lu.to, value > 0) %>%
  distinct(ns) %>% pull(ns)

message("active_cells: ", length(active_cells))

betas_list <- list()

for (cl in LU_CLASSES) {

  Yraw   <- brazil_luc %>% filter(lu.from == cl & Ts == 2015)
  Y_full <- Yraw %>%
    pivot_wider(id_cols = ns, names_from = lu.to, values_from = value,
                values_fill = 0) %>%
    tibble::column_to_rownames("ns")
  for (mc in setdiff(LU_CLASSES, colnames(Y_full))) Y_full[[mc]] <- 0
  Y_full <- Y_full[, LU_CLASSES]

  baseline <- which(colnames(Y_full) == cl)
  keep     <- rownames(Y_full) %in% active_cells & rowSums(Y_full) > 0
  Y        <- Y_full[keep, , drop = FALSE]

  if (nrow(Y) == 0) {
    betas_list[[cl]] <- expand.grid(ks = colnames(xmat_full),
                                     lu.to = setdiff(LU_CLASSES, cl),
                                     stringsAsFactors = FALSE) %>%
      mutate(lu.from = cl, value = 0) %>% select(ks, lu.to, value, lu.from)
    message(sprintf("%-10s: 0 active cells -> placeholder betas", cl)); next
  }

  Y <- Y / rowSums(Y)
  X <- xmat_std[rownames(Y), , drop = FALSE]
  stopifnot(all(is.finite(X)), all(is.finite(as.matrix(Y))))
  message(sprintf("%-10s: %d training cells", cl, nrow(Y)))

  set.seed(42)
  A0_vec <- build_A0_vector(colnames(X))
  res_mnl <- downscalr::mnlogit(as.matrix(X), as.matrix(Y), baseline = baseline,
                      niter = 100, nburn = 50, A0 = A0_vec)

  pred_coeff <- apply(res_mnl$postb[, -baseline, , drop = FALSE], c(1, 2), mean)
  betas_list[[cl]] <- as.data.frame(pred_coeff) %>%
    tibble::rownames_to_column("ks") %>%
    pivot_longer(-ks, names_to = "lu.to", values_to = "value") %>%
    mutate(lu.from = cl)
}

betas_all <- bind_rows(betas_list) %>% select(ks, lu.to, value, lu.from)
stopifnot(all(is.finite(betas_all$value)))
message("betas_all: ", nrow(betas_all), " rows")

# ── Post-estimation covariate reweighting ────────────────────────────────────
# Scales each covariate's posterior-mean beta by a fixed factor, independent
# of the A0 prior-variance shrinkage applied during MNL estimation above.
# A0 controls how much the DATA can move the beta away from 0 during
# estimation; this reweighting instead rescales the final fitted beta itself,
# e.g. to damp or boost specific covariates' influence on the downscale()
# allocation step without re-running mnlogit(). Any `ks` in betas_all not
# named in covariate_weights gets NA from the lookup -- checked below.
covariate_weights <- c(
  # LU baseline shares (otherland is the dropped reference class in
  # fdr_build_priors()'s Xmat_LU, so it has no beta column to reweight)
  forest             = 2.0,
  cropland           = 2.0,
  pasture            = 2.0,
  urban              = 2.0,
  newforest          = 2.0,   # 5th LU baseline column present in this dataset
  # Travel time
  log_MeanTravelTime = 0.3,
  # Population
  log_pop            = 0.3,
  # Livestock -- no single aggregated "log_livestock" column exists here;
  # fdr_build_priors() keeps each species (+ "alllivestock") as its own raw
  # (non-logged) column, so the 0.3 weight is applied to each of them
  alllivestock       = 0.3,
  buffalo            = 0.3,
  cattle             = 0.3,
  chicken            = 0.3,
  goats              = 0.3,
  pigs               = 0.3,
  sheep              = 0.3,
  # Terrain
  log_altitude       = 1.0,
  log_slope          = 1.0
  # No "log_crop" equivalent in this dataset's xmat_full (closest is
  # yield_level_mean, left at default 1.0 below -- confirm if that's
  # what was intended and I'll add it explicitly).
)

unmatched_ks <- setdiff(unique(betas_all$ks), names(covariate_weights))
if (length(unmatched_ks) > 0) {
  warning("covariate_weights has no entry for: ",
          paste(unmatched_ks, collapse = ", "),
          " -- defaulting their weight to 1.0 (no reweighting). ",
          "Add them to covariate_weights to control this explicitly.")
}

w_lookup <- setNames(rep(1.0, length(unmatched_ks)), unmatched_ks)
covariate_weights_full <- c(covariate_weights, w_lookup[setdiff(names(w_lookup), names(covariate_weights))])

betas_all <- betas_all %>%
  mutate(value = value * covariate_weights_full[ks])

stopifnot(all(is.finite(betas_all$value)))
message("betas_all (reweighted): ", nrow(betas_all), " rows")

# ── Diagnostics: find what's feeding a NaN into downscale()'s objective ─────
# "objective in x0 returns NA" from nloptr means the objective is already
# NaN at the initial point -- almost always exp(beta * x) overflowing to Inf
# somewhere (Inf/Inf = NaN in the multinomial-logit share calc), or NA/Inf in
# xmat or start.areas feeding straight into it. betas_all being all-finite
# does NOT rule out an overflow: a finite but large beta times a standardised
# covariate can still send exp() to Inf.
message("max |beta| in betas_all: ", max(abs(betas_all$value)))
print(betas_all %>% dplyr::arrange(dplyr::desc(abs(value))) %>% head(15))

country_start_areas <- priors$lu_levels %>%
  dplyr::mutate(times = min(brazil_FABLE$times)) %>%
  dplyr::bind_rows(
    tibble::tibble(ns = unique(as.character(priors$lu_levels$ns)),
                   lu.from = "newforest", value = 0,
                   times = min(brazil_FABLE$times))
  )

stopifnot("X_long has non-finite values" = all(is.finite(X_long$value)))
stopifnot("country_start_areas has non-finite values" =
            all(is.finite(country_start_areas$value)))

# Safety net: clip extreme betas regardless of which covariate they come from,
# rather than chasing each sparse/dummy covariate's A0 one at a time. Pick
# BETA_CLIP based on what the max |beta| print above actually shows -- 10 is
# a reasonable starting point (exp(10) ~ 22000, exp(20) ~ 4.9e8, both still
# representable; it's betas in the 30s-100s, as seen with the untightened
# A0_DEFAULT=1e4 covariates, that overflow to Inf in the solver).
BETA_CLIP <- 10
n_clipped <- sum(abs(betas_all$value) > BETA_CLIP)
if (n_clipped > 0) {
  message(n_clipped, " beta(s) exceed |", BETA_CLIP, "| and will be clipped -- ",
          "consider tightening A0_PER_COVARIATE for the covariate(s) shown above instead.")
}
betas_all <- betas_all %>%
  mutate(value = pmin(pmax(value, -BETA_CLIP), BETA_CLIP))

# ── OPTION 2: flat/uniform priors, blended with the econometric (beta) mu ────
# PRIOR_WEIGHT is the fraction of each cell's predicted transition propensity
# that comes from a UNIFORM prior instead of the MNL/beta model, blended as:
#   priors.mu = (1 - weight) * econometric_mu + weight * flat_prior
# (flat_prior is rescaled internally by solve_biascorr.mnl to preserve each
# class's total propensity mass, so this REDISTRIBUTES the same implied pull
# more evenly across cells rather than concentrating it in a handful of
# extreme-covariate "winner" cells). Applied as a constant weight across all
# rows -- same intent as a global "prior_weights" scalar, implemented via the
# per-row `weight` column that complete_priors()/solve_biascorr.mnl support.
PRIOR_WEIGHT <- 0.7

flat_priors <- betas_all %>%
  distinct(lu.from, lu.to) %>%
  tidyr::crossing(ns = unique(X_long$ns)) %>%
  mutate(value = 1, weight = PRIOR_WEIGHT) %>%
  select(ns, lu.from, lu.to, value, weight)

message("flat_priors: ", nrow(flat_priors), " rows | weight = ", PRIOR_WEIGHT)

# ── OPTION 3: existing-area prior -- bias allocation toward cells that ──────
# ALREADY hold more of the destination class, instead of spreading it evenly
# (flat_priors) or leaving it purely to the econometric betas (priors=NULL).
# This is the literal opposite of flat_priors: value = each cell's CURRENT
# area of lu.to (from priors$lu_levels at the baseline year), not a uniform 1.
# Blended the same way inside solve_biascorr.mnl:
#   priors.mu = (1 - weight) * econometric_mu + weight * rescaled_existing_area
# so a bigger EXISTING_AREA_WEIGHT concentrates more of each class's target
# growth onto cells that already have more of that class (agglomeration /
# persistence effect), rather than letting new area appear anywhere the
# covariates happen to favor.
EXISTING_AREA_WEIGHT <- 0.7

baseline_areas_by_class <- priors$lu_levels %>%
  dplyr::mutate(ns = as.character(ns)) %>%
  dplyr::group_by(ns, lu.from) %>%
  dplyr::summarise(area = sum(value, na.rm = TRUE), .groups = "drop")

existing_area_priors <- betas_all %>%
  distinct(lu.from, lu.to) %>%
  tidyr::crossing(ns = unique(X_long$ns)) %>%
  dplyr::left_join(
    baseline_areas_by_class %>% dplyr::rename(lu.to = lu.from, value = area),
    by = c("ns", "lu.to")
  ) %>%
  dplyr::mutate(
    value  = tidyr::replace_na(value, 0) + 1e-6,  # solve_biascorr.mnl requires priors >= 0;
                                                   # tiny epsilon avoids an all-zero column
                                                   # for a (ns, class) with literally no
                                                   # existing area, which would otherwise
                                                   # zero out that class everywhere
    weight = EXISTING_AREA_WEIGHT
  ) %>%
  select(ns, lu.from, lu.to, value, weight)

message("existing_area_priors: ", nrow(existing_area_priors),
        " rows | weight = ", EXISTING_AREA_WEIGHT)

# ── NewForest allocation prior, driven by pasture vigor ─────────────────────
# betas_all has essentially NO estimable signal for "-> newforest": no
# historical cell ever transitioned to newforest (it didn't exist as a class
# before this policy target), so Y is all-zero for that outcome in every
# per-class MNL fit above, and the fitted beta is a degenerate near-separation
# value (clipped to -BETA_CLIP) carrying no real information about WHERE
# afforestation should happen. To make high_vigor/medium_vigor/low_vigor
# actually drive the spatial pattern of new forest (rather than the clipped,
# uninformative econometric betas), this builds a dedicated prior for
# (pasture|cropland|otherland) -> newforest rows, scored by vigor composition,
# with a HIGH blend weight -- since there's no real econometric signal to
# weigh it against for this specific transition.
VIGOR_SCORE_WEIGHTS   <- c(high_vigor = 3, medium_vigor = 1, low_vigor = 0)
NEWFOREST_PRIOR_WEIGHT <- 0.95  # close to 1: the prior does nearly all the
                                 # within-class allocation work for newforest

newforest_prior <- tidyr::expand_grid(
  lu.from = c("pasture", "cropland", "otherland"),
  lu.to   = "newforest"
) %>%
  tidyr::crossing(ns = unique(X_long$ns)) %>%
  dplyr::left_join(
    vigor_wide %>%
      dplyr::mutate(vigor_score = VIGOR_SCORE_WEIGHTS["high_vigor"]   * high_vigor +
                                   VIGOR_SCORE_WEIGHTS["medium_vigor"] * medium_vigor +
                                   VIGOR_SCORE_WEIGHTS["low_vigor"]    * low_vigor) %>%
      dplyr::select(ns, vigor_score),
    by = "ns"
  ) %>%
  dplyr::mutate(
    value  = tidyr::replace_na(vigor_score, 0) + 1e-6,  # see existing_area_priors note on
                                                          # the epsilon -- avoids an all-zero
                                                          # column for cells with no vigor data
    weight = NEWFOREST_PRIOR_WEIGHT
  ) %>%
  dplyr::select(ns, lu.from, lu.to, value, weight)

message("newforest_prior: ", nrow(newforest_prior),
        " rows | weight = ", NEWFOREST_PRIOR_WEIGHT,
        " | vigor weights: high=", VIGOR_SCORE_WEIGHTS["high_vigor"],
        " medium=", VIGOR_SCORE_WEIGHTS["medium_vigor"],
        " low=", VIGOR_SCORE_WEIGHTS["low_vigor"])

# Splice into flat_priors / existing_area_priors: drop their generic
# (*, newforest) rows (value=1 or value=baseline-area=0, neither reflects
# vigor) and replace with the vigor-scored prior above. The "none" mode
# (priors = NULL) deliberately stays untouched by this -- it means no prior
# at all, including for newforest, which is a legitimate choice but WILL
# mean newforest allocation ignores vigor entirely in that mode.
flat_priors <- dplyr::bind_rows(
  flat_priors %>% dplyr::filter(lu.to != "newforest"),
  newforest_prior
)
existing_area_priors <- dplyr::bind_rows(
  existing_area_priors %>% dplyr::filter(lu.to != "newforest"),
  newforest_prior
)

# ── Toggle: which prior feeds the PRIMARY run (the one that continues on
#    into GHG post-processing, `results`, and fdr_save_outputs below) ────────
# "flat"     -> econometric betas + uniform blend (spreads allocation evenly),
#               newforest rows use the vigor prior instead of the uniform one
# "existing" -> econometric betas + existing-area blend (concentrates
#               allocation onto cells that already hold more of that class),
#               newforest rows use the vigor prior instead (baseline newforest
#               area is 0 everywhere, so the generic existing-area version
#               would be uninformative for it anyway)
# "none"     -> econometric betas only, priors = NULL (no blending at all,
#               INCLUDING for newforest -- vigor has no effect in this mode)
PRIOR_MODE <- "existing"

country_start_areas <- priors$lu_levels %>%
  dplyr::mutate(times = min(brazil_FABLE$times)) %>%
  dplyr::bind_rows(
    tibble::tibble(ns = unique(as.character(priors$lu_levels$ns)),
                    lu.from = "newforest", value = 0,
                    times = min(brazil_FABLE$times))
  )

downscale_out <- downscalr::downscale(
  targets      = brazil_FABLE,
  start.areas  = country_start_areas,
  xmat         = X_long,
  betas        = betas_all,
  priors       = flat_priors,
  restrictions = restrictions
)

# ── Comparison run: SAME betas, but WITHOUT the flat-prior blend ────────────
# This isolates the effect of PRIOR_WEIGHT on the spatial allocation itself:
# same econometric betas_all, same X_long/restrictions/targets, just
# priors = NULL instead of flat_priors, so solve_biascorr.mnl falls back to
# its normal (econometric-mu-only) behaviour with no uniform-prior blending.
downscale_out_normal <- downscalr::downscale(
  targets      = brazil_FABLE,
  start.areas  = country_start_areas,
  xmat         = X_long,
  betas        = betas_all,
  priors       = NULL,
  restrictions = restrictions
)

# ── Existing-area-prior run: SAME betas, blended with existing_area_priors ──
downscale_out_existing <- downscalr::downscale(
  targets      = brazil_FABLE,
  start.areas  = country_start_areas,
  xmat         = X_long,
  betas        = betas_all,
  priors       = existing_area_priors,
  restrictions = restrictions
)

# Pick which of the three feeds everything downstream (GHG calc, `results`,
# fdr_save_outputs). The other two stay available for comparison maps.
downscale_primary <- switch(PRIOR_MODE,
  flat     = downscale_out,
  existing = downscale_out_existing,
  none     = downscale_out_normal,
  stop("PRIOR_MODE must be 'flat', 'existing', or 'none'")
)
message("Primary downscale run uses PRIOR_MODE = '", PRIOR_MODE, "'")

# ── Carried over from fdr_run_downscaling_restricted(): GHG due to LUC, by
#    ecoregion, computed on the downscaled transitions above ─────────────────
EF_Pools_transition_Ecoregion <- inputs$EF_LUC %>%
  mutate(to = (ifelse(to == "forest", "newforest", to))) %>%
  filter(!is.na(ECO_NAME))

grid50_carbon_stocks <- Carbon_LUC %>%
  select(
    id_c,
    biomass_total_pasture,
    biomass_total_forest,
    biomass_total_otherland,
    biomass_total_urban,
    biomass_total_cropland
  )

downscale_primary$out.res <- downscale_primary$out.res %>%
  left_join(
    EF_Pools_transition_Ecoregion %>%
      select(id_c, from, to, ef_biomass),
    by = c("ns" = "id_c", "lu.from" = "from", "lu.to" = "to")
  ) %>%
  left_join(
    grid50_carbon_stocks,
    by = c("ns" = "id_c")
  ) %>%
  mutate(
    ef_biomass = case_when(
      lu.to == "newforest" &
        lu.from %in% c("cropland", "pasture", "otherland") ~
        -(biomass_total_forest / 50),
      lu.to == "otherland" &
        lu.from %in% c("cropland", "pasture", "forest") ~
        -(biomass_total_otherland / 20),
      TRUE ~ ef_biomass
    ),
    GHG_biomass = ef_biomass * value * 3.667 / 1000
  ) %>%
  arrange(ns, lu.from, lu.to, times) %>%
  group_by(ns, lu.from, lu.to) %>%
  mutate(
    GHG_biomass = if_else(
      (lu.to == "newforest" &
         lu.from %in% c("cropland", "pasture", "otherland")) |
        (lu.to == "otherland" &
           lu.from %in% c("cropland", "pasture", "forest")),
      cumsum(GHG_biomass),
      GHG_biomass
    )
  ) %>%
  ungroup() %>%
  select(
    -biomass_total_pasture,
    -biomass_total_forest,
    -biomass_total_otherland,
    -biomass_total_urban,
    -biomass_total_cropland
  )

results <- list(
  out.res             = downscale_out $out.res,
  downscaled_LUC      = downscale_out $out.res,
  X_long              = X_long,
  pred_coeff_long     = betas_all,
  country_start_areas = country_start_areas
)

results_normal <- list(
  out.res             = downscale_out_normal$out.res,
  downscaled_LUC      = downscale_out_normal$out.res,
  X_long              = X_long,
  pred_coeff_long     = betas_all,
  country_start_areas = country_start_areas
)


results_existing <- list(
  out.res             = downscale_out_existing $out.res,
  downscaled_LUC      = downscale_out_existing $out.res,
  X_long              = X_long,
  pred_coeff_long     = betas_all,
  country_start_areas = country_start_areas
)
# =============================================================================
# Comparison: flat-prior-blended downscale vs normal (econometric-only) ------
# =============================================================================
# Per-cell, per-year, per-class NET land use (sum of lu.to == class across all
# lu.from, minus lu.from == class across all lu.to, cumulated from start areas)
# is heavier than needed here -- what matters for a map is simply "how much of
# class k does cell ns hold in year t", i.e. the allocated AREA, which the
# downscale() out.res already gives as flows; so we reconstruct per-class
# standing area by class at the final period directly from out.res + starting
# areas. The quickest robust way: take the LAST period's lu.to inflow total per
# (ns, class) plus un-converted "stayed" area, which downscalr itself returns
# as lu.from == lu.to rows (self-transitions) in out.res, so simply filtering
# lu.to == class and summing value, per ns & time, gives the END-of-period
# area of that class in that cell (sum over all lu.from, incl. self).

compare_class_area <- function(out.res, scenario_label, cls) {
  out.res %>%
    dplyr::filter(lu.to == cls) %>%
    dplyr::group_by(ns, times) %>%
    dplyr::summarise(area = sum(value, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(scenario = scenario_label)
}

MAP_CLASS <- "forest"                       # change to compare a different class
MAP_YEAR  <- max(brazil_FABLE$times)        # final period; change as needed

area_flat   <- compare_class_area(downscale_out$out.res,        "flat_prior",   MAP_CLASS)
area_normal <- compare_class_area(downscale_out_normal$out.res, "normal_prior", MAP_CLASS)

area_compare <- dplyr::bind_rows(area_flat, area_normal) %>%
  tidyr::pivot_wider(names_from = scenario, values_from = area, values_fill = 0) %>%
  dplyr::mutate(diff = flat_prior - normal_prior)

message(sprintf(
  "%s area at t=%s | flat_prior total: %.1f | normal_prior total: %.1f | mean |diff| per cell: %.3f",
  MAP_CLASS, MAP_YEAR,
  sum(area_compare$flat_prior), sum(area_compare$normal_prior),
  mean(abs(area_compare$diff))
))

# ── Rasterize + map, using terra + discrete breaks + class palettes ─────────
# Same convention as 04_generate_maps.R: terra::classify() onto an id raster
# keyed by the numeric cell id, a fixed discrete break scale (not a continuous
# gradient) so area classes are visually comparable across scenarios/runs, and
# a per-class RColorBrewer palette for the two allocation maps, with a
# diverging RdBu palette (symmetric around 0) for the difference map.
library(terra)
library(RColorBrewer)

scale_breaks <- c(0, 0.001, 5.9, 16.8, 37.3, 71.1, 130.4, 211.5, 260, 310)

class_palettes <- list(
  forest    = brewer.pal(9, "Greens"),
  cropland  = brewer.pal(9, "Reds"),
  pasture   = brewer.pal(9, "Purples"),
  otherland = brewer.pal(9, "RdPu"),
  urban     = brewer.pal(9, "Greys"),
  newforest = brewer.pal(9, "Greens")
)

diff_ceil <- 310
diff_pal  <- colorRampPalette(brewer.pal(11, "RdBu"))(1001)

# id_raster: a SpatRaster whose cell values are the numeric cell id (ns/id_c),
# built in step (6) (fdr_build_id_maps). Mirrors 04_generate_maps.R's
# rast("data/luc/id_raster.tif").
id_raster <- id$rasterized_layer

to_raster <- function(df_ns_value) {
  if (nrow(df_ns_value) == 0) {
    return(classify(id_raster, cbind(-Inf, Inf, 0)))
  }
  # Same "_a" suffix convention as 04_generate_maps.R's load_luc(): some ns
  # values are like "12345_a" (not plain numeric), which as.numeric() alone
  # turns into NA (the warning seen above). Those get id_c = <numeric part> +
  # 1e6, matching the id_raster's own convention for those cells.
  # ifelse() evaluates BOTH branches for every row, so as.numeric(ns) still
  # runs (and warns) on the "_a" rows even though that result is discarded,
  # and as.numeric(sub("_a$","",ns)) still runs on the non-"_a" rows likewise
  # -- cosmetic noise, not a sign of bad data. suppressWarnings() here is safe
  # because genuinely unparseable ns values are still caught by the
  # is.finite(df_id$id) check and reported via message() right below.
  df_id <- df_ns_value %>%
    dplyr::mutate(
      id = suppressWarnings(
        ifelse(grepl("_a$", ns),
               as.numeric(sub("_a$", "", ns)) + 1e6,
               as.numeric(ns))
      )
    )
  n_bad <- sum(!is.finite(df_id$id))
  if (n_bad > 0) {
    message(n_bad, " ns value(s) didn't parse to a numeric id and were dropped: ",
            paste(head(df_id$ns[!is.finite(df_id$id)], 5), collapse = ", "),
            if (n_bad > 5) ", ..." else "")
  }
  reclass_mat <- as.matrix(
    df_id %>% dplyr::filter(is.finite(id)) %>% dplyr::transmute(id, value)
  )
  classify(id_raster, reclass_mat, others = NA)
}

area_compare_yr <- area_compare %>% dplyr::filter(times == MAP_YEAR)

r_flat   <- to_raster(area_compare_yr %>% dplyr::select(ns, value = flat_prior))
r_normal <- to_raster(area_compare_yr %>% dplyr::select(ns, value = normal_prior))
r_diff   <- r_flat - r_normal

pal <- class_palettes[[MAP_CLASS]]
if (is.null(pal)) {
  warning("No class_palettes entry for '", MAP_CLASS, "' -- falling back to Greens.")
  pal <- brewer.pal(9, "Greens")
}

# id$grid_sf as a boundary overlay -- same role as brazil_states/brazil_biomes
# in 04_generate_maps.R, and doubles as a visual check here that the raster
# (id$rasterized_layer) and the polygons actually line up: if the class map
# comes out blank but the grid outlines still show, the raster's id values
# don't match what to_raster()'s reclass_mat is keying on. Converted to
# SpatVector once, reused by both plot functions below.
grid_vect <- tryCatch(terra::vect(id$grid_sf), error = function(e) NULL)
if (is.null(grid_vect)) {
  warning("Could not convert id$grid_sf to a SpatVector for overlay -- ",
          "plots will show the raster without cell boundaries.")
}

plot_class_map <- function(r, title_text) {
  plot(r,
       col    = pal,
       breaks = scale_breaks,
       legend = "bottomright",
       mar    = c(1.5, 0.5, 4, 0.5),
       plg    = list(title = "1000 ha", cex = 0.8))
  if (!is.null(grid_vect)) terra::plot(grid_vect, border = "lightgray", lwd = 0.1, add = TRUE)
  title(main = title_text, adj = 1, cex.main = 0.9)
}

plot_diff_map <- function(r, title_text) {
  mm <- as.numeric(terra::minmax(r, compute = TRUE))
  if (is.finite(mm[1]) && mm[1] == mm[2]) {
    message("No difference to plot for ", title_text)
    return(invisible(NULL))
  }
  plot(r,
       col    = diff_pal,
       range  = c(-diff_ceil, diff_ceil),
       legend = "bottomright",
       mar    = c(1.5, 0.5, 4, 0.5),
       plg    = list(title = "1000 ha", cex = 0.8))
  if (!is.null(grid_vect)) terra::plot(grid_vect, border = "lightgray", lwd = 0.1, add = TRUE)
  title(main = title_text, adj = 1, cex.main = 0.9)
}

plot_class_map(r_flat,
               sprintf("%s area, %s -- flat prior (weight=%s)", MAP_CLASS, MAP_YEAR, PRIOR_WEIGHT))
plot_class_map(r_normal,
               sprintf("%s area, %s -- normal prior (econometric only)", MAP_CLASS, MAP_YEAR))
plot_diff_map(r_diff,
              sprintf("%s area diff, %s: flat_prior - normal_prior", MAP_CLASS, MAP_YEAR))

# To save as PNGs instead of/in addition to on-screen plotting, wrap each
# plot_*_map() call in png(...)/dev.off(), same as 04_generate_maps.R:
#   png("flat.png", width = 820, height = 780, res = 100)
#   plot_class_map(r_flat, "...")
#   dev.off()

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
    X_long             = results_existing$X_long,
    betas              = results_existing$pred_coeff_long,
    country_start_areas= results_existing$country_start_areas,
    downscaled_LUC     = results_existing$downscaled_LUC,
    luc_hist           = luc$country_luc
  )
)

downscaled_LUC <- results_existing$out.res
# =============================================================================
# M.  Verification: AreaEnd(period t) == AreaStart(period t+1), per class
#     (same continuity logic as the FABLE verification workbook:
#      LandCoverInit | YearStart | YearEnd | AreaStart | AreaEnd | Check)
# =============================================================================
# In downscaled_LUC, each row is a (ns, lu.from, lu.to, times, value) transition
# for the period ENDING in `times`. For a given class c and period t:
#   AreaStart(t, c) = sum(value | lu.from == c, times == t)  -- everything c
#                     had at the start of the period, however it ends up
#   AreaEnd(t, c)   = sum(value | lu.to   == c, times == t)  -- everything
#                     that is class c at the end of the period, from any source
# Continuity requires AreaEnd(t, c) == AreaStart(t+1, c) for consecutive
# periods, i.e. the closing stock downscale() produced for a period must be
# exactly the opening stock it fed into the next period.

area_start <- downscaled_LUC %>%
  mutate(times = as.integer(as.character(times))) %>%
  group_by(lu.from, times) %>%
  summarise(AreaStart = sum(value, na.rm = TRUE), .groups = "drop") %>%
  rename(LandCoverInit = lu.from, YearEnd = times)

area_end <- downscaled_LUC %>%
  mutate(times = as.integer(as.character(times))) %>%
  group_by(lu.to, times) %>%
  summarise(AreaEnd = sum(value, na.rm = TRUE), .groups = "drop") %>%
  rename(LandCoverInit = lu.to, YearEnd = times)

# TotalGains(t, c)  = area that flowed INTO c from every OTHER class this period
# TotalLosses(t, c) = area that flowed OUT of c into every OTHER class this period
# (mirrors the FABLE workbook's TotalGains/TotalLosses columns; the diagonal
#  lu.from == lu.to entries are "stayed as c" and are excluded from both, since
#  they're neither a gain nor a loss for c)
total_gains <- downscaled_LUC %>%
  mutate(times = as.integer(as.character(times))) %>%
  filter(lu.from != lu.to) %>%
  group_by(lu.to, times) %>%
  summarise(TotalGains = sum(value, na.rm = TRUE), .groups = "drop") %>%
  rename(LandCoverInit = lu.to, YearEnd = times)

total_losses <- downscaled_LUC %>%
  mutate(times = as.integer(as.character(times))) %>%
  filter(lu.from != lu.to) %>%
  group_by(lu.from, times) %>%
  summarise(TotalLosses = sum(value, na.rm = TRUE), .groups = "drop") %>%
  rename(LandCoverInit = lu.from, YearEnd = times)

verify_tbl <- full_join(area_start, area_end, by = c("LandCoverInit", "YearEnd")) %>%
  full_join(total_gains,  by = c("LandCoverInit", "YearEnd")) %>%
  full_join(total_losses, by = c("LandCoverInit", "YearEnd")) %>%
  mutate(
    TotalGains  = tidyr::replace_na(TotalGains, 0),
    TotalLosses = tidyr::replace_na(TotalLosses, 0)
  ) %>%
  # downscaled_LUC / out.res is in raw ha throughout the script; convert to
  # 1000 ha (kha) here so this table matches the FABLE workbook's units.
  mutate(
    AreaStart   = AreaStart   / 1000,
    AreaEnd     = AreaEnd     / 1000,
    TotalGains  = TotalGains  / 1000,
    TotalLosses = TotalLosses / 1000
  ) %>%
  arrange(LandCoverInit, YearEnd) %>%
  group_by(LandCoverInit) %>%
  mutate(
    YearStart      = YearEnd - 5L,
    NextAreaStart  = lead(AreaStart),                       # AreaStart (kha) of the following period
    CheckContinuity = NextAreaStart - AreaEnd,              # should be 0: this period's close == next period's open
    CheckBalance    = AreaEnd - (AreaStart + TotalGains - TotalLosses)  # should be 0: gains/losses reconcile start->end
  ) %>%
  ungroup() %>%
  select(LandCoverInit, YearStart, YearEnd, AreaStart, AreaEnd,
         TotalGains, TotalLosses, NextAreaStart, CheckContinuity, CheckBalance)

# tolerance is in kha now (1e-6 kha = 1e-3 ha) -- still far tighter than the
# floating-point residuals (~1e-10 ha) seen in the FABLE workbook itself
TOL <- 1e-6
continuity_breaks <- verify_tbl %>% filter(abs(CheckContinuity) > TOL & !is.na(CheckContinuity))
balance_breaks    <- verify_tbl %>% filter(abs(CheckBalance)    > TOL & !is.na(CheckBalance))

cat("\n=== Area verification (per class, per period) ===\n")
print(as.data.frame(verify_tbl), row.names = FALSE)

if (nrow(continuity_breaks) == 0) {
  cat("\nOK: all", nrow(verify_tbl), "class-period joins are continuous (|diff| <=", TOL, ").\n")
} else {
  cat("\nFAILED:", nrow(continuity_breaks), "class-period joins are NOT continuous (|diff| >", TOL, "):\n")
  print(as.data.frame(continuity_breaks), row.names = FALSE)
}

if (nrow(balance_breaks) == 0) {
  cat("OK: AreaStart + TotalGains - TotalLosses == AreaEnd for all", nrow(verify_tbl), "rows (|diff| <=", TOL, ").\n")
} else {
  cat("FAILED:", nrow(balance_breaks), "rows where AreaStart + TotalGains - TotalLosses != AreaEnd:\n")
  print(as.data.frame(balance_breaks), row.names = FALSE)
}

# =============================================================================
# N.  Total area by land use, per time step (wide format, 1000 ha)
# =============================================================================
# One row per time step (2000, 2005, ..., 2050), one column per land-use class,
# plus a Total column (sum across classes -- should be constant/near-constant
# across time steps, since downscale() conserves total area).
#
# For 2000 (base year) we use AreaStart of the first period (2000-2005), since
# downscaled_LUC itself only carries transitions for periods ending 2005+.
# For 2005-2050 we use AreaEnd of the period ending in that year, which is the
# total stock of each class at that point in time.

area_2000 <- verify_tbl %>%
  filter(YearStart == min(YearStart)) %>%
  select(LandCoverInit, Year = YearStart, Area = AreaStart)

area_rest <- verify_tbl %>%
  select(LandCoverInit, Year = YearEnd, Area = AreaEnd)

area_by_timestep <- bind_rows(area_2000, area_rest) %>%
  distinct(LandCoverInit, Year, .keep_all = TRUE) %>%
  arrange(Year, LandCoverInit) %>%
  tidyr::pivot_wider(names_from = LandCoverInit, values_from = Area) %>%
  mutate(Total = rowSums(across(-Year), na.rm = TRUE)) %>%
  arrange(Year)

cat("\n=== Total area by land use, per time step (1000 ha) ===\n")
print(as.data.frame(area_by_timestep), row.names = FALSE)



