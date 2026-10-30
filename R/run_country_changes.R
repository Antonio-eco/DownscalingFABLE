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

brazil_FABLE <- readxl::read_excel("Data/BRA/FABLE.xlsx", sheet = "CT") %>%
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

LU_CLASSES <- c("forest", "pasture", "cropland", "otherland", "urban")

brazil_luc <- luc$country_luc %>% dplyr::mutate(ns = as.character(ns))

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
  log_MeanTravelTime     = 10,
  alllivestock  = 10,
  log_pop        = 10,
  yield_level_mean        = 10
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

# ── Run DownscalR spatial allocation with econometric betas + flat priors ───
country_start_areas <- priors$lu_levels %>%
  dplyr::mutate(times = min(brazil_FABLE$times))

downscale_out <- downscalr::downscale(
  targets      = brazil_FABLE,
  start.areas  = country_start_areas,
  xmat         = X_long,
  betas        = betas_all,
  priors       = flat_priors,
  restrictions = restrictions
)

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

downscale_out$out.res <- downscale_out$out.res %>%
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
  out.res             = downscale_out$out.res,
  downscaled_LUC      = downscale_out$out.res,
  X_long              = X_long,
  pred_coeff_long     = betas_all,
  country_start_areas = country_start_areas
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

downscaled_LUC <- results$out.res
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



