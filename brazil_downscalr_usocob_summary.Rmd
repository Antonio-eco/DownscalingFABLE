---
title: "Brazil DownscaleR — UsoCobertura (MapBiomas) Downscaling Summary"
subtitle: "Summary of `brazil_downscalr_usocob.R`"
date: "`r Sys.Date()`"
output:
  html_document:
    toc: true
    toc_depth: 3
    toc_float: true
    number_sections: true
    theme: flatly
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(eval = FALSE, echo = TRUE)
```

# Purpose

This document summarizes the workflow implemented in `brazil_downscalr_usocob.R`,
which downscales Brazil's FABLE land-use-change (LUC) targets from national/state
level down to the grid-cell level, using the **`downscalr`** package. The key
change from earlier versions of the pipeline is that the starting land-use stock
now comes from **MapBiomas UsoCobertura** (`UsoCobertura__fable_v2.xlsx`, sheet
`usoecob2020`) instead of `hilda_br_2015.xlsx`. All other inputs (HILDA+
transition flows, FABLE targets, and covariates) are unchanged.

This is **not executable code** — it is a narrative walkthrough of the logic,
inputs, outputs, and known issues in the original script. Code chunks are shown
for reference (`eval = FALSE`) rather than run.

# Inputs

All files are expected in the working directory:

| File | Role |
|---|---|
| `UsoCobertura__fable_v2.xlsx` (sheet `usoecob2020`) | MapBiomas land-use/land-cover stock |
| `hildaluc_br_2015_2019.xlsx` | HILDA+ transition flows, 2015–2019 |
| `brazil_fable_ct.xlsx` | FABLE Calculator targets |
| `altitude_br.xlsx`, `slope_br.xlsx`, `travel_time_br.xlsx` | Terrain/accessibility covariates |
| `livestock_br.xlsx`, `pop_2020.xlsx`, `crop_yield_br.xlsx` | Socioeconomic/production covariates |
| `Areas_Protegidas_FABLE.xls` (sheet `FABLE`) | Protected-area (PA) polygons/areas used to build restrictions |
| `br_states.shp`, `br_biomes.shp` | State boundaries and biome polygons for a geographic exclusion rule |

# Package Setup and a Runtime Patch

```{r libs}
library(dplyr); library(tidyr); library(tibble)
library(readxl); library(downscalr)
```

The script begins by **monkey-patching** `downscalr::solve_biascorr.mnl()` at
runtime via `assignInNamespace()`. This fixes a real bug: when exactly one
`(lu.from, lu.to)` target triggers the bias-correction fallback, indexing
`restr.mat[, error_ind]` silently drops to a 1-D vector, which then breaks a
later `[, ccc]` subset with a dimensions error. The patch adds `drop = FALSE`
so the restriction matrix always stays 2-D. The patched function is installed
into the `downscalr` namespace, so it applies for the rest of the session
without needing to reinstall the package.

# Section-by-Section Walkthrough

## A. Land-Use Class Mapping (`col_to_lu`)

Defines the **five target land-use classes** used throughout the pipeline:

- `Forest`, `OtherLand`, `Cropland`, `Pasture`, `Urban`

A named vector `col_to_lu` maps each raw MapBiomas class column (e.g.
`Savanna`, `Mangrove`, `Soybean`, `Mining`, `PhotovoltaicPP`) onto one of these
five classes. Columns not listed (`NoData`, `Aquaculture`, `ComForestry`) are
deliberately excluded. This is flagged in the header as the single place to
edit if the classification needs to change.

## B. Building the Land-Use Stock (`lu_stock`)

- Reads `usoecob2020` from the UsoCobertura workbook.
- Keeps `id_c` (renamed to `ns`, the spatial cell ID) and all columns whose
  name contains `(ha)`.
- Strips the `(ha)` suffix, then **sums the mapped columns per cell** into
  five aggregated columns (`lu_Forest`, `lu_OtherLand`, `lu_Cropland`,
  `lu_Pasture`, `lu_Urban`) — this is the per-cell starting stock, in hectares.
- Emits a warning if any expected mapped column is missing from the sheet, and
  logs total hectares per class.

## C. FABLE Targets (`brazil_FABLE`)

- Reads `brazil_fable_ct.xlsx`, drops `.geo`/`system.index` housekeeping
  columns, and renames `LandCoverInit → lu.from`, `YearEnd → times`.
- Pivots the wide `To*` columns into long `(lu.to, value)` pairs.
- **Unit conversion**: FABLE Calculator values are in 1000 ha (kha); this
  multiplies by 1000 so everything downstream is in raw ha, matching the
  stock, covariates, and betas.
- Drops `NewForest` transitions and restricts to `times >= 2000`.

## D. Historical Transitions (`brazil_luc`)

- Reads HILDA+ 2015→2019 transition flows (unchanged from earlier pipeline
  versions).
- Two-digit `from`/`to` HILDA codes are mapped to the five FABLE classes via
  `hilda_labels` (e.g. `11→Urban`, `22→Cropland`, `33→Pasture`, `44→Forest`,
  `55`–`99→OtherLand`).
- Pivoted to long format and aggregated to `(ns, lu.from, lu.to, Ts=2015)`.
  This dataset is used later as the **observed training data** for the MNL
  (multinomial logit) transition model.

## E. Covariates (`brazil_xmat`)

Builds the covariate matrix used to predict transition propensities:

- Joins `lu_stock` (as covariates, i.e. current land-use shares) with:
  altitude, slope, travel time, livestock, population, and crop yield —
  each read per-cell and renamed to `ns`.
- **Median-imputes** any missing values per covariate (`impute_median()`),
  logging how many values were imputed and with what median.
- Log-transforms three highly skewed covariates: `log1p(livestock)`,
  `log1p(pop_2020)`, `log1p(crop_total)`, replacing the raw columns.
- Drops zero-variance covariates, then **standardizes globally**
  (mean-centered, SD-scaled) to produce `xmat_std`.
- Pivots to long format (`X_long`: `ns, ks, value`) — the format `downscalr`
  expects.
- Saves diagnostic intermediates (`xmat_wide_diag.rds`, `X_long_diag.rds`,
  `lu_stock_diag.rds`) for later inspection of area-conservation issues.

## F. Restrictions (`restrictions_br`)

Builds cell-level restrictions preventing expansion of **Cropland, Pasture,
or Urban** into protected or ecologically sensitive land:

1. **Protected areas**: reads `Areas_Protegidas_FABLE.xls`, sums the three PA
   category columns (Proteção Integral, Terra Indígena, Uso Sustentável) per
   cell, and flags any cell with `pa_total > 0` as protected.
2. **Geographic exclusion (deep-interior Amazon)**: a documented, deliberate
   addition. The PA layer alone under-restricts certain deep-interior Amazon
   cells (in Acre, Amazonas, Roraima) that have negligible official PA
   overlap but get assigned implausibly large Cropland stocks by 2050 due to
   extreme standardized `travel_time` values pulling the MNL model. A spatial
   join (via `br_states.shp` state polygons and `br_biomes.shp` Amazon biome
   polygon) identifies cells that are **both** in `AC`, `AM`, or `RR` **and**
   inside the Amazon biome, and adds them to the restriction set. The
   rationale explicitly excludes Rondônia, Pará, and Mato Grosso, since those
   contain the genuine, actively monitored "arc of deforestation" frontier
   that should *not* be blocked.
3. The union of PA-protected and geographically-excluded cells forms
   `all_restricted_cells`. For every such cell, transitions from any class
   into `Pasture`, `Urban`, or `Cropland` are set to `value = 1` (blocked),
   for all `lu.from != lu.to` combinations.

## G. MNL Estimation (`betas_all`)

For each of the five LU classes as an origin (`lu.from`):

- Builds an observed transition-share matrix `Y` from `brazil_luc` (2015),
  restricted to "active" cells (cells with any observed transition and
  `value > 0`), normalized to row-sum-1 shares.
- If a class has zero active training cells, placeholder betas of `0` are
  assigned instead of running the model.
- Otherwise fits `downscalr::mnlogit()` (Bayesian MNL via Gibbs/MCMC,
  100 iterations, 50 burn-in) with a class as baseline, using standardized
  covariates `xmat_std` as predictors.
- **Custom per-covariate priors (`A0`)**: rather than the default flat scalar
  prior variance (`A0 = 1e4`), four covariates prone to near-separation
  (`travel_time`, `log_livestock`, `log_pop`, `log_crop`) are given a tighter
  prior variance (`100`) via `build_A0_vector()`. This is a deliberate fix
  for observed extreme posterior betas (e.g. `travel_time → Urban ≈ 9.19`)
  driven by a handful of extreme-covariate cells, which — because the
  covariate matrix is static across all chained periods in `downscale()` —
  would otherwise cause implausible land-use cascades through the same cells
  every period.
- Posterior mean coefficients (`postb`, excluding baseline) are extracted per
  class and combined into a single long-format `betas_all` table
  (`ks, lu.to, value, lu.from`).

## H. Starting Areas and Flat-Prior Blending

- `br_start_areas`: reshapes `lu_stock` into long format as the `downscale()`
  starting stock.
- **Flat/uniform prior blending (`flat_priors`)**: a second deliberate
  modeling choice. A uniform prior (`value = 1`) is defined for every
  observed `(lu.from, lu.to)` combination across all cells, with a constant
  blend weight `PRIOR_WEIGHT = 0.7`. Inside `solve_biascorr.mnl()`, this
  blends 70% flat prior with 30% econometric (MNL) prediction
  (`priors.mu = (1 - weight) * econometric + weight * flat`), which
  redistributes the pull implied by the model more evenly across cells
  rather than concentrating change in a small number of extreme-covariate
  "winner" cells.

## I. Downscaling (`downscale()`)

Runs `downscalr::downscale()` with:

- `targets = brazil_FABLE` (national/state targets by period)
- `start.areas = br_start_areas` (UsoCobertura 2020 stock)
- `xmat = X_long`, `betas = betas_all` (MNL model)
- `priors = flat_priors` (0.7-weighted flat blend)
- `restrictions = restrictions_br` (PA ∪ deep-interior exclusion)

This is where the runtime-patched `solve_biascorr.mnl()` is actually invoked
across all chained periods.

## J. Save

Persists results as RDS/CSV:

- `results_DS_usocob_states.rds` — full downscale() output (incl. solver info)
- `downscaled_LUC_usocob_states.rds` / `.csv` — long-format downscaled
  transitions (`ns, lu.from, lu.to, times, value`)
- `betas_all_usocob_states.rds`, `start_areas_usocob_states.rds`

## K. Plotting

Produces a reference-style choropleth map of a single class/period slice:

- Rebuilds cell geometries from the `.geo` column in `pop_2020.xlsx`,
  rasterizes to a 0.05° grid via `terra`.
- Applies a fixed, proportionally-rescaled set of breaks (`BREAKS_KHA`) and a
  red sequential palette, with optional state/biome shapefile overlays.
- Saves `landuse_year_pathway_states_breaks.png` and
  `landuse_year_pathway_states_kha.rds`.

## L. Diagnostics

Ad hoc diagnostic block investigating specific "offender" cells
(e.g. `68690_a`) known to oscillate or show implausible area growth, printing
per-period totals and the ratio of 2050-to-2020 total area for a hard-coded
list of cell IDs, plus a global count of cells with `ratio > 2`.

## M. Area-Conservation Verification

The script ends with **two nearly identical verification blocks** (first in
raw ha, then repeated in kha — see notes below) that:

- Compute `AreaStart` (sum of `value` where `lu.from == class`) and `AreaEnd`
  (sum of `value` where `lu.to == class`) per class per period.
- Compute `TotalGains`/`TotalLosses` (excluding the `lu.from == lu.to`
  "stayed as" diagonal).
- Check two invariants per class/period, each against a tolerance of `1e-6`:
  - **Continuity**: this period's `AreaEnd` equals the next period's
    `AreaStart`.
  - **Balance**: `AreaStart + TotalGains − TotalLosses == AreaEnd`.
- Prints pass/fail summaries and writes results to
  `area_continuity_check_usocob_states.csv`.

# Outputs Summary

| File | Description |
|---|---|
| `results_DS_usocob_states.rds` | Full `downscale()` result object |
| `downscaled_LUC_usocob_states.rds`/`.csv` | Downscaled LUC transitions, long format |
| `betas_all_usocob_states.rds` | Estimated MNL coefficients |
| `start_areas_usocob_states.rds` | Starting (2000) land-use areas per cell |
| `landuse_year_pathway_states_breaks.png` | Choropleth map (see note on naming above) |
| `landuse_year_pathway_states_kha.rds` | Data behind the map |
| `area_continuity_check_usocob_states.csv` | Area conservation verification table |
