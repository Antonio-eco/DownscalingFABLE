# 🌍 FABLE Spatial Downscaling — Brazil Run Script (`run_country_changes.R`)

This document explains, in plain language, what `run_country_changes.R` actually does, and how it differs from the standard workflow described in the project's main `README.md`.

> **Bottom line:** this is not a generic "run any country" script anymore. It is a **Brazil-specific working script** that bolts a bug fix, several manual overrides, and two new analysis steps (carbon emissions + results verification) onto the standard FABLE downscaling workflow. Several of its inputs point to files on one person's computer rather than the shared `Data/<COUNTRY>/` folder, so as written it will not run as-is on another machine.

------------------------------------------------------------------------

## 1. What the script does, step by step

| \# | Step | Plain-language description |
|------------|------------|------------------------------------------------|
| 0 | **Bug fix (patch)** | Before anything else runs, the script repairs a known bug in the underlying `downscalr` package (an error that occurs in rare edge cases during bias correction) by replacing the faulty function in memory. No files are permanently changed — the fix only applies while this script is running. |
| 1 | **Load packages** | Loads the FABLE downscaling package plus several mapping/formatting packages (for country borders, pretty charts, and country-code lookups). |
| 2 | **Read configuration** | Reads settings from a config file — but the file is hardcoded to `BRA.yml` (Brazil) rather than being easily swapped. |
| 3 | **Prepare country information** | Looks up Brazil's official name and border shape, used later for mapping/context. |
| 4–6 | **Load & harmonize data** | Loads Brazil's spatial and FABLE Calculator data, builds the historical land-cover-change table, and adjusts ("harmonizes") the starting land-cover maps so they match the official FABLE national totals. |
| 7 | **Build the prediction grid** | Converts the country's grid cells into a standard raster layer at the resolution set in the config file. |
| 8 | **Build model inputs ("priors")** | Builds the table of landscape characteristics (travel time, population, livestock, crop yield, etc.) used to predict where land-use change is likely to happen, keeping only valid grid cells. |
| 9 | **Set land-use restrictions** | Loads a separate Protected Areas spreadsheet and flags any grid cell overlapping a protected area so that Forest cannot convert to Cropland or Pasture there. |
| 10 | **Manual target harmonization** | Rescales certain land-use totals (forest, pasture, cropland, other land, urban) to match manually specified target figures. |
| 11 | **Fit the statistical model** | For each land-use class, fits a model predicting where that land is likely to change, using Brazil's historical data and the landscape characteristics from step 8. |
| 12 | **Blend in a "flat" prior** | Mixes in a simple, uniform assumption (land change is equally likely everywhere) with the statistical model's predictions, to prevent a handful of unusual grid cells from dominating the results. |
| 13 | **Run the spatial allocation ("downscaling")** | Combines the national FABLE targets, the statistical model, the flat-prior blend, and the protected-area restrictions to spatially distribute Brazil's projected land-use change across the grid, year by year. |
| 14 | **Estimate greenhouse-gas emissions from land-use change** | Using separate carbon-stock data, calculates the CO₂-equivalent emissions (or removals) associated with each land conversion, for example when forest is cleared or new forest grows back. |
| 15 | **Save outputs** | Writes all key results (maps, model coefficients, downscaled projections, emissions) to the standard `Output/<COUNTRY>/` folder. |
| 16 | **Verify the results** | Automatically checks that the areas reported at the end of one time period match the areas reported at the start of the next, and that gains/losses add up correctly — flags any mismatches. |
| 17 | **Summarize total area over time** | Produces a simple year-by-year table of total area by land-use class, so results can be sanity-checked at a glance. |

------------------------------------------------------------------------

## 2. Key differences from the main project README

The original `README.md` describes a **generic, reusable workflow**: any country team drops its files into `Data/<COUNTRY_CODE>/`, fills in a config file, and runs `source("R/run_country.R")`. This script deviates from that in several important ways:

### 🔧 Added a package bug fix

- The script includes a **runtime patch** for a bug in the `downscalr` package's bias-correction function (it can crash in a specific edge case). This fix is not mentioned anywhere in the original README and is not part of a package release — it's applied fresh each time the script runs.

### 🇧🇷 Hardcoded to one country, not configurable "out of the box"

- `config_file <- "BRA.yml"` is hardcoded, rather than something easily swapped per the README's "replace IND.yml with your configuration file" instruction.
- A Brazil-specific FABLE sheet name (`"CT"`) and class-name translations (Forest → forest, etc.) are written directly into the script.

### 📁 Uses files outside the standard `Data/<COUNTRY>/` folder

- Two inputs are read from **absolute paths on one person's computer** (`C:/Users/User/Desktop/...`), not from the shared `Data/BRA/` folder the README specifies:
  - A carbon-stock dataset (`grid50_carbon_stocks.rds`)
  - A Protected Areas spreadsheet (`Areas_Protegidas_FABLE.xls`)
- This means the script **will not run for another team member** without adjusting those paths — it is not yet "packaged" the way the README describes for general use.

### 🎯 Manual, hardcoded overrides not documented in the README

- A set of **manual target totals** for forest, pasture, cropland, other land, and urban areas is hardcoded into the script (`HARMONIZE_TARGETS`) and used to rescale the model's starting values. This kind of manual adjustment isn't part of the documented workflow.
- Custom **model-sensitivity settings** (`A0_PER_COVARIATE`) are hardcoded for specific landscape variables (travel time, livestock, population, crop yield) to prevent a few extreme data points from distorting results — a refinement added to address a specific accuracy issue observed during testing.
- A **"flat prior" blend** (`PRIOR_WEIGHT = 0.7`) is introduced, mixing the statistical model's predictions with a simple uniform assumption. This is a new modeling choice not mentioned in the README.

### 🌱 New step: greenhouse-gas emissions calculation

- The script adds an entirely new analysis step that was **not described in the original README at all**: estimating CO₂-equivalent emissions or removals from each land-use change, using an external carbon-stock dataset.

### ✅ New step: built-in results verification

- The script adds **automatic consistency checks** after the run (do areas balance between periods? do gains and losses add up?) and prints a clear OK/FAILED summary. The original README has no equivalent quality-check step — previously, this would need to be checked manually.

### 📊 New step: summary table of total area by year

- A simple, easy-to-read table showing total area per land-use class for every year of the projection is generated at the end — a convenience addition not present in the original workflow.

### 📦 Extra package requirements not listed in the README

- The script relies on several packages the README's "Requirements" section does not mention, including `here`, `ggnewscale`, `ggpattern`, `rnaturalearth`, `rnaturalearthdata`, `countrycode`, `stringr`, `readxl`, and `nloptr`.

------------------------------------------------------------------------

## 3. What this means for other country teams

Because this script hardcodes Brazil-specific values, file paths, and target numbers, it currently functions more as a **working/debugging version of the Brazil run** than a drop-in replacement for the generic `R/run_country.R` described in the main README. Before another team can reuse it, it would need:

1.  The two hardcoded absolute file paths replaced with files placed in `Data/<COUNTRY>/`, following the standard structure.
2.  The hardcoded config file name, FABLE sheet name, and manual target totals turned into configurable parameters (e.g., read from the `.yml` config).
3.  The new GHG-emissions and verification steps documented as optional/standard parts of the pipeline, with their own input requirements added to the "Data Preparation" section of the main README.

------------------------------------------------------------------------

## 📩 Support

For questions about this script or the underlying workflow, please contact: António Santos: [acasantos\@ciencias.ulisboa.pt](mailto:acasantos@ciencias.ulisboa.pt){.email} Or open an issue on GitHub.
