# 🌍 FABLE Spatial Downscaling — Brazil Run Script (`run_country_newforest.R`)

This document explains, in plain language, what `run_country_newforest.R` actually does, and how it differs from the standard workflow described in the project's main `README.md`.

> **Bottom line:** this is not a generic "run any country" script anymore. It is a **Brazil-specific working script** that bolts a bug fix, several manual overrides, a new "afforestation" land-use class, and extra analysis steps (carbon emissions + results verification) onto the standard FABLE downscaling workflow. Several of its inputs point to files on one person's computer rather than the shared `Data/<COUNTRY>/` folder, so as written it will not run as-is on another machine.

------------------------------------------------------------------------

## 1. What the script does, step by step

| \# | Step | Plain-language description |
|------------|------------|------------------------------------------------|
| 0 | **Bug fix (patch)** | Before anything else runs, the script repairs a known bug in the underlying `downscalr` package (an error that occurs in rare edge cases during bias correction) by replacing the faulty function in memory. No files are permanently changed — the fix only applies while this script is running. |
| 1 | **Load packages** | Loads the FABLE downscaling package plus several mapping/formatting/data-reading packages (for country borders, pretty charts, country-code lookups, and reading the new pasture-condition file format described in step 8). |
| 2 | **Read configuration** | Reads settings from a config file — but the file is hardcoded to `BRA.yml` (Brazil) rather than being easily swapped. |
| 3 | **Prepare country information** | Looks up Brazil's official name and border shape, used later for mapping/context. |
| 4–6 | **Load & harmonize data** | Loads Brazil's spatial and FABLE Calculator data, builds the historical land-cover-change table, and adjusts ("harmonizes") the starting land-cover maps so they match the official FABLE national totals. |
| 7 | **Build the prediction grid** | Converts the country's grid cells into a standard raster layer at the resolution set in the config file. |
| 8 | **Build model inputs ("priors")** | Builds the table of landscape characteristics (travel time, population, livestock, crop yield, forest-management type, etc.) used to predict where land-use change is likely to happen, keeping only valid grid cells. This version adds one more landscape characteristic: a measure of **pasture condition** (what share of each cell's pasture is in good, medium, or poor condition), read from a separate national dataset. |
| 9 | **Set land-use restrictions** | Loads a separate Protected Areas spreadsheet and flags any grid cell overlapping a protected area so that **forest cannot convert to cropland or pasture** there. (Forest-to-urban conversion is not restricted in this version.) |
| 10 | **Manual target harmonization** | Rescales certain land-use totals (forest, pasture, cropland, other land, urban) to match manually specified target figures. |
| 11 | **Add "New Forest" as a sixth land-use class** *(new)* | Alongside forest, cropland, pasture, other land, and urban, the script now tracks a sixth category: **new forest** — land being actively afforested (planted or restored), separate from pre-existing forest. Its national target schedule comes from a different spreadsheet tab than everything else, and the script automatically double-checks that the numbers reach 12 million hectares by 2030 and stay flat afterward through 2050, stopping with an error if that stops being true. New forest can only come from other land, cropland, or pasture — never from forest or urban. |
| 12 | **Fit the statistical model** | For each land-use class, fits a model predicting where that land is likely to change, using Brazil's historical data and the landscape characteristics from step 8. This version fixes a naming mismatch that had been silently disabling the intended "how confident are we in this factor" settings for travel time and livestock, adds a safety cap so an extreme, poorly estimated factor can't break the calculation, and — because **new forest has no historical examples to learn from** — builds a dedicated, separate rule for where new forest should go: cells with healthier/less-degraded pasture are favored as afforestation candidates. |
| 13 | **Blend the model's predictions with a spreading assumption** | The model's raw predictions are mixed with an assumption about how evenly change should be spread across the landscape, to avoid a handful of unusual cells dominating the results. This version runs that blend **three different ways** and lets the team pick which one counts as "official": spread evenly everywhere (the original approach), use the statistical model with no spreading assumption at all, or — the new, **currently official** option — favor cells that **already hold more** of a given land type, on the logic that land-use patterns tend to persist and grow from where they already exist rather than appear evenly or at random. New forest always uses its own pasture-condition-based rule from step 12, regardless of which of the three is chosen as official. A comparison map is generated showing how different the "spread evenly" and "model only" results are from one another, for a chosen class and year. |
| 14 | **Run the spatial allocation ("downscaling")** | Combines the national FABLE targets, the statistical model, the chosen spreading assumption, and the protected-area restrictions to spatially distribute Brazil's projected land-use change across the grid, year by year — run once for each of the three options in step 13. |
| 15 | **Estimate greenhouse-gas emissions from land-use change** | Using separate carbon-stock data, calculates the CO₂-equivalent emissions (or removals) associated with each land conversion. This version adds proper handling for **new-forest and other-land conversions specifically**: instead of a generic emissions factor, it uses each cell's own measured carbon stock, spread out over a realistic number of years for forest to grow back (50 years) or other land to settle (20 years), accumulated over time rather than calculated period-by-period in isolation. |
| 16 | **Save outputs** | Writes all key results (maps, model coefficients, downscaled projections, emissions) to the standard `Output/<COUNTRY>/` folder. |
| 17 | **Verify the results** | Automatically checks that the areas reported at the end of one time period match the areas reported at the start of the next, and that gains/losses add up correctly — flags any mismatches. |
| 18 | **Summarize total area over time** | Produces a simple year-by-year table of total area by land-use class, so results can be sanity-checked at a glance. |

------------------------------------------------------------------------

## 2. Key differences from the main project README

The original `README.md` describes a **generic, reusable workflow**: any country team drops its files into `Data/<COUNTRY_CODE>/`, fills in a config file, and runs `source("R/run_country.R")`. This script deviates from that in several important ways:

### 🔧 Added a package bug fix

- The script includes a **runtime patch** for a bug in the `downscalr` package's bias-correction function (it can crash in a specific edge case). This fix is not mentioned anywhere in the original README and is not part of a package release — it's applied fresh each time the script runs.

### 🇧🇷 Hardcoded to one country, not configurable "out of the box"

- `config_file <- "BRA.yml"` is hardcoded, rather than something easily swapped per the README's "replace IND.yml with your configuration file" instruction.
- Brazil-specific FABLE sheet names and class-name translations (Forest → forest, etc.) are written directly into the script. New forest's target schedule is read from a **second, different sheet** than the rest of the land-use targets.

### 📁 Uses files outside the standard `Data/<COUNTRY>/` folder

- Three inputs are read from **absolute paths on one person's computer** (`C:/Users/User/Desktop/...`), not from the shared `Data/BRA/` folder the README specifies:
  - A carbon-stock dataset (`grid50_carbon_stocks.rds`)
  - A Protected Areas spreadsheet (`Areas_Protegidas_FABLE.xls`)
  - A pasture-condition dataset (`Long_Fable_ago.geoparquet`), used to decide where new forest should be placed
- This means the script **will not run for another team member** without adjusting those paths — it is not yet "packaged" the way the README describes for general use.

### 🎯 Manual, hardcoded overrides not documented in the README

- A set of **manual target totals** for forest, pasture, cropland, other land, and urban areas is hardcoded into the script (`HARMONIZE_TARGETS`) and used to rescale the model's starting values. This kind of manual adjustment isn't part of the documented workflow.
- Custom **model-sensitivity settings** are hardcoded for specific landscape variables (travel time, livestock, population, crop yield, forest-management type) to prevent a few extreme data points from distorting results. A naming mismatch that had silently disabled the intended settings for travel time and livestock has been fixed in this version.
- A **second, separate adjustment** is applied after the model is fit: each landscape factor's influence on the final result is additionally scaled up or down by a fixed amount (e.g., travel time and livestock influence dialed down, each land-use type's own baseline share dialed up). This is on top of, not instead of, the model-sensitivity settings above.
- A **safety cap** prevents an extreme, poorly estimated factor from breaking the calculation.

### 🌳 A new "New Forest" land-use class, with its own location-finding rule

- The model now tracks **six** land-use classes instead of five, adding "new forest" (actively afforested land) alongside forest, cropland, pasture, other land, and urban.
- Because no historical data exists for this class, a **dedicated rule** decides where it goes, based on a new pasture-condition dataset: cells with healthier pasture are favored for afforestation over degraded pasture. This rule is given very high weight (95%) in the calculation, since there's no real historical evidence to weigh it against.
- A **built-in sanity check** verifies the new-forest schedule reaches 12 million hectares cumulatively by 2030 and then holds flat through 2050, and stops the script with an error if the underlying spreadsheet ever stops matching that expectation.

### 🔀 Three alternative "spreading" assumptions instead of one, including a new one

- The original flat/uniform spreading assumption (mixing in "change is equally likely everywhere") is still available.
- A **brand-new option** is added: favor cells that **already hold more** of a given land-use class — the same kind of blend as the flat option, but using each cell's current area of that class as the steering signal (with a 70% weight) instead of a flat value, plus a tiny floor so a class with zero existing presence isn't treated as permanently impossible there. **This new option is the one currently used as the official result.**
- A third option — the statistical model alone, with no spreading assumption at all — is also run for comparison.
- New forest is excluded from all three of these and always uses its own dedicated pasture-condition rule instead (see above).
- A **comparison map** is generated automatically, showing exactly how much the "spread evenly" and "model only" results differ for a chosen class and year, so the effect of this choice can be visually inspected.

### 🌱 New step: greenhouse-gas emissions calculation, with special handling for new forest

- The script adds an analysis step not described in the original README at all: estimating CO₂-equivalent emissions or removals from each land-use change, using an external carbon-stock dataset.
- This version adds **specific handling for new-forest and other-land conversions**: rather than one generic emissions factor, it uses each cell's own carbon-stock measurement, divided across a realistic number of recovery years (50 years for new forest, 20 for other land), and accumulated over time.

### 🛡️ Adjustment to protected-area restrictions

- Protected-area restrictions apply specifically to **forest-to-pasture and forest-to-cropland** conversion. Forest-to-urban conversion is not restricted.

### ✅ New step: built-in results verification

- The script adds **automatic consistency checks** after the run (do areas balance between periods? do gains and losses add up?) and prints a clear OK/FAILED summary. The original README has no equivalent quality-check step — previously, this would need to be checked manually.

### 📊 New step: summary table of total area by year

- A simple, easy-to-read table showing total area per land-use class for every year of the projection is generated at the end — a convenience addition not present in the original workflow.

### 📦 Extra package requirements not listed in the README

- The script relies on several packages the README's "Requirements" section does not mention, including `here`, `ggnewscale`, `ggpattern`, `rnaturalearth`, `rnaturalearthdata`, `countrycode`, `stringr`, `readxl`, `nloptr`, `arrow` (for the new pasture-condition dataset), and `RColorBrewer` (for the new comparison maps).

------------------------------------------------------------------------

## 3. Two new mapping scripts: `luc_change_map.R` and `ghg_emissions_map.R`

Alongside the main run script, two standalone R scripts turn its outputs into maps. Neither is mentioned in the original project README — both are new additions that pick up where `run_country_newforest.R` leaves off, using the files it saves to `Output/BRA/`.

### 🗺️ `luc_change_map.R` — Land-use change map, 2020–2050

**What it shows:** For each land-use class (forest, cropland, pasture, other land, urban), one map panel of Brazil showing which grid cells **gained** that land use and which **lost** it between 2020 and 2050 — red for loss, blue for gain — all classes side by side in a single grid image, with state and biome boundaries drawn on top.

**How it works, simply:** 1. Loads the downscaled land-use results (`...downscaled_LUC.rds`) and the matching grid raster (`...rasterized_layer.tif`), both produced by `run_country_newforest.R`. 2. Matches each grid cell's ID in the results table to its ID in the raster. 3. For each class, adds up the total area in that class in 2020 and again in 2050, per grid cell, and subtracts to get the net change. 4. Drops the temporary "new forest" category, since it is not of interest for this particular map. 5. Colors every class on the **same** red-to-blue scale, so the panels are directly comparable to one another. 6. Overlays Brazilian state borders (thin, light grey) and biome borders (thick, black) and saves one combined image, `luc_change_2020_2050.png`.

### 🌱 `ghg_emissions_map.R` — Cumulative GHG (biomass) emissions map, 2020–2050

**What it shows:** For each land-use class, a map of Brazil showing where land-use change led to a net **release** of greenhouse gases (red) versus a net **removal/sequestration** (blue/green), added up across every 5-year period from 2020 to 2050 — plus one extra panel summing across all classes together.

**How it works, simply:** 1. Loads the same downscaled results file as above, but this time uses its `GHG_biomass` column (the emissions/sequestration estimate calculated inside `run_country_newforest.R`). 2. Adds up `GHG_biomass` for every grid cell, land-use class, and time period, treating any missing values as zero. 3. Matches each cell to the raster the same way as the land-use-change map. 4. Sums each cell's GHG value across **all** periods (2020–2050), separately for each class, producing one "running total" map per class. 5. Adds a combined "Total (all classes)" panel summing across all classes. 6. Because cropland's emissions are far larger in magnitude than other classes, each class panel gets its **own** color scale (so small signals aren't washed out), while the "Total" panel uses a special compressed scale to make both the widespread small sequestration signal and the rare large emission hot-spots visible on the same map. 7. Overlays the same state and biome borders and saves one combined image, `ghg_emissions_cumulative_by_class.png`.

### ⚠️ Same portability caveat as the main script

Both mapping scripts share the same limitation noted above for `run_country_newforest.R`: all of their inputs — the downscaled results, the rasterized grid, and the Brazilian state/biome boundary shapefiles — are read from **hardcoded absolute paths on one person's computer** (`C:/Users/User/Desktop/...`), rather than from a shared, documented location. They are also **Brazil-specific** (state/biome boundary files, "BRA" in the file names), so as written they cannot be reused by another country team without editing the paths and swapping in that country's boundary files.

------------------------------------------------------------------------

## 4. What this means for other country teams

Because this script hardcodes Brazil-specific values, file paths, and target numbers, it currently functions more as a **working/debugging version of the Brazil run** than a drop-in replacement for the generic `R/run_country.R` described in the main README. Before another team can reuse it, it would need:

1.  The three hardcoded absolute file paths replaced with files placed in `Data/<COUNTRY>/`, following the standard structure.
2.  The hardcoded config file name, FABLE sheet names, and manual target totals turned into configurable parameters (e.g., read from the `.yml` config).
3.  The new GHG-emissions, new-forest, and verification steps documented as optional/standard parts of the pipeline, with their own input requirements (including the pasture-condition dataset) added to the "Data Preparation" section of the main README.
4.  A documented, teamwide decision on which of the three spreading assumptions (Section 1, step 13) should be the default for any country, rather than it being a hardcoded choice in this script.

------------------------------------------------------------------------

## 📩 Support

For questions about this script or the underlying workflow, please contact: António Santos: [acasantos\@ciencias.ulisboa.pt](mailto:acasantos@ciencias.ulisboa.pt){.email} Or open an issue on GitHub.
