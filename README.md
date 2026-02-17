# 🌍 FABLE Spatial Downscaling Workflow

This repository contains the workflow used to run the spatial downscaling of land-use scenarios from the FABLE Calculator using the [FABLEDownscalR](https://github.com/FABLE-consortium/FABLEDownscalR) package.

It is designed so that country teams can:
  * Run the full downscaling for their country
  * Use their own national spatial data
  * Produce harmonized, spatially explicit land-use projections
    
No modification of the package code is required.


## 📦 Related Repository

This workflow relies on the companion package:

👉 FABLEDownscalR
https://github.com/FABLE-consortium/FABLEDownscalR

You do not need to edit this package.


## ✅ Requirements

Before starting, make sure you have:
  * R (≥ 4.2)
    https://cran.r-project.org
  * RStudio (recommended)
    https://posit.co/download/rstudio-desktop/

## 🚀 Installation

1️⃣ Download the Workflow

In RStudio, run:
```r
install.packages("usethis")
usethis::create_from_github("FABLE-consortium/DownscalingFABLE")
```
Or download directly from GitHub.

2️⃣ Install the Package

Install the analysis package once:
```r
install.packages("remotes")
remotes::install_github("FABLE-consortium/FABLEDownscalR")
```
## 📁 Data Preparation

Each country team must prepare a standardized data folder before running the downscaling workflow.

All inputs must follow the structure below so that the FABLEDownscalR package can load them automatically.

1️⃣ Folder Structure

Inside the repository, create:  ```Data/<COUNTRY_CODE>/```

Example:  ```Data/IND/```

2️⃣ Required Files

Each country folder must contain:

| File                                    | Description                      |
| --------------------------------------- | -------------------------------- |
| `FABLE.xlsx`                            | FABLE Calculator outputs         |
| `Population2020.geojson`                | Population data                  |
| `ProtectedAreas.geojson`                | Protected areas                  |
| `TravelTime.geojson`                    | Travel time                      |
| `LandCoverHILDA2015.geojson`            | Baseline land cover (HILDA)      |
| `LandCoverCopernicus2019.geojson`       | Baseline land cover (Copernicus) |
| `LandCoverChangeHILDA2015_2019.geojson` | Land-use transitions             |
| `ForestManagement.geojson`              | Forest management                |
| `Altitude.geojson`                      | Altitude                         |
| `Slope.geojson`                         | Slope                            |
| `GAEZCropDistribution2015.geojson`      | Crop yields                      |
| `GLW4WorldGriddedLivestock2020.geojson` | Livestock                        |

⚠️ Filenames must match exactly.

### Spatial files

All spatial files must contain a column called:
 id_c 
which uniquely identifies grid cells.

### FABLE Calculator Outputs (FABLE.xlsx)

Each country must provide national land-use projections.

Location: ```Data/<COUNTRY>/FABLE.xlsx```

Example: ````Data/IND/FABLE.xlsx````

#### Required Sheets
1) Baseline

Contains national baseline areas.

| LandCover     | value |
| ------------- | ----- |
| forest        | 25000 |
| cropland      | 18000 |
| pasture       | 16000 |
| urban         | 2000  |
| newforest     | 5000  |
| not relevant  | 4000  |

2) Pathway Sheets

One sheet per pathway:
````
CurrentTrends
Sustainable
HighAmbition
````

Each sheet must contain:

| Column        | Description        |
| ------------- | ------------------ |
| LandCoverInit | Initial land cover |
| YearStart     | Start year         |
| YearEnd       | End year           |
| ToForest      | → Forest           |
| ToOtherLand   | → Other land       |
| ToCropland    | → Cropland         |
| ToPasture     | → Pasture          |
| ToUrban       | → Urban            |
| ToNewForest   | → New forest       |

These are copied from:
````
FABLE-C → Sheet: 4_calc_land → Table: calc_landmatrix
````
## ⚙️ Configuration

1️⃣ Create a Config File

Go to: ```config/```

2️⃣Copy the template:
```
template.yml
```
3️⃣Rename it, for example:
```
IND.yml
```
4️⃣Edit the Config File

Example:
```r
country: "IND"
pathway: "CurrentTrends"
start_map_source: "HILDA"

data_root: "Data"
output_root: "Output"

seed: 1234
mnl_niter: 100
mnl_nburn: 50

crs_equal_area: 6933
pixel_res_m: 50000
```
Main Parameters
| Parameter          | Description                                |
| ------------------ | ------------------------------------------ |
| `country`          | ISO3 country code                          |
| `pathway`          | FABLE scenario                             |
| `start_map_source` | `HILDA` or `Copernicus`                    |
| `pixel_res_m`      | Raster resolution (50000 / 10000 / 1000 m) |
Dates and output tags are generated automatically.

## 📂 Setting the Working Directory 

Before running the workflow, R must know where your project folder is.
This is called setting the working directory.

All paths in the scripts assume that your working directory is the root of the DownscalingFABLE repository.

✅ Option 1 - Use RStudio Projects

Open RStudio

Click File → Open Project

Select:
```
DownscalingFABLE/DownscalingFABLE.Rproj
```

RStudio will automatically set the working directory to the project root.

You can check with:
```r
getwd()
```

It should show something like:
```
.../DownscalingFABLE
```
✅ Option 2 — Manually Set the Directory in R

If you are not using an .Rproj file:

Step 1: Find your project folder

```
Example:
C:/Users/Name/Documents/GitHub/DownscalingFABLE
```
Step 2: Set it in R
```r
setwd("C:/Users/Name/Documents/GitHub/DownscalingFABLE")
```

Then check:
```r
getwd()
```

## ▶️ Running the Workflow

Once the working directory is correctly set:

```r
source("scripts/run_country.R")
```
Replace IND.yml with your configuration file.

📤 Outputs
Results are saved in:
```
Output/<COUNTRY_CODE>/
```
Example:
```
Output/IND/
```
This folder contains:
  * Harmonized baseline maps
  * Rasterized ID layers
  * Model coefficients
  * Downscaled land-use projections
  * Intermediate datasets
These outputs are ready for mapping and analysis.

## ⚠️ Common Issues

1) **Missing Files**
    
```Error: file not found```
  
➡️ Check filenames in Data/\<country>/
  
2) **Missing id_c**

```Error: missing join key```

➡️ Ensure all GeoJSON files contain id_c

3) **Package Not Found**

```Error: FABLEDownscalR not installed```
  
➡️ Re-run installation step

4) **Model Fails to Converge**

```Warning: MNL did not converge    ```

➡️ Usually due to limited data for some land uses
  
➡️ The workflow will automatically apply regularization or skip unstable origins

## 📄 Documentation

Detailed methods are available in the package vignette (forthcoming)

The full technical pipeline is documented in the FABLEDownscalR repository

## 📩 Support

For questions or issues, please contact:
Clara Douzal: clara.douzal@unsdsn.org
Or open an issue on GitHub.
