# SpatialSPEI: Spatial Calculation of SPEI and Evapotranspiration in R

This repository provides optimized R functions for the spatial calculation of Potential Evapotranspiration (PET / ETo) and the Standardised Precipitation Evapotranspiration Index (SPEI) from raster-based time series.

It is specifically designed for spatial ecology, forest modeling, and large-scale environmental data science. The workflow leverages C++ parallel processing via the `terra` package and the spatiotemporal data cube architecture of the `stars` package to ensure safe execution.

## Dependencies

To use these functions, the following R packages must be installed:

```R
install.packages(c("terra", "stars", "SPEI"))
```

## Direct Usage (Sourcing)

There is no need to install a formal package. You can load the functions directly into your R session by running:

```R
source("https://raw.githubusercontent.com/pereGelabert/spatialSPEI/main/spei_funcs.R")
```

---

## Function Documentation

### 1. `calc_spatial_pet()`

Calculates Potential Evapotranspiration (Thornthwaite) or Reference Evapotranspiration (Classic Hargreaves / Modified Hargreaves) by applying the equations cell-by-cell in parallel.

**Arguments:**
*   `tmed`: `SpatRaster`. Mean monthly temperature (ºC). Required for Thornthwaite.
*   `tmin`: `SpatRaster`. Minimum monthly temperature (ºC). Required for Hargreaves.
*   `tmax`: `SpatRaster`. Maximum monthly temperature (ºC). Required for Hargreaves.
*   `prec`: `SpatRaster`. Monthly precipitation (mm). Optional. If provided alongside Hargreaves, the modified Droogers & Allen (2002) method is applied to adjust radiation calculations.
*   `method`: Character. Calculation method: `"thornthwaite"` or `"hargreaves"`.
*   `start_year`: Numeric. Starting year of the time series (essential for calculating extraterrestrial radiation).
*   `start_month`: Numeric. Starting month of the time series (1-12). Defaults to `1`.
*   `cores`: Numeric. Number of cores for matrix parallelization. Defaults to `1`.

**Returns:**
A `SpatRaster` object containing the monthly PET/ETo values (mm).

---

### 2. `calc_spatial_spei()`

Calculates the SPEI based on the climatic water balance (Precipitation - Demand). It utilizes a data cube architecture (`stars`) to guarantee the integrity of spatiotemporal dimensions and prevent matrix coercion errors.

**Arguments:**
*   `prec`: `SpatRaster`. Monthly precipitation (mm).
*   `pet`: `SpatRaster`. Monthly evapotranspiration (mm), obtained via `calc_spatial_pet()`.
*   `scale`: Numeric. The time scale of the SPEI accumulation (e.g., 3, 6, 12, 24 months).
*   `start_year`: Numeric. Starting year of the time series (ensures correct monthly grouping for probabilistic fitting).
*   `start_month`: Numeric. Starting month of the time series (1-12). Defaults to `1`.
*   `kernel`: List. Smoothing kernel parameters, e.g., `list(type = 'rectangular', shift = 0)`.
*   `distribution`: Character. Probability distribution to fit: `'log-Logistic'` (default), `'Gamma'`, or `'PearsonIII'`.
*   `fit`: Character. Parameter fitting method: `'ub-pwm'` (default), `'pp-pwm'`, or `'max-lik'`.

**Returns:**
A `SpatRaster` object with the same spatial dimensions containing the standardized SPEI anomaly values.

---

## Applied Workflow Example

The following example demonstrates the complete pipeline using simulated raster inputs. It can be directly adapted to real gridded climate datasets (e.g., ERA5-Land, WorldClim time-series).

```R
library(terra)
library(stars)
library(SPEI)

# Load the functions
source("https://raw.githubusercontent.com/pereGelabert/spatialSPEI/main/spei_funcs.R")

# 1. Load Data (Pre-aggregated monthly SpatRasters)
# r_prec <- rast("monthly_precipitation.tif")
# r_tmin <- rast("monthly_tmin.tif")
# r_tmax <- rast("monthly_tmax.tif")

# 2. Atmospheric Demand Calculation (ETo using Modified Hargreaves)
# Starting the time series at the hydrological year (October 2015)
monthly_eto <- calc_spatial_pet(
  tmin = r_tmin, 
  tmax = r_tmax, 
  prec = r_prec, 
  method = "hargreaves", 
  start_year = 2020, 
  start_month = 10, 
  cores = 4
)

# 3. SPEI Calculation at a 6-month scale (Agricultural/Ecological drought)
spei_6 <- calc_spatial_spei(
  prec = r_prec, 
  pet = monthly_eto, 
  scale = 6, 
  start_year = 2020, 
  start_month = 10,
  distribution = "log-Logistic"
)

# 4. Export the result for downstream spatial modeling (e.g., Random Forest via ranger)
writeRaster(spei_6, "SPEI_6_2020_2025.tif", overwrite = TRUE)
```

## Methodological Considerations

*   **Time Series Length:** Robust fitting of probability distributions (like the log-logistic using L-moments) ideally requires a historical climate baseline exceeding 30 years. In shorter time series (e.g., < 10 years), the fitting process may fail in extreme pixels, resulting in `NA` values. In such data-limited scenarios, evaluating alternative fitting methods (e.g., `fit = 'max-lik'`) is recommended.
*   **NA Handling:** The internal architecture is programmed to automatically ignore continuous `NA` cells (e.g., oceans or out-of-bounds areas) using `all(is.na(x))`. However, it strictly requires temporal continuity (no missing months) within valid terrestrial cells.
