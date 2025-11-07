# InSAR 3D Surface Modeling of Mt. Etna

## Overview

This repository was developed as a **course project** for the *Geospatial Data Analysis* class (A.Y. 2024/2025) at Politecnico di Milano, supervised by **Prof. Giovanna Venuti** and **Prof. Mirko Reguzzoni**.

This repository contains MATLAB scripts for the 3D reconstruction and **spatial correlation analysis** of ground deformation fields over **Mount Etna (Italy)** using **Persistent Scatterer Interferometric SAR (PS-InSAR)** observations.

The workflow integrates **statistical spatial modeling**, surface reconstruction, and **covariance-based correlation analysis** to jointly estimate ground deformation and its uncertainty structure.

It combines two complementary components:

- **Spline interpolation** (via the geoSplinter package, distributed in Prof. Mirko Reguzzoni’s Geospatial Data Analysis course) to derive a continuous 3D deformation surface and to prepare residual fields for later statistical modeling, in `Etna_PS_InSAR_Spline_Interp.m`.
- **Least Squares Collocation (LSC)** a geodetic spatial estimation method based on covariance modeling and spatial correlation analysis — to model the covariance of the residual field, perform optimal filtering of PS-InSAR data, and predict deformation values on a regular grid together with uncertainty estimates, in `Etna_PS_InSAR_LSCollocation_Interp.m`.

## Repository structure

```
├── Etna_PS_InSAR_Spline_Interp.m         # Main spline-processing script
├── Etna_PS_InSAR_LSCollocation_Interp.m  # Main LSC-processing script
├── geoSplinter.m                         # Plotting utility for geoSplinter outputs
├── Matlab_Functions/                     # Helper functions for job creation & statistics
│   ├── lambdaSplines2D.m                 # Regularization parameter estimator for 2D splines
│   ├── jobFile_analysis2D.m              # geoSplinter job builder (2D)
│   ├── jobFile_execution.m               # Cross-platform geoSplinter runner
│   └── ...
├── SAR_Data/                             # Input PS data set and QGIS project
├── data_input/, data_output/, job/       # Created by the spline script to interface geoSplinter
└── Etna_LSC_final_results.mat            # Final collocation results (created after running LSC)
```

## Requirements

- MATLAB R2021a or newer.
- Mapping Toolbox (used for coordinate transformations through `projcrs`, `projinv`, and `deg2utm`).
- geoSplinter executable (`geoSplinter_analysis.exe` for Windows or `geoSplinter_analysis_macOS`), provided to the course by Prof. Mirko and included in the repository root.
- The helper functions under `Matlab_Functions/` must remain on the MATLAB path; the main scripts add them automatically.

## Input data

The folder `SAR_Data/` contains **Persistent Scatterer (PS) InSAR** products downloaded from the [European Ground Motion Service (EGMS)](https://egms.land.copernicus.eu)
 for the southeastern flank of Mount Etna (2019–2023).
The data represent mean vertical ground velocities derived from Copernicus Sentinel-1 interferometric SAR observations, provided under the Copernicus free and open data policy.

## Getting started

1. **Open MATLAB and set the repository as the current folder.**
2. **Run the spline interpolation stage:**

   ```matlab
   Etna_PS_InSAR_Spline_Interp
   ```

   This script ingests the PS data set, performs polynomial trend fitting, prepares geoSplinter input files in `data_input/`, triggers spline interpolation, and exports diagnostic plots and cleaned displacement fields required by the LSC workflow.

3. **Run the least squares collocation stage:**

   ```matlab
   Etna_PS_InSAR_LSCollocation_Interp
   ```

   Performs empirical covariance estimation, fitting of an exponential-cosine spatial correlation model, and collocation filtering + prediction, producing deformation grids and uncertainty maps.

Intermediate and final products (plots, intermediate `.mat` files, and the final collocation grids) are written to `data_output/` and the repository root. Inspect `Etna_LSC_final_results.mat` to reuse the modeled fields in downstream analyses.

## Methodological Background

This workflow implements 3D spatial correlation modeling and statistical surface reconstruction for InSAR-derived deformation fields.

- Deterministic component: a 2nd-order polynomial trend capturing large-scale displacement gradients.
- Stochastic component: a zero-mean random field with spatial correlation characterized by the empirical covariance function **𝐶𝑣(𝜏)**, estimated from PS-InSAR residuals.
- Covariance modeling: an exponential-cosine model describes how correlation decays with distance, revealing the deformation’s spatial continuity and anisotropy.
- Least Squares Collocation: computes the Best Linear Unbiased Prediction (BLUP) of deformation and its posterior variance, fusing deterministic and stochastic information into a unified 3D deformation model.

## geoSplinter integration helpers

The `Matlab_Functions/` directory offers utilities that automate the preparation of geoSplinter jobs:

- `jobFile_analysis2D.m` formats 2D spline jobs, including automatic padding for bicubic splines and directory-aware paths to input/output files.
- `jobFile_execution.m` selects the appropriate geoSplinter executable for macOS or Windows and runs the generated job file.
- `lambdaSplines2D.m` estimates a regularization parameter from the noise variance, grid spacing, and local derivatives between neighboring PS observations.

These helpers are invoked inside the main scripts; you can also call them manually to prototype different spline configurations.

## Outputs and diagnostics

- **Collocation results:** `Etna_LSC_final_results.mat` contains the collocated displacement field at observation points, the modeled grid, and uncertainty summaries for further analysis.
- **Spatial correlation analysis:** empirical covariance plots, modeled covariance curves, and correlation-length diagnostics that quantify the deformation field’s spatial dependence.
- **Visualization:** 3D scatter plots, spline and LSC surface comparisons, and posterior-standard-deviation maps.

## Citation and acknowledgments

- **Data source:**
European Ground Motion Service (EGMS), part of the Copernicus Land Monitoring Service operated by the European Environment Agency (EEA).
https://egms.land.copernicus.eu
Data © European Union, 2023 — contains modified Copernicus Sentinel data (2019–2023), processed by the EGMS service.

- **Methodology:**
The workflow integrates the geoSplinter package and the Least Squares Collocation (LSC) framework developed by Giovanna Venuti and Mirko Reguzzoni in Geospatial Data Analysis course, Politecnico di Milano.