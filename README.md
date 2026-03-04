# ch4-data-filtering

Minimal detectable flux, quality filtering, and data transformation methods for chamber-based CH4 flux measurements. Companion code and data for a standalone methods paper.

## Overview

This repository evaluates how quality filtering, detection limit methodology, instrument precision, and data transformation choices affect CH4 flux estimates from closed-chamber measurements. The analysis uses three datasets:

- **Harvard Forest canopy fluxes** (136 measurements across stem, branch, and leaf surfaces; 3 ABB GLA131/UGGA instruments)
- **Harvard Forest stem fluxes** (1640 measurements; GLA131/UGGA 2023-24 + LI-COR LI-7810 2025)
- **Yale-Myers Forest wetland soil fluxes** (288 soil + 369 tree stem measurements; GLA131/UGGA)

## Key findings

1. **Empirical instrument precision differs from manufacturer specs by up to 40x** between instruments (GLA131 median Allan deviation 4.0 ppb vs LI-7810 0.095 ppb)
2. **Removing below-detection measurements inflates mean flux by 17-135%** depending on threshold stringency
3. **The right filtering approach depends on the scientific question** — retain all values for budgets/scaling, filter for detection/biology questions
4. **Log transforms silently drop 10-83% of CH4 data** depending on the proportion of negative (uptake) fluxes; arcsinh is preferred for variance stabilization

## Repository structure

```
scripts/           Analysis pipeline (run in order, or use run_all.R)
  00_setup.R         Libraries, paths, constants, helper functions
  01_load_data.R     Load all datasets
  02_allan_deviation.R  Per-measurement Allan deviation (all instruments)
  03_mdf_computation.R  MDF (3 methods) + quality filter definitions
  04_canopy_filter_sensitivity.R  HF canopy: component ridges, budget scaling
  05_instrument_comparison.R  UGGA vs LI-7810 precision & distributions
  06_seasonal_analysis.R  Monthly/seasonal negative flux patterns
  07_treatment_comparison.R  Remove vs zero vs keep-as-is bias analysis
  08_transformation_comparison.R  Raw vs log vs arcsinh
  run_all.R          Master script

data/              Input data (self-contained)
  canopy/            HF canopy compiled results + RData (manID, imp)
  stem/              HF stem flux CSV, chamber volumes, field logs
  stem_raw/          Raw 1-Hz instrument data (~650 MB, for Allan deviation)
    7810/              LI-7810 pre-segmented CSVs
    LGR/               LGR1/2/3 continuous recordings
  soil/              YMF wetland soil/stem fluxes
  ymf_canopy/        YMF Black Oak canopy fluxes

figures/           Output figures (PNG + PDF)
outputs/           Intermediate computed results (RData, CSV)
manuscript/        Writeup (.md, .qmd, .pdf, .docx)
```

## Usage

```r
# From the repo root:
source("scripts/run_all.R")

# Or run individual scripts in order:
source("scripts/00_setup.R")
source("scripts/01_load_data.R")
# ... etc
```

Scripts 03-08 can be re-run independently after scripts 01-02 have been run once (they load intermediate results from `outputs/`).

## Requirements

R packages: `dplyr`, `ggplot2`, `ggridges`, `tidyr`, `scales`, `patchwork`

Optional (for Allan deviation from raw data): `readxl` (for field log parsing)

## Data provenance

| Dataset | Source project | Instruments | N |
|---------|---------------|-------------|---|
| HF canopy | whole_tree_flux | GLA131 (LGR1, LGR2, LGR3) | 136 |
| HF stem | tree-flux-2025 | GLA131 (2023-24) + LI-7810 (2025) | 1640 |
| YMF soil+stem | tree-methanogens | GLA131 | 657 |
| YMF canopy | whole_tree_flux | GLA131 (LGR3) | 30 |
