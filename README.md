# ch4-data-filtering

Code and data for:

> Gewirtzman, J. and Rosentreter, J. A.: Consequential choices and hidden bias: data processing for chamber greenhouse-gas fluxes near detection limits (in preparation).

The paper follows the processing chain after a chamber flux is fitted (instrument precision, detection limit, closure validity, quality filters, treatment of values below detection, and transformation) and measures how each choice changes reported means, medians, budgets, comparisons and tests, using CH₄ fluxes from tree stems, an upland–wetland soil gradient and an agricultural soil, N₂O from the agricultural soil, and simulations with known true fluxes.

## Repository layout

```
scripts/        Analysis pipeline (run_all.R runs it in order)
analysis/       Precision-estimator comparison on simulated records (feeds Fig. S2)
results/        Closure-level results for the stems (stem_closures_processed.csv; columns in results/README.md)
data/           Inputs: closure-level fluxes, chamber geometry, field logs, processed traces
figures/        Created by the pipeline (not tracked)
outputs/        Intermediate results created by the pipeline (not tracked)
```

## Pipeline

Run from the repository root:

```r
source("scripts/run_all.R")
```

| Script | What it does | Feeds |
|---|---|---|
| `00_setup.R` | Packages, paths, constants (analyzer specifications, chamber volumes), helpers | all |
| `01_load_data.R` | Loads stem, soil and canopy datasets | 02– |
| `02_allan_deviation.R` | Per-closure precision from the raw 1 Hz records (second-difference estimator), logging interval, closure length, slope p-value and residual autocorrelation | 03, 17 |
| `02b_goflux_fluxes.R` | Refits every stem closure with a raw trace using goFlux (linear and Hutchinson–Mosier, best-flux selection), detection class and trace screens | 03, 16, 17, 21 |
| `03_mdf_computation.R` | MDFs (datasheet; campaign σ at z = 1.645, 1.96, 2.576) and filters | 07–22 |
| `07_treatment_comparison.R` | Remove / set to zero / keep below-MDF closures | Fig. 7a |
| `08_transformation_comparison.R` | Distributions under log₁₀ and arcsinh | Fig. 8a |
| `09_mdf_figures.R` | Detection rates by MDF method; MDF against closure length | Fig. 3b, Fig. S5 |
| `10_seasonal_confounding.R` | Exclusion by temperature and season; cumulative budgets | Figs. 6, 7c–d |
| `13_scaling_uncertainty.R` | Uncertainty of the mean under each treatment | text (confidence intervals) |
| `15_co2_validity_figure.R` | CO₂-tracer example closures | Fig. 5 |
| `16_pvalue_vs_mdf.R` | Regression-based detection against the MDF | Fig. 4 |
| `17_sigma_figures.R` | How σ is estimated; worked example; estimator comparison | Fig. 2, Figs. S1–S2 |
| `18_known_truth_simulation.R` | Simulated field years with known fluxes; detection calibration; noise sweep | Figs. 7b, 9, S4 |
| `20_generality.R` | The same quantities for every dataset | Fig. 9, Table S4 |
| `21_logging_interval.R` | Re-logged 1 Hz records at 2, 5 and 10 s | Fig. S3 |
| `22_transform_inference.R` | What each transformation does to tests, on the stems and in simulation | Fig. 8c, Fig. S6 |
| `23_export_closures.R` | Closure-level results table for the stems | `results/` |
| `19_main_figures.R` | Assembles the main-text figures (style in `style.R`) | Figs. 1, 3, 4, 6–9 |
| `11_organize_figures.R` | Copies figures to `figures/main_text/` and `figures/SI/` with paper numbering | |

Scripts 03 onward read intermediate results from `outputs/` and can be re-run individually once 01–02b have run.

## Data

The small inputs are copied here so that scripts 01–23 run from the repository; the authoritative versions are the sources listed.

| Folder | Contents | Source |
|---|---|---|
| `data/stem/` | Harvard Forest stem fluxes 2023–2025 (source table), chamber volumes, field logs | J. Gewirtzman et al., in prep.; to be archived with the raw records (below) |
| `data/stem_raw/` | Raw 1 Hz analyzer records (≈650 MB; **not in the repository**) | Harvard Forest Data Archive package of the stem study (DOI to follow) |
| `data/canopy/`, `data/ymf_canopy/` | Canopy closures (stems, branches, leaves) and processed traces, used for the precision of two analyzer units | Gewirtzman et al. (in review) |
| `data/soil/` | Yale-Myers Forest upland–wetland soil fluxes | Gewirtzman (2026), doi:10.5281/zenodo.18779715 |
| `data/external/` | Raw soil closure windows (75 MB; **not in the repository**) | same Zenodo record |
| `data/agsoil/` | De-identified agricultural-soil CH₄ and N₂O closures (**not in the repository** until their source study is published; see `data/agsoil/README.md`) | J. Gewirtzman et al., in prep. |

To rerun from the raw records, place the analyzer files in `data/stem_raw/` (layout in `data/stem_raw/README.md`) and the soil closure windows (`lgr_manual_identification_results_soil.csv`) in `data/external/`, or set the environment variable `SOIL_WINDOWS` to that file. Scripts 03 onward run from the intermediate results of 01–02b. The paper's closure-level results for the stems are in `results/stem_closures_processed.csv`.

## Requirements

R 4.4 with: dplyr, tidyr, tibble, ggplot2, patchwork, scales, ggrepel, ggh4x, ggtext, ggridges, data.table, lme4, metafor, readxl, jsonlite, magick, ragg, and the goFlux fork **v0.5.0.9002** (goFlux, Rheault et al. 2024, with additions for empirical precision, detection classes and trace screens; doi:10.5281/zenodo.23256675). `.Rprofile` adds a project-only library so that other projects keep their own goFlux:

```r
remotes::install_github("jgewirtzman/goFlux@v0.5.0.9002",
                        lib = "~/Library/R/ch4-data-filtering-lib", dependencies = FALSE)
```

Versions used for the paper: R 4.4.3; dplyr 1.2.0, tidyr 1.3.2, ggplot2 4.0.2, patchwork 1.3.2, scales 1.4.0, ggrepel 0.9.6, ggh4x 0.3.1, ggtext 0.1.2, ggridges 0.5.6, data.table 1.18.2.1, lme4 1.1.35.5, metafor 4.8.0, readxl 1.4.5, magick 2.9.0.

## Licence

Code: MIT (`LICENSE`). Data tables and results in `data/` and `results/`: CC BY 4.0 (`LICENSE-data`).
