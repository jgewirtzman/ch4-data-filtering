# ch4-data-filtering

Code and data for:

> Gewirtzman, J. and Rosentreter, J. A.: Consequential choices and hidden bias: data processing for chamber greenhouse-gas fluxes near detection limits (in preparation).

The paper follows the processing chain after a chamber flux is fitted (instrument precision, detection limit, closure validity, quality filters, treatment of values below detection, and transformation) and measures how each choice changes reported means, medians, budgets, comparisons and tests, using CH₄ fluxes from tree stems, an upland–wetland soil gradient and an agricultural soil, N₂O from the agricultural soil, and simulations with known true fluxes.

## Repository layout

```
scripts/        Analysis pipeline (run_all.R runs it in order)
analysis/       Precision-estimator comparison on simulated records (feeds Fig. S2)
results/        Closure-level results for the stems (stem_closures_processed.csv; columns in results/README.md)
data/           Inputs: minimal copies of the source datasets (NOT new data; data/README.md)
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
| `02_allan_deviation.R` | Per-closure precision from the raw 1 Hz records (second-difference estimator), logging interval, closure length, slope p-value and residual autocorrelation (needs the raw records; otherwise its saved results in `data/derived_from_stem_raw/` are used) | 03, 17 |
| `02b_goflux_fluxes.R` | Refits every stem closure with a raw trace using goFlux (linear and Hutchinson–Mosier, best-flux selection), detection class and trace screens (needs the raw records; as for 02) | 03, 16, 17, 21 |
| `03_mdf_computation.R` | MDFs (datasheet; campaign σ at z = 1.645, 1.96, 2.576) and below-MDF flags | 07–23 |
| `07_treatment_comparison.R` | Remove / set to zero / keep below-MDF closures | Fig. 7a |
| `08_transformation_comparison.R` | Distributions under log₁₀ and arcsinh | Fig. 8a |
| `09_mdf_figures.R` | Detection rates by MDF method; MDF against closure length | Fig. 3b, Fig. S5 |
| `10_seasonal_confounding.R` | Exclusion by temperature and season; cumulative budgets | Figs. 6, 7c–d |
| `13_scaling_uncertainty.R` | Uncertainty of the mean under each treatment | text (confidence intervals) |
| `14_figure_closures.R` | Copies the raw readings of the closures drawn in Figs. 2, 5 and S1 to `data/from_hf_stem_study/figure_closures/` (skipped without the raw records) | 15, 17 |
| `15_co2_validity_figure.R` | CO₂-tracer example closures | Fig. 5 |
| `16_pvalue_vs_mdf.R` | Regression-based detection against the MDF | Fig. 4 |
| `17_sigma_figures.R` | How σ is estimated; worked example; estimator comparison | Fig. 2, Figs. S1–S2, Table S1 |
| `18_known_truth_simulation.R` | Simulated field years with known fluxes; detection calibration; noise sweep | Figs. 7b, 9, S4 |
| `20_generality.R` | The same quantities for every dataset; soil σ and MDF from the soil closure windows | Figs. 3a, 6d, 8a, 9, Table S4 |
| `21_logging_interval.R` | Re-logged 1 Hz records at 2, 5 and 10 s | Fig. S3 |
| `22_transform_inference.R` | What each transformation does to tests, on the stems and in simulation | Fig. 8c, Fig. S6 |
| `23_export_closures.R` | Closure-level results table for the stems | `results/` |
| `19_main_figures.R` | Assembles the main-text figures (style in `style.R`) | Figs. 1, 3, 4, 6–9 |
| `11_organize_figures.R` | Copies figures to `figures/main_text/` and `figures/SI/` with paper numbering | |

`run_all.R` runs end to end from this repository alone. The raw stem records (`data/stem_raw/`, ≈650 MB) are not included; without them, `run_all.R` skips 02, 02b and 14, and `00_setup.R` copies the saved results of 02 and 02b from `data/derived_from_stem_raw/` to `outputs/` (only those missing there). From a fresh clone:

```sh
Rscript scripts/run_all.R
```

This reproduces every figure in `figures/main_text/` and `figures/SI/` and `results/stem_closures_processed.csv` exactly (only 02b's review exports, `outputs/02b_*.csv` and the QC plot PDF, need the raw records). With the raw records in `data/stem_raw/`, `run_all.R` runs 02, 02b and 14 as well and refreshes the copies in `data/derived_from_stem_raw/` and `data/from_hf_stem_study/figure_closures/`.

Scripts 03 onward read intermediate results from `outputs/` and can be re-run individually once 01 has run (and 02–02b, or their saved results). 15 and 17 read the readings of the closures they draw from `data/from_hf_stem_study/figure_closures/`, not from the raw records. Figures are written to `figures/` (or to the folder in the environment variable `FIG_DIR`).

`analysis/precision_bakeoff/synthetic.R` scores the σ estimators on simulated records (Fig. S2); its results (`synthetic_results.rds`) are tracked, and 17 reads them. To regenerate them: `Rscript analysis/precision_bakeoff/synthetic.R`.

### Re-running 02b

`02b_goflux_fluxes.R` extracts every closure window from the raw records and refits it with goFlux, which takes several minutes. A full run saves the fits to `outputs/02b_checkpoint.RData` and the windows to `outputs/02b_manID.RData`; with the environment variable `RESUME_02B=1`, 02b reloads them instead of refitting and only redoes the merge, detection classes and exports (for example `RESUME_02B=1 Rscript scripts/run_all.R`). Unset it after changing anything that affects the windows or fits.

## Data

**The data in `data/` are copies, not new data.** Each `from_*` folder is a minimal copy of another project's dataset (only the files, rows and columns the scripts read), reproduced so that the analysis can be rerun from this repository alone; `derived_from_stem_raw/` holds results computed here from the stem study's raw records. Obtain the data from the original sources and cite them; do not count these files as additional or independent observations in compilations, syntheses or meta-analyses. Sources, contents and licences are in `data/README.md` and each folder's README.

| Folder | Contents | Source |
|---|---|---|
| `data/from_hf_stem_study/` | Harvard Forest stem closures 2023–2025 (closure table), chamber volumes, field-log closure times; `figure_closures/`: readings of the few closures drawn in Figs. 2, 5 and S1, extracted from the raw records by `14_figure_closures.R` | J. Gewirtzman et al., in prep.; to be archived in the Harvard Forest Data Archive (DOI to follow) |
| `data/stem_raw/` | Raw 1 Hz analyzer records of the stem study (≈650 MB; **not in the repository**; layout in its README) | same |
| `data/derived_from_stem_raw/` | Saved results of 02 and 02b, computed from the raw records, so that 03 onward run without them | this repository (processing of the stem study) |
| `data/from_canopy_study/` | Canopy closures (stems, branches, leaves) and processed traces of three analyzer units, Harvard Forest and Yale-Myers, used for analyzer precision | Gewirtzman et al., in review |
| `data/from_ymf_soil_zenodo/` | Yale-Myers Forest upland–wetland soil fluxes, goFlux fits and closure windows (window readings only) | Gewirtzman (2026), doi:10.5281/zenodo.18779715 |
| `data/from_agsoil_study/` | De-identified agricultural-soil CH₄ and N₂O closures | J. Gewirtzman et al., in prep. |

To recompute 02 and 02b from scratch, place the analyzer files in `data/stem_raw/` (layout in `data/stem_raw/README.md`). To run 20 on the complete soil closure windows of the Zenodo record rather than the copy here, set the environment variable `SOIL_WINDOWS` to that file (results are identical). The paper's closure-level results for the stems are in `results/stem_closures_processed.csv`.

## Requirements

R 4.4 with: dplyr, tidyr, tibble, ggplot2, patchwork, scales, ggrepel, ggh4x, ggtext, data.table, R.utils (lets data.table read `.gz`), lme4, metafor, jsonlite, ragg, png (optional, to compare figures), remotes (to install the fork), and the goFlux fork **v0.5.0.9002** (goFlux, Rheault et al. 2024, with additions for empirical precision, detection classes and trace screens; doi:10.5281/zenodo.23256675). `.Rprofile` adds a project-only library so that other projects keep their own goFlux:

```r
remotes::install_github("jgewirtzman/goFlux@v0.5.0.9002",
                        lib = "~/Library/R/ch4-data-filtering-lib", dependencies = FALSE)
```

The library path in `.Rprofile` and in the command above is macOS-style; on Linux or Windows any writable folder works (use the same path in both places).

Versions used for the paper: R 4.4.3; dplyr 1.2.0, tidyr 1.3.2, tibble 3.3.1, ggplot2 4.0.2, patchwork 1.3.2, scales 1.4.0, ggrepel 0.9.6, ggh4x 0.3.1, ggtext 0.1.2, data.table 1.18.2.1, R.utils 2.13.0, lme4 1.1.35.5, metafor 4.8.0, jsonlite 2.0.0, ragg 1.3.3, remotes 2.5.0.

## Licence

Code: MIT (`LICENSE`). Data: each copy in `data/` stays under its source's licence (all CC BY 4.0); derived results in `data/derived_from_stem_raw/` and `results/`: CC BY 4.0. See `LICENSE-data`, which also repeats that these are copies, not new data.
