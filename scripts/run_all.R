# ============================================================================
# run_all.R — Master script to run the full CH4 data filtering pipeline
#
# Sources all analysis scripts in order. Run from the repo root directory.
# ============================================================================

message("============================================================")
message("CH4 Data Filtering Pipeline")
message("Working directory: ", getwd())
message("Start time: ", Sys.time())
message("============================================================\n")

# 00: Setup (libraries, paths, constants, helpers)
source("scripts/00_setup.R")

# 01: Load and prepare all flux datasets
source("scripts/01_load_data.R")

if (raw_available) {
  # 02: Per-closure precision (second-difference sigma; Allan deviation for comparison)
  #     and campaign sigma, from the raw records
  source("scripts/02_allan_deviation.R")
  # 02b: goFlux re-fit of every stem closure (RESUME_02B=1 reuses a previous fit; README)
  source("scripts/02b_goflux_fluxes.R")
  # Refresh the tracked copies of their results used when the raw records are absent
  file.copy(file.path(output_dir, DERIVED_FILES), file.path(derived_dir, DERIVED_FILES), overwrite = TRUE)
} else {
  # Without the raw stem records (data/stem_raw/, not in the repository), 03 onward use the
  # saved results of 02 and 02b in data/derived_from_stem_raw/ (copied to outputs/ by 00_setup.R)
  message("Raw stem records not found in data/stem_raw/: skipping 02, 02b and 14; using ",
          derived_dir, "/")
}

# 03: MDF thresholds and below-MDF flags (stems and canopy)
source("scripts/03_mdf_computation.R")

# 07: Below-MDF treatment comparison
source("scripts/07_treatment_comparison.R")

# 08: Transformation comparison (raw vs log vs arcsinh)
source("scripts/08_transformation_comparison.R")

# 09: MDF detection rates (plot data) and MDF vs measurement duration (Fig. S5)
source("scripts/09_mdf_figures.R")

# 10: Seasonal confounding (filter exclusion by season and temperature, budgets)
source("scripts/10_seasonal_confounding.R")

# 13: Scaling uncertainty / variance decomposition
source("scripts/13_scaling_uncertainty.R")

# 14: Readings of the closures drawn in Figs. 2, 5 and S1 -> data/from_hf_stem_study/figure_closures/
#     (from the raw records after 02b; skipped without them, and 15/17 use the tracked CSVs)
source("scripts/14_figure_closures.R")

# 15: CO2 validity example closures
source("scripts/15_co2_validity_figure.R")

# 16: Regression p-value vs MDF detection (over-detection by fit statistics; Fig. 4 data)
source("scripts/16_pvalue_vs_mdf.R")

# 17: sigma estimator figures (Fig. 2, Figs. S1-S2)
source("scripts/17_sigma_figures.R")

# 18: Known-truth simulation (filters and treatments scored against the true flux)
source("scripts/18_known_truth_simulation.R")

# 20: Same chain across gases, surfaces and analyzers (generality table)
source("scripts/20_generality.R")

# 21: Logging interval: re-logged 1 Hz records
source("scripts/21_logging_interval.R")

# 22: What each transformation does to inference (Fig. 8c, Fig. S6)
source("scripts/22_transform_inference.R")

# 23: Closure-level results table (results/stem_closures_processed.csv)
source("scripts/23_export_closures.R")

# 19: Assemble main-text figures (paper style, scripts/style.R)
source("scripts/19_main_figures.R")

# 11: Copy/rename figures into figures/main_text and figures/SI
source("scripts/11_organize_figures.R")

message("\n============================================================")
message("Pipeline complete!")
message("End time: ", Sys.time())
message("============================================================")
