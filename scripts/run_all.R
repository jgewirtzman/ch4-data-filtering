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

# 02: Compute Allan deviation (per-measurement noise)
source("scripts/02_allan_deviation.R")

# 03: Compute MDF thresholds and quality filter flags
source("scripts/03_mdf_computation.R")

# 04: Canopy filter sensitivity analysis
source("scripts/04_canopy_filter_sensitivity.R")

# 05: Instrument comparison (ridgeline + LGR vs LI-7810)
source("scripts/05_instrument_comparison.R")

# 06: Seasonal / monthly breakdown
source("scripts/06_seasonal_analysis.R")

# 07: Below-MDF treatment comparison
source("scripts/07_treatment_comparison.R")

# 08: Transformation comparison (raw vs log vs arcsinh)
source("scripts/08_transformation_comparison.R")

message("\n============================================================")
message("Pipeline complete!")
message("End time: ", Sys.time())
message("============================================================")
