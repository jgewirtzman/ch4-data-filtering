## ============================================================================
## 00_setup.R — Shared libraries, paths, constants, and helper functions
## ============================================================================

# Libraries
library(dplyr)
library(ggplot2)
library(ggridges)
library(tidyr)
library(scales)
library(patchwork)

# Paths (relative to repo root)
data_dir   <- "data"
fig_dir    <- "figures"
output_dir <- "outputs"

dir.create(fig_dir, showWarnings = FALSE)
dir.create(output_dir, showWarnings = FALSE)

# ── Instrument specifications ───────────────────────────────────────────────
# Manufacturer precision (1-sigma, 1 sec averaging)

# ABB GLA131-GGA (Microportable UGGA)
PREC_CH4_LGR  <- 0.9    # ppb
PREC_CO2_LGR  <- 0.35   # ppm

# LI-COR LI-7810
PREC_CH4_7810 <- 0.6    # ppb
PREC_CO2_7810 <- 3.5    # ppm

# ── Chamber geometry constants ──────────────────────────────────────────────
SURFAREA_M2    <- pi * 0.0508^2  # m^2, collar radius = 5.08 cm
R_hPa_L        <- 83.14472      # hPa*L/(mol*K), ideal gas constant
EXTRA_TUBE_VOL <- 0.028         # L, connecting tubing volume

# ── Helper functions ────────────────────────────────────────────────────────

#' Compute Allan deviation at tau = 1 s from a concentration timeseries.
#'
#' First-differencing removes any smooth trend (the flux signal), isolating
#' instrument white noise: sigma_Allan = SD(diff(x)) / sqrt(2).
allan_sd <- function(x) {
  x <- x[!is.na(x)]
  diffs <- diff(x)
  if (length(diffs) < 2) return(NA_real_)
  sd(diffs) / sqrt(2)
}

#' Compute Allan deviation per UniqueID from a manID-style data frame
#' (expects columns: UniqueID, flag, CO2dry_ppm, CH4dry_ppb)
compute_allan_per_uid <- function(manID_df) {
  manID_df %>%
    filter(flag == 1) %>%
    group_by(UniqueID) %>%
    summarise(
      allan_sd_CO2 = allan_sd(CO2dry_ppm),
      allan_sd_CH4 = allan_sd(CH4dry_ppb),
      n_meas_pts   = n(),
      .groups = "drop"
    )
}

message("Setup complete. Paths: data_dir=", data_dir,
        ", fig_dir=", fig_dir, ", output_dir=", output_dir)
