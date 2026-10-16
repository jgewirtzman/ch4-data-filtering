## ============================================================================
## 00_setup.R — Shared libraries, paths, constants, and helper functions
## ============================================================================

# Locale: labels carry Greek/Unicode (sigma, R-squared); Rscript inherits a
# non-UTF-8 locale on some systems, which silently breaks string matching
if (!grepl("UTF-8", Sys.getlocale("LC_CTYPE"), ignore.case = TRUE))
  suppressWarnings(Sys.setlocale("LC_CTYPE", "en_US.UTF-8"))

# Libraries
library(dplyr)
library(ggplot2)
library(tidyr)
library(scales)
library(patchwork)

# Paths (relative to repo root)
data_dir   <- "data"
fig_dir    <- Sys.getenv("FIG_DIR", "figures")
output_dir <- "outputs"

# Inputs, one folder per source (data/README.md). They are minimal COPIES of other projects'
# data, reproduced only so that this analysis can be rerun: not new observations.
stem_dir    <- file.path(data_dir, "from_hf_stem_study")    # stem study (in prep.)
canopy_dir  <- file.path(data_dir, "from_canopy_study")     # canopy study (in review)
soil_dir    <- file.path(data_dir, "from_ymf_soil_zenodo")  # doi:10.5281/zenodo.18779715
agsoil_dir  <- file.path(data_dir, "from_agsoil_study")     # agricultural soil (in prep.)
raw_dir     <- file.path(data_dir, "stem_raw")              # raw stem records: local only
derived_dir <- file.path(data_dir, "derived_from_stem_raw") # saved results of 02 and 02b

# Without the raw stem records, 02, 02b and 14 cannot run; 03 onward then use the results of
# 02 and 02b saved in data/derived_from_stem_raw/, copied to outputs/ here when missing there.
raw_available <- dir.exists(file.path(raw_dir, "LGR", "LGR1")) && dir.exists(file.path(raw_dir, "7810"))
DERIVED_FILES <- c("02_allan_deviation.RData", "02b_goflux_fluxes.RData", "02b_manID.RData")

dir.create(fig_dir, showWarnings = FALSE)

# Figure style (theme_paper, role palettes) and the ggsave wrapper used by every script
source("scripts/style.R")
dir.create(output_dir, showWarnings = FALSE)
if (!raw_available) {
  for (f in DERIVED_FILES) if (!file.exists(file.path(output_dir, f))) {
    file.copy(file.path(derived_dir, f), file.path(output_dir, f))
    message("Raw stem records not found: using ", file.path(derived_dir, f))
  }
}

# ── Instrument specifications ───────────────────────────────────────────────
# Manufacturer precision (1-sigma, 1 sec averaging)

# ABB GLA131-GGA (Microportable UGGA)
PREC_CH4_LGR  <- 0.9    # ppb
PREC_CO2_LGR  <- 0.35   # ppm

# LI-COR LI-7810
PREC_CH4_7810 <- 0.6    # ppb
PREC_CO2_7810 <- 3.5    # ppm

# LI-COR LI-7820 (agricultural-soil N2O)
PREC_N2O_7820 <- 0.4    # ppb

# ── Detection-limit convention (as in Gewirtzman et al., in prep., and the goFlux fork) ──
# Reference MDF = Z95 * sigma_campaign / t * flux_term
#   Z95            two-sided 95 % normal quantile (1.96; the "2 sigma" criterion)
#   sigma_campaign median over the day's closures of the per-closure second-difference
#                  precision (sigma_closure / sigma_campaign below); the whole-record MAD of
#                  first differences (precision_mad_runs) is kept as a comparison only
#   t              closure length in SECONDS (not nb.obs: LGR2 logs every 10 s,
#                  the YMF oak record every 5 s)
Z90 <- qnorm(0.95); Z95 <- qnorm(0.975); Z99 <- qnorm(0.995)

# Display labels for the MDF variants (single source of truth for figures/tables).
# All variants share MDF = k * sigma / t * flux_term and differ only in where
# sigma comes from and what k is:
#   manuf   datasheet sigma (1 s spec), k = 1            (goFlux default)
#   wass**  campaign sigma (closure median of the second-difference sigma),
#           k = z (90/95/99 %)                           (reference = wass95)
#   chr**   per-closure Allan deviation, k = 3 * t_crit  (conservative comparison;
#           NOT Christiansen et al. 2015, whose MDF is precision / t with k = 1)
MDF_LABELS <- c(
  manuf  = "Datasheet \u03c3 (k = 1)",
  goflux = "Datasheet \u03c3 (k = 1)",
  wass90 = "Campaign \u03c3 90%", wass95 = "Campaign \u03c3 95%", wass99 = "Campaign \u03c3 99%",
  chr90  = "Per-closure \u03c3 \u00d73t 90%", chr95 = "Per-closure \u03c3 \u00d73t 95%",
  chr99  = "Per-closure \u03c3 \u00d73t 99%"
)

# ── Chamber geometry constants ──────────────────────────────────────────────
SURFAREA_M2    <- pi * 0.0508^2  # m^2, collar radius = 5.08 cm
R_hPa_L        <- 83.14472      # hPa*L/(mol*K), ideal gas constant
# Analyzer internal volume added to the collar volume (tree_volumes.csv already holds
# collar + cap + 28.9 cm3 of tubing); 0.028 L for both analyzers (LI-7810 total sample volume).
ANALYZER_VOL_L <- 0.028

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

# ── Precision, closure duration and detection: the goFlux fork is the single source ──
# goFlux (Rheault et al. 2024) version 0.5.0.9002 with additions (Gewirtzman 2026,
# doi:10.5281/zenodo.23256675; https://github.com/jgewirtzman/goFlux, tag v0.5.0.9002),
# installed in the project-only library set by .Rprofile so other projects keep their goFlux.
# The helpers below are thin wrappers that keep this repo's call sites; all estimation is
# done by the fork (empirical.prec, closure.time, spec.at.interval, flux.class, qc.flags).
if (!requireNamespace("goFlux", quietly = TRUE) ||
    utils::packageVersion("goFlux") < "0.5.0.9002")
  stop("This repo needs the goFlux fork >= 0.5.0.9002 in its project library:\n",
       "  remotes::install_github('jgewirtzman/goFlux@v0.5.0.9002',\n",
       "                          lib = '~/Library/R/ch4-data-filtering-lib', dependencies = FALSE)")

#' Precision of ONE closure (goFlux::empirical.prec, method "hadamard"): robust median-based
#' second-difference precision, 1.4826 MAD(Δ²x)/√6, differences never crossing gaps,
#' duplicate timestamps or interval changes (needs posix for that).
sigma_closure <- function(x, posix = NULL) {
  d <- data.frame(UniqueID = "c", flag = 1L, gas = x)
  if (!is.null(posix)) d$POSIX.time <- as.POSIXct(posix)
  r <- goFlux::empirical.prec(d, "gas", warn = FALSE)
  if (nrow(r) == 0) NA_real_ else r$prec[which.max(r$n.diff)]
}

#' Campaign precision: median of the per-closure values over the closures of one group
#' (goFlux::empirical.prec). Returns sigma and the number of closures.
sigma_campaign <- function(x, id, posix = NULL) {
  d <- data.frame(UniqueID = id, flag = 1L, gas = x)
  if (!is.null(posix)) d$POSIX.time <- as.POSIXct(posix)
  r <- goFlux::empirical.prec(d, "gas", warn = FALSE)
  k <- which.max(r$n.closures)
  c(sigma = r$prec[k], n_closures = r$n.closures[k])
}

#' Superseded whole-record first-difference MAD (goFlux::empirical.prec, method
#' "mad_diff1"), per run of constant logging interval; kept as a comparison only.
precision_mad_runs <- function(x, posix, ...) {
  d <- data.frame(gas = x, POSIX.time = as.POSIXct(posix))
  r <- suppressWarnings(goFlux::empirical.prec(d, "gas", method = "mad_diff1", warn = FALSE))
  if (nrow(r) == 0) return(data.frame(dt_s = NA_real_, n = 0L, mad_sigma = NA_real_))
  data.frame(dt_s = r$dt_s, n = r$n, mad_sigma = r$prec)
}

#' Logging interval (s) of a POSIXct vector: the median positive spacing
logging_interval <- function(posix) {
  d <- diff(as.numeric(posix))
  d <- d[is.finite(d) & d > 0]
  if (length(d) < 2) return(NA_real_)
  median(d)
}

#' Closure length in seconds per UniqueID of a flagged manID-style data frame:
#' goFlux::closure.time (span of the window plus one logging interval; never nb.obs).
closure_seconds <- function(manID_df) {
  manID_df %>%
    filter(flag == 1) %>%
    group_by(UniqueID) %>%
    summarise(
      dt_s  = logging_interval(POSIX.time),
      t_sec = goFlux::closure.time(as.numeric(POSIX.time) - as.numeric(min(POSIX.time))),
      n_obs = n(),
      .groups = "drop"
    )
}

#' Datasheet precision at a logging interval (goFlux::spec.at.interval); valid only when
#' each logged value is an average over the interval (checked for the 5 s and 10 s records).
spec_at_interval <- function(spec_1s, dt_s) goFlux::spec.at.interval(spec_1s, dt_s)

#' Ordinary least-squares slope test on a closure trace: p-value of the slope,
#' its SE (concentration units per second) and the lag-1 autocorrelation of the
#' residuals (which is what makes 1 Hz p-values over-detect).
slope_fit_stats <- function(x, t_s) {
  ok <- is.finite(x) & is.finite(t_s)
  x <- x[ok]; t_s <- t_s[ok]
  if (length(x) < 5) return(c(lm_p = NA_real_, lm_se = NA_real_, ac1 = NA_real_))
  fit <- lm(x ~ t_s)
  r <- residuals(fit)
  c(lm_p  = summary(fit)$coefficients[2, 4],
    lm_se = summary(fit)$coefficients[2, 2],
    ac1   = suppressWarnings(cor(r[-1], r[-length(r)])))
}

# ── Known-truth simulations (18_known_truth_simulation.R, 22_transform_inference.R) ──
# Monthly mean air temperature of a simulated field year: ~ -3 C (Jan) to 25 C (Jul)
SIM_TEMP_M <- 11 + 14 * sin(2 * pi * ((1:12) - 4.5) / 12)

#' True stem fluxes of n closures at air temperature `temp`: lognormal emission with a
#' Q10 = 2.5 seasonal median (0.12 nmol m-2 s-1 at 15 C, scaled by `mult`), and a share
#' `p_up` of small true uptake.
sim_true_flux <- function(n, temp, mult = 1, sdlog = 1.3, p_up = 0.10) {
  f <- rlnorm(n, log(0.12 * 2.5^((temp - 15) / 10)), sdlog) * mult
  up <- runif(n) < p_up
  f[up] <- -abs(rnorm(sum(up), 0.02, 0.01))
  f
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

#' The worked example of Supplement Sect. S1 (14_figure_closures.R, 17_sigma_figures.R):
#' two near-linear closures of the noisier analyzer (LGR1) from one well-sampled day, a flux
#' just below the MDF (about 0.6-0.95 x MDF, where the MDF and a regression disagree) and a
#' large flux. Near-linear = straight-line residual SD close to the closure's own sigma_i (no
#' saturation, no placement spike), so the figure shows the estimator, not a misfit.
#' `win`: the flagged (flag == 1) rows of outputs/02b_manID.RData, ordered by UniqueID and time.
#' Returns the day and one row each (low, high) of per-closure fit statistics.
pick_worked_example <- function(win) {
  lin <- win %>% group_by(UniqueID, analyzer, date_str) %>%
    summarise(t = list(as.numeric(POSIX.time - POSIX.time[1], units = "secs")), x = list(CH4dry_ppb), .groups = "drop") %>%
    mutate(n = lengths(x), sigma_i = mapply(function(x, t) sigma_closure(x, as.POSIXct(t, origin = "1970-01-01")), x, t),
           slope = mapply(function(t, x) coef(lm(x ~ t))[2], t, x),
           resid_sd = mapply(function(t, x) sd(residuals(lm(x ~ t))), t, x),
           jump = vapply(x, function(v) max(abs(diff(v))), numeric(1))) %>%
    filter(analyzer == "LGR1", n >= 170, n <= 240, resid_sd < 1.25 * sigma_i, jump < 12 * sigma_i)
  day_pick <- lin %>% group_by(date_str) %>%
    summarise(k = n(), lo = sum(between(abs(slope) * 180 / sigma_i, 1.2, 1.9)), hi = sum(abs(slope) * 180 / sigma_i > 25), .groups = "drop") %>%
    filter(lo > 0, hi > 0) %>% arrange(desc(k)) %>% slice(1)
  cand <- lin %>% filter(date_str == day_pick$date_str)
  list(day  = day_pick$date_str,
       low  = cand %>% filter(between(abs(slope) * 180 / sigma_i, 1.2, 1.9)) %>% arrange(resid_sd / sigma_i) %>% slice(1),
       high = cand %>% filter(abs(slope) * 180 / sigma_i > 25) %>% arrange(resid_sd / sigma_i) %>% slice(1))
}

# Small extracts of the raw records shown in the figures (written by 14_figure_closures.R,
# read by 15 and 17): tracked in the repository so those figures need no raw archive
fig_closures_dir <- file.path(stem_dir, "figure_closures")

message("Setup complete. Paths: data_dir=", data_dir,
        ", fig_dir=", fig_dir, ", output_dir=", output_dir)

# ============================================================================
# Raw stem records (shared by 02_allan_deviation.R and 02b_goflux_fluxes.R)
# Both loaders return goFlux-style frames: POSIX.time, CO2dry_ppm, CH4dry_ppb,
# H2O_ppm and the *_prec columns goFlux reads at fit time (datasheet values; they
# feed only goFlux's own MDF column, never the campaign-sigma MDF).
# ============================================================================

# Parse one LGR/UGGA raw file (txt or txt.zip); timestamps are the analyzer's
# SysTime (local clock, the clock the field logs' comp_* times refer to).
parse_lgr_file <- function(filepath) {
  tryCatch({
    if (grepl("\\.zip$", filepath)) {
      txt_name <- sub("\\.zip$", "", basename(filepath))
      con <- unz(filepath, txt_name)
      lines <- readLines(con, warn = FALSE)
      close(con)
    } else {
      lines <- readLines(filepath, warn = FALSE)
    }
    if (length(lines) < 3) return(NULL)
    dat <- tryCatch(
      read.csv(text = paste(lines[-1], collapse = "\n"),
               stringsAsFactors = FALSE, strip.white = TRUE),
      error = function(e) NULL)
    if (is.null(dat) || nrow(dat) == 0) return(NULL)
    ts_col  <- names(dat)[1]
    ch4_col <- grep("CH4.*d_ppm", names(dat), value = TRUE)[1]
    co2_col <- grep("CO2.*d_ppm", names(dat), value = TRUE)[1]
    h2o_col <- grep("^X\\.?H2O\\.?_ppm$|H2O\\._ppm$", names(dat), value = TRUE)[1]
    if (is.na(ch4_col) || is.na(co2_col)) return(NULL)
    out <- data.frame(
      # analyzer clock time, labelled UTC (all stem clock times share this label)
      POSIX.time = as.POSIXct(dat[[ts_col]], format = "%m/%d/%Y %H:%M:%OS", tz = "UTC"),
      CO2dry_ppm = as.numeric(dat[[co2_col]]),
      CH4dry_ppb = as.numeric(dat[[ch4_col]]) * 1000,
      H2O_ppm    = if (is.na(h2o_col)) 0 else as.numeric(dat[[h2o_col]]),
      stringsAsFactors = FALSE)
    out[!is.na(out$POSIX.time), ]
  }, error = function(e) NULL)
}

# All _f####.txt(.zip) data files of one LGR date folder, concatenated and sorted
load_lgr_day <- function(folder_path) {
  all_files  <- list.files(folder_path, full.names = TRUE)
  data_files <- all_files[grepl("_f\\d+\\.txt(\\.zip)?$", all_files)]   # "(1)" duplicate copies excluded
  if (length(data_files) == 0) return(NULL)
  day_data <- do.call(rbind, Filter(Negate(is.null), lapply(data_files, parse_lgr_file)))
  if (is.null(day_data) || nrow(day_data) == 0) return(NULL)
  day_data <- day_data[order(day_data$POSIX.time), ]
  day_data <- day_data[!duplicated(day_data$POSIX.time), ]   # a day archived in two folders
  day_data$DATE     <- format(day_data$POSIX.time, "%Y-%m-%d")
  day_data$CO2_prec <- PREC_CO2_LGR
  day_data$CH4_prec <- PREC_CH4_LGR
  day_data$H2O_prec <- 200
  rownames(day_data) <- NULL
  day_data
}

# One LI-7810 daily .data file through goFlux's importer, keeping REMARK
# (the operator's closure label) and the instrument's local DATE/TIME.
load_li7810_day <- function(filepath, tz = "UTC") {
  d <- goFlux::import.LI7810(filepath, date.format = "ymd", timezone = tz,
                             save = FALSE, keep_all = TRUE,
                             prec = c(PREC_CO2_7810, PREC_CH4_7810, 45))
  d$REMARK[is.na(d$REMARK)] <- ""
  d$REMARK <- trimws(as.character(d$REMARK))
  d <- d[order(d$POSIX.time), ]
  rownames(d) <- NULL
  d
}

# Contiguous runs of a non-empty REMARK = one closure each (the LI-7810 operator
# types the tree tag at chamber closure and clears it at opening). Returns one
# row per run with row indices and times; k numbers repeat visits of a tag on a day.
segment_remarks <- function(rec) {
  r <- rle(rec$REMARK)
  i1 <- cumsum(r$lengths); i0 <- i1 - r$lengths + 1
  keep <- r$values != ""
  if (!any(keep)) return(NULL)
  seg <- data.frame(REMARK = r$values[keep], i0 = i0[keep], i1 = i1[keep],
                    n = r$lengths[keep], stringsAsFactors = FALSE)
  seg$DATE  <- as.character(rec$DATE[seg$i0])
  seg$start <- rec$POSIX.time[seg$i0]
  seg$end   <- rec$POSIX.time[seg$i1]
  seg <- seg[order(seg$start), ]
  seg$k <- ave(seq_len(nrow(seg)), seg$REMARK, seg$DATE, FUN = seq_along)
  seg$seg_key <- paste(seg$REMARK, seg$DATE, seg$k, sep = "_")
  rownames(seg) <- NULL
  seg
}
