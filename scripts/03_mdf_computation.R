# ============================================================================
# 03_mdf_computation.R — Compute MDF thresholds and quality filter flags
#
# All variants share MDF = k * sigma / t * flux.term (t in SECONDS) and differ
# only in the sigma source and the multiplier k:
#   1. Datasheet:    sigma = manufacturer 1-s spec, k = 1        (goFlux default;
#                    the form of Christiansen et al. 2015)
#   2. Campaign:     sigma = MAD(first differences)/sqrt(2) over the whole record
#                    of an analyzer x field day, k = z (90/95/99 %; Wassmann et
#                    al. 2018 form). z = 1.96 (95 %) is the REFERENCE method
#                    (guidelines paper / goFlux fork; Parkin et al. 2012 corroborate).
#   3. Per-closure x3t: sigma = per-closure Allan deviation, k = 3 * t_crit(df).
#                    Kept as a conservative comparison. Column suffix "chr" is
#                    historical: this is NOT Christiansen et al. (2015), whose
#                    3*SD*t99 is a GC quantification limit, not a per-closure MDF.
#
# Loads: outputs/01_loaded_data.RData, outputs/02_allan_deviation.RData
# Saves: outputs/03_mdf_results.RData
# ============================================================================

source("scripts/00_setup.R")

message("\n=== Loading data from 01 and 02 ===")
load(file.path(output_dir, "01_loaded_data.RData"))
load(file.path(output_dir, "02_allan_deviation.RData"))

# goFlux-fork re-fits of the stem closures (02b) replace the legacy fluxes in df
f02b <- file.path(output_dir, "02b_goflux_fluxes.RData")
if (file.exists(f02b)) {
  load(f02b)
  message("Stem fluxes from 02b: ", sum(df$flux_src == "goFlux"), " goFlux best.flux, ",
          sum(df$flux_src != "goFlux"), " legacy (no archived trace)")
  # closures without an archived raw record cannot be refitted or given their own
  # precision; they are dropped so that every stem closure is processed the same way
  df <- df %>% filter(flux_src == "goFlux")
  message("Kept the ", nrow(df), " stem closures with an archived raw record")
} else {
  message("outputs/02b_goflux_fluxes.RData not found: using the legacy stem fluxes")
}

# Load YMF imp object if saved separately
ymf_imp_file <- file.path(output_dir, "01_ymf_imp.RData")
if (file.exists(ymf_imp_file)) {
  load(ymf_imp_file)
}

# ============================================================================
# PART 1: STEM FLUX — COMPUTE FLUX TERM, NOISE, MDF THRESHOLDS
# ============================================================================

message("\n=== Computing stem MDF thresholds ===")

df <- df %>%
  mutate(
    # flux_term = nmol / surfarea
    flux_term = nmol / surfarea,

    # Closure length in seconds (from 02; span + logging interval); fall back
    # to 300 s if the raw record was not found
    t_sec_est = if_else(!is.na(t_sec), as.numeric(t_sec), 300),

    # Empirical noise floor in flux units (Allan SD / t * flux_term)
    CH4_noise_floor = ifelse(!is.na(allan_sd_CH4) & !is.na(flux_term),
                              allan_sd_CH4 / t_sec_est * flux_term,
                              NA_real_),
    CO2_noise_floor = ifelse(!is.na(allan_sd_CO2) & !is.na(flux_term),
                              allan_sd_CO2 / t_sec_est * flux_term,
                              NA_real_),

    # Empirical SNR (Allan deviation based)
    CH4_snr_allan = ifelse(!is.na(CH4_noise_floor) & CH4_noise_floor > 0,
                            abs(CH4_flux_nmolpm2ps) / CH4_noise_floor,
                            NA_real_),

    # Manufacturer precision (instrument-specific)
    prec_ch4 = if_else(year == 2025, PREC_CH4_7810, PREC_CH4_LGR),
    prec_co2 = if_else(year == 2025, PREC_CO2_7810, PREC_CO2_LGR),

    # Manufacturer MDF = precision / t * flux_term
    CH4_MDF_manufacturer = ifelse(!is.na(flux_term),
                                   prec_ch4 / t_sec_est * flux_term,
                                   NA_real_),

    # Datasheet spec rescaled to the record's logging interval (white noise);
    # the fair comparison for LGR2 (10 s) records
    prec_ch4_equiv = spec_at_interval(prec_ch4, dt_s),
    prec_co2_equiv = spec_at_interval(prec_co2, dt_s),

    # Per-closure x3t variant = allan_sd * 3 * t_crit / t * flux_term
    # (df from the number of samples in the fit, not seconds)
    df_meas = pmax(coalesce(n_pts, t_sec_est) - 2, 1),
    t99 = qt(0.995, df = df_meas),
    t95 = qt(0.975, df = df_meas),
    t90 = qt(0.95,  df = df_meas),

    CH4_MDF_chr99 = ifelse(!is.na(allan_sd_CH4) & !is.na(flux_term),
                            (allan_sd_CH4 * 3 * t99) / t_sec_est * flux_term,
                            NA_real_),
    CH4_MDF_chr95 = ifelse(!is.na(allan_sd_CH4) & !is.na(flux_term),
                            (allan_sd_CH4 * 3 * t95) / t_sec_est * flux_term,
                            NA_real_),
    CH4_MDF_chr90 = ifelse(!is.na(allan_sd_CH4) & !is.na(flux_term),
                            (allan_sd_CH4 * 3 * t90) / t_sec_est * flux_term,
                            NA_real_),

    # Below-MDF flags
    CH4_below_MDF_manuf = ifelse(!is.na(CH4_MDF_manufacturer),
                                   abs(CH4_flux_nmolpm2ps) < CH4_MDF_manufacturer,
                                   NA),
    CH4_below_MDF_chr99 = ifelse(!is.na(CH4_MDF_chr99),
                                   abs(CH4_flux_nmolpm2ps) < CH4_MDF_chr99,
                                   NA),
    CH4_below_MDF_chr95 = ifelse(!is.na(CH4_MDF_chr95),
                                   abs(CH4_flux_nmolpm2ps) < CH4_MDF_chr95,
                                   NA),
    CH4_below_MDF_chr90 = ifelse(!is.na(CH4_MDF_chr90),
                                   abs(CH4_flux_nmolpm2ps) < CH4_MDF_chr90,
                                   NA)
  )

# --- Campaign-sigma MDF (reference): z * sigma_camp / t * flux_term ---
# sigma_camp (from 02) = median over that analyzer x field day's closures of the
# per-closure second-difference precision. Analyzer-level summaries for reporting:
global_sd_ch4_lgr  <- median(df$sigma_camp_CH4[df$year < 2025], na.rm = TRUE)
global_sd_ch4_7810 <- median(df$sigma_camp_CH4[df$year == 2025], na.rm = TRUE)

message("Campaign CH4 precision (median over days of the closure-median sigma):")
message("  LGR:     ", round(global_sd_ch4_lgr, 4), " ppb")
message("  LI-7810: ", round(global_sd_ch4_7810, 4), " ppb")

z99 <- Z99; z95 <- Z95; z90 <- Z90

df <- df %>%
  mutate(
    # campaign sigma used in the z-form MDF (kept under the old name for
    # downstream scripts)
    global_sd_ch4 = sigma_camp_CH4,

    CH4_MDF_wass99 = ifelse(!is.na(global_sd_ch4) & !is.na(flux_term),
                             (z99 * global_sd_ch4) / t_sec_est * flux_term,
                             NA_real_),
    CH4_MDF_wass95 = ifelse(!is.na(global_sd_ch4) & !is.na(flux_term),
                             (z95 * global_sd_ch4) / t_sec_est * flux_term,
                             NA_real_),
    CH4_MDF_wass90 = ifelse(!is.na(global_sd_ch4) & !is.na(flux_term),
                             (z90 * global_sd_ch4) / t_sec_est * flux_term,
                             NA_real_),

    CH4_below_MDF_wass99 = ifelse(!is.na(CH4_MDF_wass99),
                                    abs(CH4_flux_nmolpm2ps) < CH4_MDF_wass99,
                                    NA),
    CH4_below_MDF_wass95 = ifelse(!is.na(CH4_MDF_wass95),
                                    abs(CH4_flux_nmolpm2ps) < CH4_MDF_wass95,
                                    NA),
    CH4_below_MDF_wass90 = ifelse(!is.na(CH4_MDF_wass90),
                                    abs(CH4_flux_nmolpm2ps) < CH4_MDF_wass90,
                                    NA),

    # Reference method and detection class (retain-and-flag; nothing deleted)
    CH4_MDF_ref       = CH4_MDF_wass95,
    CH4_below_MDF_ref = CH4_below_MDF_wass95,
    CH4_det_class     = case_when(
      is.na(CH4_MDF_ref)                    ~ NA_character_,
      CH4_flux_nmolpm2ps >  CH4_MDF_ref     ~ "emission",
      CH4_flux_nmolpm2ps < -CH4_MDF_ref     ~ "uptake",
      TRUE                                  ~ "below detection")
  )

# ============================================================================
# STEM MDF SUMMARY
# ============================================================================

message("\n=== Stem MDF Summary ===")
# Helper to report filter stats accounting for NAs
report_mdf <- function(flag_col, label) {
  vals <- df[[flag_col]]
  n_eval <- sum(!is.na(vals))
  n_below <- sum(vals, na.rm = TRUE)
  pct <- if (n_eval > 0) round(100 * n_below / n_eval, 1) else NA
  message(sprintf("  %-25s %d/%d evaluated, %d below (%.1f%%)",
                  label, n_eval, n_total, n_below,
                  ifelse(is.na(pct), 0, pct)))
}

report_mdf("CH4_below_MDF_manuf", "Manufacturer MDF")
for (m in c("wass90", "wass95", "wass99", "chr90", "chr95", "chr99"))
  report_mdf(paste0("CH4_below_MDF_", m), MDF_LABELS[[m]])
message("  Detection class (reference, campaign sigma 95%):")
print(table(df$inst_label <- ifelse(df$year == 2025, "Precise analyzer", "Noisier analyzer"),
            df$CH4_det_class, useNA = "ifany"))

# ============================================================================
# PART 2: DEFINE ALL STEM QUALITY FILTERS
# ============================================================================

r2_sym <- "\u00B2"

quality_filters <- list()
quality_filters[["No filter"]] <- function(d) rep(FALSE, nrow(d))

# Manufacturer MDF (instrument-specific)
quality_filters[["Manufacturer MDF"]] <- function(d) {
  ifelse(!is.na(d$CH4_below_MDF_manuf), d$CH4_below_MDF_manuf, FALSE)
}

# CO2-based filters
quality_filters[["CO2 flux > 0"]]                    <- function(d) d$CO2_flux_umolpm2ps <= 0
quality_filters[["CO2 SNR (SE) > 2"]]                <- function(d) d$CO2_snr_se <= 2
quality_filters[[paste0("CO2 R", r2_sym, " > 0.7")]] <- function(d) d$CO2_r2 <= 0.7
quality_filters[[paste0("CO2 R", r2_sym, " > 0.8")]] <- function(d) d$CO2_r2 <= 0.8

# CH4 SE-based SNR
quality_filters[["CH4 SNR (SE) > 2"]] <- function(d) d$CH4_snr_se <= 2
quality_filters[["CH4 SNR (SE) > 3"]] <- function(d) d$CH4_snr_se <= 3

# CH4 Allan deviation SNR
# For measurements without Allan deviation, they pass by default (not excluded)
quality_filters[["CH4 SNR (Allan) > 2"]] <- function(d) {
  ifelse(!is.na(d$CH4_snr_allan), d$CH4_snr_allan <= 2, FALSE)
}
quality_filters[["CH4 SNR (Allan) > 3"]] <- function(d) {
  ifelse(!is.na(d$CH4_snr_allan), d$CH4_snr_allan <= 3, FALSE)
}

# Campaign-sigma (z) MDF thresholds and the per-closure x3t variant; labels
# from MDF_LABELS in 00_setup.R
for (m in c("wass90", "wass95", "wass99", "chr90", "chr95", "chr99")) {
  local({
    col <- paste0("CH4_below_MDF_", m)
    quality_filters[[MDF_LABELS[[m]]]] <<- function(d) ifelse(!is.na(d[[col]]), d[[col]], FALSE)
  })
}

# CH4 R^2 thresholds
quality_filters[[paste0("CH4 R", r2_sym, " > 0.5")]] <- function(d) d$CH4_r2 <= 0.5
quality_filters[[paste0("CH4 R", r2_sym, " > 0.7")]] <- function(d) d$CH4_r2 <= 0.7
quality_filters[[paste0("CH4 R", r2_sym, " > 0.9")]] <- function(d) d$CH4_r2 <= 0.9

filter_names <- names(quality_filters)

# ============================================================================
# SORT FILTERS BY STRINGENCY
# ============================================================================

n_retained <- sapply(filter_names, function(fn) {
  fails <- quality_filters[[fn]](df)
  sum(!fails)
})
filter_order <- names(sort(n_retained, decreasing = TRUE))

# ============================================================================
# BUILD LONG DATA FOR STEM RIDGELINE
# ============================================================================

# Count how many measurements have Allan deviation for annotation
n_with_allan <- sum(!is.na(df$allan_sd_CH4))

pass_labels <- sapply(filter_order, function(fn) {
  fails <- quality_filters[[fn]](df)
  n_pass <- sum(!fails)
  pct <- round(100 * n_pass / n_total, 1)
  n_neg <- sum(df$CH4_flux_nmolpm2ps[!fails] < 0)
  pct_neg <- round(100 * n_neg / n_pass, 1)
  paste0(fn, "\n(", n_pass, "/", n_total, ", ", pct, "%",
         " | ", pct_neg, "% neg)")
})

rows_list <- lapply(seq_along(filter_order), function(i) {
  fn <- filter_order[i]
  fails <- quality_filters[[fn]](df)
  data.frame(
    filter   = pass_labels[i],
    CH4_flux = df$CH4_flux_nmolpm2ps,
    passes   = !fails,
    stringsAsFactors = FALSE
  )
})

ridge_df <- do.call(rbind, rows_list)
ridge_df$filter <- factor(ridge_df$filter, levels = rev(pass_labels))

# Per-filter summary stats
filter_stats <- ridge_df %>%
  filter(passes) %>%
  group_by(filter) %>%
  summarise(
    mean_flux   = mean(CH4_flux, na.rm = TRUE),
    median_flux = median(CH4_flux, na.rm = TRUE),
    n_neg       = sum(CH4_flux < 0),
    n_total     = n(),
    pct_neg     = round(100 * n_neg / n_total, 1),
    .groups     = "drop"
  )

cat("\n=== Filter summary (CH4 fluxes) ===\n")
print(as.data.frame(filter_stats), row.names = FALSE)

# Retention % for fill color
ridge_df <- ridge_df %>%
  group_by(filter) %>%
  mutate(pct_retained = 100 * sum(passes) / n()) %>%
  ungroup()

# ============================================================================
# PART 3: CANOPY FLUX — COMPUTE THREE MDF APPROACHES
# ============================================================================

message("\n=== Computing canopy MDF thresholds (three approaches) ===")

# GLA131 manufacturer precision (1 sigma, 1s)
prec_co2 <- PREC_CO2_LGR  # 0.35 ppm
prec_ch4 <- PREC_CH4_LGR  # 0.9 ppb

# Wassmann multipliers (z-scores for normal distribution)
z_mult <- c("99" = qnorm(0.995), "95" = qnorm(0.975), "90" = qnorm(0.95))

compute_mdf_all <- function(compiled, allan_df, rw_df, mad_df, tsec_df) {

  # Merge Allan deviation per measurement
  compiled <- merge(compiled, allan_df[, c("UniqueID", "allan_sd_CO2",
                                            "allan_sd_CH4", "n_meas_pts")],
                    by = "UniqueID", all.x = TRUE)

  # Merge closure length (s) and logging interval per measurement
  compiled <- merge(compiled, tsec_df[, c("UniqueID", "t_sec", "dt_s", "n_obs")],
                    by = "UniqueID", all.x = TRUE)

  # Merge rolling window precision per instrument (comparison only)
  compiled <- merge(compiled,
                    rw_df[, c("instrument", "flat_sd_CO2", "flat_sd_CH4")],
                    by = "instrument", all.x = TRUE)

  # Merge campaign MAD sigma per instrument record (enters the z-form MDF)
  compiled <- merge(compiled,
                    mad_df[, c("instrument", "mad_sd_CO2", "mad_sd_CH4")],
                    by = "instrument", all.x = TRUE)

  compiled <- compiled %>%
    mutate(
      # Closure length in seconds (nb.obs only equals seconds at 1 Hz; LGR2 logs
      # every 10 s and the YMF record every 5 s). Fall back to nb.obs * dt.
      t_sec = if_else(is.na(t_sec), CO2_nb.obs * coalesce(dt_s, 1), t_sec),
      prec_co2_equiv = spec_at_interval(prec_co2, coalesce(dt_s, 1)),
      prec_ch4_equiv = spec_at_interval(prec_ch4, coalesce(dt_s, 1)),

      # --- Approach 1: goFlux / manufacturer precision (no confidence level) ---
      CO2_MDF_goflux = prec_co2 / t_sec * CO2_flux.term,
      CH4_MDF_goflux = prec_ch4 / t_sec * CH4_flux.term,

      # --- Approach 2: campaign sigma (z * MAD sigma) at 99/95/90 % ---
      CO2_MDF_wass99 = (z_mult["99"] * mad_sd_CO2) / t_sec * CO2_flux.term,
      CH4_MDF_wass99 = (z_mult["99"] * mad_sd_CH4) / t_sec * CH4_flux.term,
      CO2_MDF_wass95 = (z_mult["95"] * mad_sd_CO2) / t_sec * CO2_flux.term,
      CH4_MDF_wass95 = (z_mult["95"] * mad_sd_CH4) / t_sec * CH4_flux.term,
      CO2_MDF_wass90 = (z_mult["90"] * mad_sd_CO2) / t_sec * CO2_flux.term,
      CH4_MDF_wass90 = (z_mult["90"] * mad_sd_CH4) / t_sec * CH4_flux.term,

      # --- Approach 3: per-closure x3t variant (Allan SD * 3 * t_alpha) ---
      df_meas = pmax(n_meas_pts - 2, 1),
      t99 = qt(0.995, df = df_meas),
      t95 = qt(0.975, df = df_meas),
      t90 = qt(0.95,  df = df_meas),
      CO2_MDF_chr99 = (allan_sd_CO2 * 3 * t99) / t_sec * CO2_flux.term,
      CH4_MDF_chr99 = (allan_sd_CH4 * 3 * t99) / t_sec * CH4_flux.term,
      CO2_MDF_chr95 = (allan_sd_CO2 * 3 * t95) / t_sec * CO2_flux.term,
      CH4_MDF_chr95 = (allan_sd_CH4 * 3 * t95) / t_sec * CH4_flux.term,
      CO2_MDF_chr90 = (allan_sd_CO2 * 3 * t90) / t_sec * CO2_flux.term,
      CH4_MDF_chr90 = (allan_sd_CH4 * 3 * t90) / t_sec * CH4_flux.term
    )

  # --- Flag below-detection fluxes ---
  compiled <- compiled %>%
    mutate(
      # CO2 flags -- all methods and confidence levels
      CO2_below_MDF_goflux = abs(CO2_best.flux) < CO2_MDF_goflux,
      CO2_below_MDF_wass99 = abs(CO2_best.flux) < CO2_MDF_wass99,
      CO2_below_MDF_wass95 = abs(CO2_best.flux) < CO2_MDF_wass95,
      CO2_below_MDF_wass90 = abs(CO2_best.flux) < CO2_MDF_wass90,
      CO2_below_MDF_chr99  = abs(CO2_best.flux) < CO2_MDF_chr99,
      CO2_below_MDF_chr95  = abs(CO2_best.flux) < CO2_MDF_chr95,
      CO2_below_MDF_chr90  = abs(CO2_best.flux) < CO2_MDF_chr90,
      # Summary flags use 99% (most conservative)
      CO2_below_any_MDF = CO2_below_MDF_goflux | CO2_below_MDF_wass99 |
                           CO2_below_MDF_chr99,
      CO2_below_all_MDF = CO2_below_MDF_goflux & CO2_below_MDF_wass99 &
                           CO2_below_MDF_chr99,

      # CH4 flags -- all methods and confidence levels
      CH4_below_MDF_goflux = abs(CH4_best.flux) < CH4_MDF_goflux,
      CH4_below_MDF_wass99 = abs(CH4_best.flux) < CH4_MDF_wass99,
      CH4_below_MDF_wass95 = abs(CH4_best.flux) < CH4_MDF_wass95,
      CH4_below_MDF_wass90 = abs(CH4_best.flux) < CH4_MDF_wass90,
      CH4_below_MDF_chr99  = abs(CH4_best.flux) < CH4_MDF_chr99,
      CH4_below_MDF_chr95  = abs(CH4_best.flux) < CH4_MDF_chr95,
      CH4_below_MDF_chr90  = abs(CH4_best.flux) < CH4_MDF_chr90,
      # Summary flags use 99% (most conservative)
      CH4_below_any_MDF = CH4_below_MDF_goflux | CH4_below_MDF_wass99 |
                           CH4_below_MDF_chr99,
      CH4_below_all_MDF = CH4_below_MDF_goflux & CH4_below_MDF_wass99 &
                           CH4_below_MDF_chr99,

      # Reference method (campaign sigma, 95 %) and detection class
      CH4_MDF_ref       = CH4_MDF_wass95,
      CH4_below_MDF_ref = CH4_below_MDF_wass95,
      CH4_det_class     = case_when(
        is.na(CH4_MDF_ref)               ~ NA_character_,
        CH4_best.flux >  CH4_MDF_ref     ~ "emission",
        CH4_best.flux < -CH4_MDF_ref     ~ "uptake",
        TRUE                             ~ "below detection")
    )

  # Clean up intermediate columns
  compiled <- compiled %>% select(-df_meas, -t99, -t95, -t90)

  compiled
}

# Apply to HF
hf_out <- compute_mdf_all(hf_compiled, allan_hf, rw_results, mad_results, tsec_hf)

# Apply to YMF
ymf_out <- compute_mdf_all(ymf_compiled, allan_ymf, rw_results, mad_results, tsec_ymf)

# --- Add best-model R-squared and SNR ---
# SNR uses per-measurement Allan deviation as empirical noise floor:
#   noise_floor_flux = allan_sd / t_sec * flux.term  (1 sigma in flux units)
#   SNR = |flux| / noise_floor_flux
add_quality_cols <- function(df) {
  df %>%
    mutate(
      CH4_best.r2 = ifelse(CH4_model == "LM", CH4_LM.r2, CH4_HM.r2),
      CO2_best.r2 = ifelse(CO2_model == "LM", CO2_LM.r2, CO2_HM.r2),
      CH4_noise_floor = allan_sd_CH4 / t_sec * CH4_flux.term,
      CH4_SNR     = abs(CH4_best.flux) / CH4_noise_floor
    )
}
hf_out  <- add_quality_cols(hf_out)
ymf_out <- add_quality_cols(ymf_out)

# ============================================================================
# CANOPY MDF SUMMARY
# ============================================================================

# Helper to print detection summary for a dataset
print_mdf_summary <- function(out, dataset_name) {
  message("\n=== MDF Summary: ", dataset_name, " ===")

  for (gas in c("CO2", "CH4")) {
    n_valid <- sum(!is.na(out[[paste0(gas, "_best.flux")]]))
    message("\n  ", gas, " (n=", n_valid, " measurements):")

    # MDF ranges
    for (method in c("goflux", "wass99", "wass95", "wass90",
                     "chr99", "chr95", "chr90")) {
      col <- paste0(gas, "_MDF_", method)
      if (col %in% names(out)) {
        label <- MDF_LABELS[[method]]
        message(sprintf("    MDF range (%-20s): %s - %s", label,
                round(min(out[[col]], na.rm = TRUE), 6),
                round(max(out[[col]], na.rm = TRUE), 6)))
      }
    }

    # Detection rates
    message("    --- % below detection ---")
    below_g <- paste0(gas, "_below_MDF_goflux")
    message(sprintf("    goFlux (manufacturer):    %d / %d (%.1f%%)",
            sum(out[[below_g]], na.rm = TRUE), n_valid,
            100 * mean(out[[below_g]], na.rm = TRUE)))

    for (conf in c("99", "95", "90")) {
      wass_col <- paste0(gas, "_below_MDF_wass", conf)
      chr_col  <- paste0(gas, "_below_MDF_chr", conf)
      if (wass_col %in% names(out)) {
        message(sprintf("    Campaign sigma %s%%:        %d / %d (%.1f%%)", conf,
                sum(out[[wass_col]], na.rm = TRUE), n_valid,
                100 * mean(out[[wass_col]], na.rm = TRUE)))
      }
      if (chr_col %in% names(out)) {
        message(sprintf("    Per-closure x3t %s%%:      %d / %d (%.1f%%)", conf,
                sum(out[[chr_col]], na.rm = TRUE), n_valid,
                100 * mean(out[[chr_col]], na.rm = TRUE)))
      }
    }

    below_any <- paste0(gas, "_below_any_MDF")
    below_all <- paste0(gas, "_below_all_MDF")
    message(sprintf("    Below ANY (99%%):          %d / %d (%.1f%%)",
            sum(out[[below_any]], na.rm = TRUE), n_valid,
            100 * mean(out[[below_any]], na.rm = TRUE)))
    message(sprintf("    Below ALL (99%%):          %d / %d (%.1f%%)",
            sum(out[[below_all]], na.rm = TRUE), n_valid,
            100 * mean(out[[below_all]], na.rm = TRUE)))
  }
}

print_mdf_summary(hf_out, "Harvard Forest")
print_mdf_summary(ymf_out, "Yale-Myers Forest")

# --- Save precision summary ---

precision_summary <- data.frame(
  source = c("GLA131 manufacturer (1\u03C3, 1s)",
             rw_results$instrument),
  CO2_precision_ppm = c(prec_co2, rw_results$flat_sd_CO2),
  CH4_precision_ppb = c(prec_ch4, rw_results$flat_sd_CH4),
  method = c("Datasheet",
             rep("Rolling window (bottom 5%)", nrow(rw_results))),
  stringsAsFactors = FALSE
)

# Campaign sigma (closure median of second-difference sigma; enters the reference MDF),
# the earlier whole-record MAD (comparison), and the datasheet
# value rescaled to each record's logging interval (white-noise assumption)
precision_summary <- rbind(
  precision_summary,
  data.frame(source = mad_results$instrument,
             CO2_precision_ppm = mad_results$mad_sd_CO2,
             CH4_precision_ppb = mad_results$mad_sd_CH4,
             method = "Campaign sigma (closure median)", stringsAsFactors = FALSE),
  data.frame(source = mad_results$instrument,
             CO2_precision_ppm = mad_results$wr_sd_CO2,
             CH4_precision_ppb = mad_results$wr_sd_CH4,
             method = "Whole-record MAD (comparison)", stringsAsFactors = FALSE),
  data.frame(source = mad_results$instrument,
             CO2_precision_ppm = spec_at_interval(prec_co2, mad_results$dt_s),
             CH4_precision_ppb = spec_at_interval(prec_ch4, mad_results$dt_s),
             method = paste0("Datasheet at logging interval (",
                             round(mad_results$dt_s), " s)"), stringsAsFactors = FALSE)
)

# Add median Allan deviation per instrument
allan_median_hf <- manID_all_hf %>%
  filter(flag == 1) %>%
  group_by(UniqueID) %>%
  summarise(
    allan_co2 = allan_sd(CO2dry_ppm),
    allan_ch4 = allan_sd(CH4dry_ppb),
    .groups = "drop"
  )

# Map UniqueID -> instrument
allan_median_hf$instrument <- NA_character_
allan_median_hf$instrument[allan_median_hf$UniqueID %in% lgr1_uids] <- "LGR1"
allan_median_hf$instrument[allan_median_hf$UniqueID %in% lgr2_uids] <- "LGR2"
allan_median_hf$instrument[allan_median_hf$UniqueID %in% lgr3_uids] <- "LGR3"

allan_by_inst <- allan_median_hf %>%
  group_by(instrument) %>%
  summarise(
    median_allan_CO2 = median(allan_co2, na.rm = TRUE),
    median_allan_CH4 = median(allan_ch4, na.rm = TRUE),
    .groups = "drop"
  )

# YMF Allan median
allan_ymf_median <- data.frame(
  instrument = "YMF",
  median_allan_CO2 = median(allan_ymf$allan_sd_CO2, na.rm = TRUE),
  median_allan_CH4 = median(allan_ymf$allan_sd_CH4, na.rm = TRUE)
)
allan_by_inst <- rbind(allan_by_inst, allan_ymf_median)

allan_summary <- data.frame(
  source = allan_by_inst$instrument,
  CO2_precision_ppm = allan_by_inst$median_allan_CO2,
  CH4_precision_ppb = allan_by_inst$median_allan_CH4,
  method = "Allan deviation (median)",
  stringsAsFactors = FALSE
)

precision_summary <- rbind(precision_summary, allan_summary)

# Rename for display (anonymized: GLA131 -> Noisier analyzer, individual units -> A-1/A-2/A-3)
precision_summary$source <- sub("^GLA131.*", "Noisier analyzer spec", precision_summary$source)
precision_summary$source <- sub("^LGR1$", "A-1", precision_summary$source)
precision_summary$source <- sub("^LGR2$", "A-2", precision_summary$source)
precision_summary$source <- sub("^LGR3$", "A-3", precision_summary$source)
precision_summary$source <- sub("^YMF$", "A-3 (YMF)", precision_summary$source)

write.csv(precision_summary,
          file.path(output_dir, "precision_comparison.csv"),
          row.names = FALSE)

message("\n--- Precision comparison table ---")
print(precision_summary)

# ============================================================================
# SAVE ALL RESULTS
# ============================================================================

message("\n=== Saving MDF results ===")

save(df, n_total, n_with_allan,
     quality_filters, filter_names, filter_order, pass_labels,
     ridge_df, filter_stats,
     global_sd_ch4_lgr, global_sd_ch4_7810,
     hf_out, ymf_out,
     mad_results, mad_days, mad_inst_median, tsec_hf, tsec_ymf,
     precision_summary, allan_by_inst,
     file = file.path(output_dir, "03_mdf_results.RData"))

message("Saved: ", file.path(output_dir, "03_mdf_results.RData"))
message("\n=== MDF computation complete ===")
