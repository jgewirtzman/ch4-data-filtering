# ============================================================================
# 07_treatment_comparison.R — Below-MDF treatment comparison
#
# Three approaches for handling fluxes below detection:
#   1. REMOVE:      exclude below-MDF measurements entirely
#   2. SET TO ZERO: replace below-MDF fluxes with 0
#   3. KEEP AS-IS:  retain original values (unfiltered baseline)
#
# Stringent filtering causes a rightward shift (positive bias)
# by selectively removing near-zero fluxes. Setting to zero
# preserves sample size and avoids this bias.
#
# Tables (console): summary stats by filter x treatment (and x instrument),
#   and % change in mean/median relative to the unfiltered baseline.
# The treatment panels of Fig. 7 are drawn in 19_main_figures.R.
#
# Loads: outputs/03_mdf_results.RData
# Saves: outputs/07_treatment.RData
# ============================================================================

source("scripts/00_setup.R")

message("\n=== Loading MDF results from 03 ===")
load(file.path(output_dir, "03_mdf_results.RData"))

# ── Ensure inst_label exists ─────────────────────────────────────────────────
if (!"inst_label" %in% names(df)) {
  df <- df %>%
    mutate(
      inst_label = case_when(
        year == 2025 ~ "Precise analyzer",
        year < 2025  ~ "Noisier analyzer"
      )
    )
}

# ============================================================
# BELOW-MDF TREATMENT COMPARISON
# ============================================================

message("\n=== Below-MDF treatment comparison ===")

# Focus on MDF-based filters only (these have clear thresholds)
mdf_filters <- list(
  "Manufacturer MDF" = "CH4_below_MDF_manuf",
  "Campaign \u03c3 90%"     = "CH4_below_MDF_wass90",
  "Campaign \u03c3 95%"     = "CH4_below_MDF_wass95",
  "Campaign \u03c3 99%"     = "CH4_below_MDF_wass99",
  "Per-closure \u03c3 \u00d73t 90%" = "CH4_below_MDF_chr90",
  "Per-closure \u03c3 \u00d73t 95%" = "CH4_below_MDF_chr95",
  "Per-closure \u03c3 \u00d73t 99%" = "CH4_below_MDF_chr99"
)

# --- Build long data frame with three treatments per filter ---
treat_rows <- list()
k <- 0

for (filt_name in names(mdf_filters)) {
  flag_col <- mdf_filters[[filt_name]]
  below <- df[[flag_col]]
  below[is.na(below)] <- FALSE  # no MDF (no precision estimate): keep the measurement

  n_below <- sum(below)
  n_above <- sum(!below)

  # Treatment 1: Remove below-MDF
  k <- k + 1
  treat_rows[[k]] <- data.frame(
    filter    = filt_name,
    treatment = "Remove below MDF",
    CH4_flux  = df$CH4_flux_nmolpm2ps[!below],
    instrument = df$inst_label[!below],
    stringsAsFactors = FALSE
  )

  # Treatment 2: Set below-MDF to zero
  flux_zeroed <- df$CH4_flux_nmolpm2ps
  flux_zeroed[below] <- 0
  k <- k + 1
  treat_rows[[k]] <- data.frame(
    filter    = filt_name,
    treatment = "Set below MDF to zero",
    CH4_flux  = flux_zeroed,
    instrument = df$inst_label,
    stringsAsFactors = FALSE
  )

  # Treatment 3: Keep original (unfiltered)
  k <- k + 1
  treat_rows[[k]] <- data.frame(
    filter    = filt_name,
    treatment = "Keep all (unfiltered)",
    CH4_flux  = df$CH4_flux_nmolpm2ps,
    instrument = df$inst_label,
    stringsAsFactors = FALSE
  )
}

treat_df <- do.call(rbind, treat_rows)

# Order filters from least to most stringent
filter_order_mdf <- c("Per-closure \u03c3 \u00d73t 99%", "Per-closure \u03c3 \u00d73t 95%", "Per-closure \u03c3 \u00d73t 90%",
                       "Campaign \u03c3 99%", "Campaign \u03c3 95%", "Campaign \u03c3 90%",
                       "Manufacturer MDF")
treat_df$filter <- factor(treat_df$filter, levels = filter_order_mdf)

# Treatment ordering: unfiltered as baseline, then the two alternatives
treat_df$treatment <- factor(treat_df$treatment,
                              levels = c("Keep all (unfiltered)",
                                         "Set below MDF to zero",
                                         "Remove below MDF"))

# --- TABLE: Summary stats by filter x treatment ---
cat("\n=== Summary stats: filter x treatment ===\n")
treat_stats <- treat_df %>%
  group_by(filter, treatment) %>%
  summarise(
    n          = n(),
    mean_flux  = round(mean(CH4_flux), 4),
    median_flux = round(median(CH4_flux), 4),
    pct_neg    = round(100 * mean(CH4_flux < 0), 1),
    pct_zero   = round(100 * mean(CH4_flux == 0), 1),
    .groups    = "drop"
  )
print(as.data.frame(treat_stats), row.names = FALSE)

# --- TABLE: Same split by instrument ---
cat("\n=== Summary stats: filter x treatment x instrument ===\n")
treat_inst_stats <- treat_df %>%
  group_by(filter, treatment, instrument) %>%
  summarise(
    n           = n(),
    mean_flux   = round(mean(CH4_flux), 4),
    median_flux = round(median(CH4_flux), 4),
    pct_neg     = round(100 * mean(CH4_flux < 0), 1),
    .groups     = "drop"
  )
print(as.data.frame(treat_inst_stats), row.names = FALSE)

# --- TABLE: % change in mean/median relative to the unfiltered baseline ---
baseline <- treat_stats %>%
  filter(treatment == "Keep all (unfiltered)") %>%
  select(filter, base_mean = mean_flux, base_median = median_flux)

pct_change <- treat_stats %>%
  filter(treatment != "Keep all (unfiltered)") %>%
  left_join(baseline, by = "filter") %>%
  mutate(
    pct_change_mean   = round(100 * (mean_flux - base_mean) / abs(base_mean), 1),
    pct_change_median = round(100 * (median_flux - base_median) / abs(base_median), 1)
  )

cat("\n=== % change in mean/median relative to unfiltered baseline ===\n")
print(as.data.frame(pct_change %>%
  select(filter, treatment, mean_flux, base_mean, pct_change_mean,
         median_flux, base_median, pct_change_median)), row.names = FALSE)

message("\n=== Treatment comparison complete ===")

# Plot-ready data for the assembled main figures (19_main_figures.R)
save(treat_stats, pct_change, treat_df, file = file.path(output_dir, "07_treatment.RData"))
