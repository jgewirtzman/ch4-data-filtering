# ============================================================================
# 09_mdf_figures.R — MDF detection rate comparison and MDF vs measurement
#                    duration
#
#   1. MDF detection rates by method and instrument (plot data for
#      19_main_figures.R)
#   2. MDF vs measurement duration (theoretical curves; Fig. S5)
#   Summary statistics printed to the console.
#
# Loads: outputs/03_mdf_results.RData
# Saves: outputs/09_detection_rates.RData; figures/mdf_vs_duration (PNG + PDF)
# ============================================================================

source("scripts/00_setup.R")

message("\n=== Loading MDF results ===")
load(file.path(output_dir, "03_mdf_results.RData"))

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)


# ============================================================================
# MDF DETECTION RATES BY METHOD AND INSTRUMENT
# ============================================================================

message("\n=== MDF detection rates ===")

# --- Canopy detection rates ---

canopy_methods <- data.frame(
  method = c("Manufacturer\n(goFlux)",
             "Campaign \u03c3\n90%", "Campaign \u03c3\n95%", "Campaign \u03c3\n99%",
             "Per-closure \u03c3 \u00d73t\n90%", "Per-closure \u03c3 \u00d73t\n95%", "Per-closure \u03c3 \u00d73t\n99%"),
  flag_col = c("CH4_below_MDF_goflux",
               "CH4_below_MDF_wass90", "CH4_below_MDF_wass95", "CH4_below_MDF_wass99",
               "CH4_below_MDF_chr90", "CH4_below_MDF_chr95", "CH4_below_MDF_chr99"),
  approach = c("Manufacturer", rep("Campaign \u03c3 (z)", 3), rep("Per-closure \u03c3 \u00d73t", 3)),
  stringsAsFactors = FALSE
)

canopy_rates <- lapply(seq_len(nrow(canopy_methods)), function(i) {
  col <- canopy_methods$flag_col[i]
  vals <- hf_out[[col]]
  data.frame(
    method   = canopy_methods$method[i],
    approach = canopy_methods$approach[i],
    pct_below = 100 * mean(vals, na.rm = TRUE),
    n_below  = sum(vals, na.rm = TRUE),
    n_eval   = sum(!is.na(vals)),
    dataset  = "HF canopy",
    instrument = "Noisier analyzer (canopy)",
    stringsAsFactors = FALSE
  )
})
canopy_rate_df <- do.call(rbind, canopy_rates)

# --- Stem detection rates by instrument ---

stem_methods <- data.frame(
  method = c("Manufacturer",
             "Campaign \u03c3\n90%", "Campaign \u03c3\n95%", "Campaign \u03c3\n99%",
             "Per-closure \u03c3 \u00d73t\n90%", "Per-closure \u03c3 \u00d73t\n95%", "Per-closure \u03c3 \u00d73t\n99%"),
  flag_col = c("CH4_below_MDF_manuf",
               "CH4_below_MDF_wass90", "CH4_below_MDF_wass95", "CH4_below_MDF_wass99",
               "CH4_below_MDF_chr90", "CH4_below_MDF_chr95", "CH4_below_MDF_chr99"),
  approach = c("Manufacturer", rep("Campaign \u03c3 (z)", 3), rep("Per-closure \u03c3 \u00d73t", 3)),
  stringsAsFactors = FALSE
)

stem_rates <- lapply(seq_len(nrow(stem_methods)), function(i) {
  col <- stem_methods$flag_col[i]

  # GLA131
  lgr_vals <- df[[col]][df$year < 2025]
  lgr_row <- data.frame(
    method   = stem_methods$method[i],
    approach = stem_methods$approach[i],
    pct_below = 100 * mean(lgr_vals, na.rm = TRUE),
    n_below  = sum(lgr_vals, na.rm = TRUE),
    n_eval   = sum(!is.na(lgr_vals)),
    dataset  = "HF stem",
    instrument = "Noisier analyzer (stem)",
    stringsAsFactors = FALSE
  )

  # LI-7810
  li_vals <- df[[col]][df$year == 2025]
  li_row <- data.frame(
    method   = stem_methods$method[i],
    approach = stem_methods$approach[i],
    pct_below = 100 * mean(li_vals, na.rm = TRUE),
    n_below  = sum(li_vals, na.rm = TRUE),
    n_eval   = sum(!is.na(li_vals)),
    dataset  = "HF stem",
    instrument = "Precise analyzer (stem)",
    stringsAsFactors = FALSE
  )

  rbind(lgr_row, li_row)
})
stem_rate_df <- do.call(rbind, stem_rates)

# Combine
rate_df <- rbind(canopy_rate_df, stem_rate_df)

# Method ordering (increasing stringency)
method_levels <- c("Manufacturer\n(goFlux)", "Manufacturer",
                   "Campaign \u03c3\n90%", "Campaign \u03c3\n95%", "Campaign \u03c3\n99%",
                   "Per-closure \u03c3 \u00d73t\n90%", "Per-closure \u03c3 \u00d73t\n95%", "Per-closure \u03c3 \u00d73t\n99%")
rate_df$method <- factor(rate_df$method, levels = method_levels)

# Instrument ordering
rate_df$instrument <- factor(rate_df$instrument,
                              levels = c("Noisier analyzer (canopy)", "Noisier analyzer (stem)", "Precise analyzer (stem)"))

# MDF value columns (stem), summarized at the end
mdf_value_cols <- c("CH4_MDF_manufacturer",
                     "CH4_MDF_wass90", "CH4_MDF_wass95", "CH4_MDF_wass99",
                     "CH4_MDF_chr90", "CH4_MDF_chr95", "CH4_MDF_chr99")

# ============================================================================
# FIGURE: MDF VS MEASUREMENT DURATION
# ============================================================================

message("\n=== Figure: MDF vs measurement duration ===")

# Use median flux_term from stem data as representative chamber geometry
median_flux_term <- median(df$flux_term, na.rm = TRUE)
message("  Median flux_term: ", round(median_flux_term, 3))

# Precision values to compare: datasheet (1 s) vs the campaign sigma that
# enters the reference MDF (median over field days of the closure-median
# second-difference sigma, per analyzer)
precision_values <- data.frame(
  label = c(
    paste0("Noisier analyzer datasheet (", PREC_CH4_LGR, " ppb)"),
    paste0("Noisier analyzer campaign \u03c3 (", round(global_sd_ch4_lgr, 1), " ppb)"),
    paste0("Precise analyzer datasheet (", PREC_CH4_7810, " ppb)"),
    paste0("Precise analyzer campaign \u03c3 (", round(global_sd_ch4_7810, 3), " ppb)")
  ),
  precision = c(PREC_CH4_LGR, global_sd_ch4_lgr, PREC_CH4_7810, global_sd_ch4_7810),
  linetype = c("dashed", "solid", "dashed", "solid"),
  instrument = c("Noisier analyzer", "Noisier analyzer", "Precise analyzer", "Precise analyzer"),
  stringsAsFactors = FALSE
)

# Generate curves: MDF = Z95 * sigma / t × flux_term (reference form; the
# datasheet curve uses the same multiplier so only sigma differs)
t_seq <- seq(60, 600, by = 5)

mdf_curves <- lapply(seq_len(nrow(precision_values)), function(i) {
  data.frame(
    t_sec     = t_seq,
    mdf       = Z95 * precision_values$precision[i] / t_seq * median_flux_term,
    label     = precision_values$label[i],
    precision = precision_values$precision[i],
    linetype_grp = precision_values$linetype[i],
    instrument   = precision_values$instrument[i],
    stringsAsFactors = FALSE
  )
})
mdf_curve_df <- do.call(rbind, mdf_curves)

# Color by instrument, linetype by manufacturer vs empirical
mdf_curve_df$source <- ifelse(grepl("datasheet", mdf_curve_df$label),
                               "Datasheet (1 s spec)", "Campaign \u03c3")

p_mdf_duration <- ggplot(mdf_curve_df,
                          aes(x = t_sec, y = mdf, color = instrument,
                              linetype = source)) +
  geom_line(linewidth = 0.9) +
  geom_vline(xintercept = 300, linetype = "dotted", color = "grey50", linewidth = 0.4) +
  annotate("text", x = 310, y = max(mdf_curve_df$mdf) * 0.9,
           label = "5 min", hjust = 0, size = 3, color = "grey50") +
  scale_color_manual(
    values = c("Noisier analyzer" = COL$A, "Precise analyzer" = COL$B),
    name = "Analyzer"
  ) +
  scale_linetype_manual(
    values = c("Datasheet (1 s spec)" = "dashed", "Campaign \u03c3" = "solid"),
    labels = c("Datasheet (1 s spec)" = "Datasheet \u03c3", "Campaign \u03c3" = "Measured \u03c3"),
    name = NULL
  ) +
  scale_y_log10(
    labels = scales::label_number(drop0trailing = TRUE)
  ) +
  scale_x_continuous(breaks = seq(60, 600, by = 60),
                     labels = function(x) paste0(x / 60, " min")) +
  labs(
    x = "Measurement duration",
    y = expression(MDF~(nmol~m^{-2}~s^{-1}))
  ) +
  theme_paper(base_size = 9) +
  theme(
    legend.position = "inside", legend.position.inside = c(0.98, 0.98), legend.justification = c(1, 1),
    legend.background = element_rect(fill = scales::alpha("white", 0.85), colour = NA),
    legend.box = "vertical",
    legend.key.width = unit(28, "pt")
  ) +
  guides(color = guide_legend(order = 1),
         linetype = guide_legend(order = 2))

ggsave(file.path(fig_dir, "mdf_vs_duration.png"),
       p_mdf_duration, width = 9, height = 6, dpi = 300)
ggsave(file.path(fig_dir, "mdf_vs_duration.pdf"),
       p_mdf_duration, width = 9, height = 6)
message("Saved: mdf_vs_duration")

# ============================================================================
# PRINT SUMMARY
# ============================================================================

cat("\n=== MDF Summary Statistics ===\n\n")

cat("Canopy (HF, n=", nrow(hf_out), "):\n")
cat("  Manufacturer (goFlux):  ", round(100 * mean(hf_out$CH4_below_MDF_goflux, na.rm = TRUE), 1), "% below\n")
cat("  Campaign sigma 90-99%:  ",
    round(100 * mean(hf_out$CH4_below_MDF_wass90, na.rm = TRUE), 1), "-",
    round(100 * mean(hf_out$CH4_below_MDF_wass99, na.rm = TRUE), 1), "%\n")
cat("  Per-closure x3t 90-99%: ",
    round(100 * mean(hf_out$CH4_below_MDF_chr90, na.rm = TRUE), 1), "-",
    round(100 * mean(hf_out$CH4_below_MDF_chr99, na.rm = TRUE), 1), "%\n")

cat("\nStem (n=", nrow(df), "):\n")
for (inst in c("Noisier analyzer", "Precise analyzer")) {
  mask <- if (inst == "Noisier analyzer") df$year < 2025 else df$year == 2025
  n_inst <- sum(mask)
  cat("  ", inst, " (n=", n_inst, "):\n", sep = "")
  cat("    Manufacturer:         ", round(100 * mean(df$CH4_below_MDF_manuf[mask], na.rm = TRUE), 1), "%\n")
  cat("    Campaign sigma 90-99%:",
      round(100 * mean(df$CH4_below_MDF_wass90[mask], na.rm = TRUE), 1), "-",
      round(100 * mean(df$CH4_below_MDF_wass99[mask], na.rm = TRUE), 1), "%\n")
  cat("    Per-closure x3t 90-99%:",
      round(100 * mean(df$CH4_below_MDF_chr90[mask], na.rm = TRUE), 1), "-",
      round(100 * mean(df$CH4_below_MDF_chr99[mask], na.rm = TRUE), 1), "%\n")
}

cat("\nMDF range (stem, all methods):\n")
for (col in mdf_value_cols) {
  vals <- df[[col]]
  cat(sprintf("  %-25s min=%.4f  median=%.4f  max=%.4f\n",
    col, min(vals, na.rm = TRUE), median(vals, na.rm = TRUE), max(vals, na.rm = TRUE)))
}

cat("\nMedian flux_term used for duration curves:", round(median_flux_term, 3), "\n")
cat("MDF (95%) at t=300s, Noisier analyzer datasheet sigma:",
    round(Z95 * PREC_CH4_LGR / 300 * median_flux_term, 4), "nmol/m2/s\n")
cat("MDF (95%) at t=300s, Noisier analyzer campaign sigma:",
    round(Z95 * global_sd_ch4_lgr / 300 * median_flux_term, 4), "nmol/m2/s\n")
cat("MDF (95%) at t=300s, Precise analyzer campaign sigma:",
    round(Z95 * global_sd_ch4_7810 / 300 * median_flux_term, 6), "nmol/m2/s\n")

message("\n=== MDF figures complete ===")

# Plot-ready data for the assembled main figures (19_main_figures.R)
save(rate_df, file = file.path(output_dir, "09_detection_rates.RData"))
