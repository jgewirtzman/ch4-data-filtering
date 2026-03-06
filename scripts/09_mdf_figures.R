# ============================================================================
# 09_mdf_figures.R — MDF detection rate comparison, threshold distributions,
#                    and MDF vs measurement duration
#
# Three figures to support Section 3.2 (and Section 2.4 note) of the outline:
#   A. MDF detection rates by method and instrument
#   B. MDF threshold distributions by method and instrument
#   C. MDF vs measurement duration (theoretical curves)
#
# Loads: outputs/03_mdf_results.RData
# Saves: figures to fig_dir (PNG + PDF)
# ============================================================================

source("scripts/00_setup.R")

message("\n=== Loading MDF results ===")
load(file.path(output_dir, "03_mdf_results.RData"))

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)


# ============================================================================
# FIGURE A: MDF DETECTION RATES BY METHOD AND INSTRUMENT
# ============================================================================

message("\n=== Figure A: MDF detection rates ===")

# --- Canopy detection rates ---

canopy_methods <- data.frame(
  method = c("Manufacturer\n(goFlux)",
             "Wassmann\n90%", "Wassmann\n95%", "Wassmann\n99%",
             "Christiansen\n90%", "Christiansen\n95%", "Christiansen\n99%"),
  flag_col = c("CH4_below_MDF_goflux",
               "CH4_below_MDF_wass90", "CH4_below_MDF_wass95", "CH4_below_MDF_wass99",
               "CH4_below_MDF_chr90", "CH4_below_MDF_chr95", "CH4_below_MDF_chr99"),
  approach = c("Manufacturer", rep("Wassmann", 3), rep("Christiansen", 3)),
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
    instrument = "Analyzer A (canopy)",
    stringsAsFactors = FALSE
  )
})
canopy_rate_df <- do.call(rbind, canopy_rates)

# --- Stem detection rates by instrument ---

stem_methods <- data.frame(
  method = c("Manufacturer",
             "Wassmann\n90%", "Wassmann\n95%", "Wassmann\n99%",
             "Christiansen\n90%", "Christiansen\n95%", "Christiansen\n99%"),
  flag_col = c("CH4_below_MDF_manuf",
               "CH4_below_MDF_wass90", "CH4_below_MDF_wass95", "CH4_below_MDF_wass99",
               "CH4_below_MDF_chr90", "CH4_below_MDF_chr95", "CH4_below_MDF_chr99"),
  approach = c("Manufacturer", rep("Wassmann", 3), rep("Christiansen", 3)),
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
    instrument = "Analyzer A (stem)",
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
    instrument = "Analyzer B (stem)",
    stringsAsFactors = FALSE
  )

  rbind(lgr_row, li_row)
})
stem_rate_df <- do.call(rbind, stem_rates)

# Combine
rate_df <- rbind(canopy_rate_df, stem_rate_df)

# Method ordering (increasing stringency)
method_levels <- c("Manufacturer\n(goFlux)", "Manufacturer",
                   "Wassmann\n90%", "Wassmann\n95%", "Wassmann\n99%",
                   "Christiansen\n90%", "Christiansen\n95%", "Christiansen\n99%")
rate_df$method <- factor(rate_df$method, levels = method_levels)

# Instrument ordering
rate_df$instrument <- factor(rate_df$instrument,
                              levels = c("Analyzer A (canopy)", "Analyzer A (stem)", "Analyzer B (stem)"))

# Approach colors
approach_colors <- c("Manufacturer" = "#999999",
                      "Wassmann" = "#E69F00",
                      "Christiansen" = "#0072B2")

# --- Main figure: stem instruments only (most relevant for Section 3.2) ---

stem_only <- rate_df %>% filter(dataset == "HF stem")

# For the stem plot, unify manufacturer label
stem_only$method <- as.character(stem_only$method)
stem_only$method[stem_only$method == "Manufacturer\n(goFlux)"] <- "Manufacturer"
stem_method_levels <- c("Manufacturer",
                         "Wassmann\n90%", "Wassmann\n95%", "Wassmann\n99%",
                         "Christiansen\n90%", "Christiansen\n95%", "Christiansen\n99%")
stem_only$method <- factor(stem_only$method, levels = stem_method_levels)

p_rates_stem <- ggplot(stem_only,
                        aes(x = method, y = pct_below, fill = instrument)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.85) +
  geom_text(aes(label = paste0(round(pct_below, 0), "%")),
            position = position_dodge(width = 0.7),
            vjust = -0.5, size = 3) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  annotate("text", x = 7.4, y = 52, label = "50%", size = 2.8, color = "grey50") +
  scale_fill_manual(
    values = c("Analyzer A (stem)" = "#E69F00", "Analyzer B (stem)" = "#56B4E9"),
    name = "Instrument"
  ) +
  scale_y_continuous(limits = c(0, 90), expand = expansion(mult = c(0, 0.05))) +
  labs(
    x = "MDF method",
    y = "% of measurements below detection",
    title = expression(CH[4]~detection~rates~by~MDF~method~and~instrument),
    subtitle = paste0("Stem flux dataset: Analyzer A (n=", sum(df$year < 2025),
                      ") vs Analyzer B (n=", sum(df$year == 2025), ")")
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(size = 9, lineheight = 0.9),
    legend.position = "bottom",
    plot.title = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 10, color = "grey40")
  )

ggsave(file.path(fig_dir, "mdf_detection_rates_stem.png"),
       p_rates_stem, width = 10, height = 6, dpi = 300)
ggsave(file.path(fig_dir, "mdf_detection_rates_stem.pdf"),
       p_rates_stem, width = 10, height = 6)
message("Saved: mdf_detection_rates_stem")

# --- Combined figure: canopy + stem ---

# For combined, use canopy Manufacturer label as-is
p_rates_all <- ggplot(rate_df,
                       aes(x = method, y = pct_below, fill = instrument)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65, alpha = 0.85) +
  geom_text(aes(label = paste0(round(pct_below, 0), "%")),
            position = position_dodge(width = 0.75),
            vjust = -0.5, size = 2.5) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  scale_fill_manual(
    values = c("Analyzer A (canopy)" = "#D55E00",
               "Analyzer A (stem)" = "#E69F00",
               "Analyzer B (stem)" = "#56B4E9"),
    name = "Instrument / dataset"
  ) +
  scale_y_continuous(limits = c(0, 90), expand = expansion(mult = c(0, 0.05))) +
  labs(
    x = "MDF method",
    y = "% of measurements below detection",
    title = expression(CH[4]~detection~rates~across~MDF~methods~","~instruments~","~and~datasets),
    subtitle = paste0("Canopy (n=", nrow(hf_out), ") | Stem Analyzer A (n=", sum(df$year < 2025),
                      ") | Stem Analyzer B (n=", sum(df$year == 2025), ")")
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(size = 8, lineheight = 0.9),
    legend.position = "bottom",
    plot.title = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 9, color = "grey40")
  )

ggsave(file.path(fig_dir, "mdf_detection_rates_all.png"),
       p_rates_all, width = 12, height = 6, dpi = 300)
ggsave(file.path(fig_dir, "mdf_detection_rates_all.pdf"),
       p_rates_all, width = 12, height = 6)
message("Saved: mdf_detection_rates_all")


# ============================================================================
# FIGURE B: MDF THRESHOLD DISTRIBUTIONS BY METHOD AND INSTRUMENT
# ============================================================================

message("\n=== Figure B: MDF threshold distributions ===")

# Reshape MDF value columns to long format
mdf_value_cols <- c("CH4_MDF_manufacturer",
                     "CH4_MDF_wass90", "CH4_MDF_wass95", "CH4_MDF_wass99",
                     "CH4_MDF_chr90", "CH4_MDF_chr95", "CH4_MDF_chr99")
mdf_labels <- c("Manufacturer",
                 "Wassmann\n90%", "Wassmann\n95%", "Wassmann\n99%",
                 "Christiansen\n90%", "Christiansen\n95%", "Christiansen\n99%")

mdf_long <- lapply(seq_along(mdf_value_cols), function(i) {
  data.frame(
    method    = mdf_labels[i],
    mdf_value = df[[mdf_value_cols[i]]],
    instrument = ifelse(df$year == 2025, "Analyzer B", "Analyzer A"),
    stringsAsFactors = FALSE
  )
})
mdf_long_df <- do.call(rbind, mdf_long)
mdf_long_df <- mdf_long_df[!is.na(mdf_long_df$mdf_value), ]

mdf_long_df$method <- factor(mdf_long_df$method, levels = mdf_labels)
mdf_long_df$instrument <- factor(mdf_long_df$instrument,
                                  levels = c("Analyzer A", "Analyzer B"))

# Approach grouping for color
mdf_long_df$approach <- ifelse(grepl("Manufacturer", mdf_long_df$method), "Manufacturer",
                        ifelse(grepl("Wassmann", mdf_long_df$method), "Wassmann",
                               "Christiansen"))

p_mdf_dist <- ggplot(mdf_long_df,
                      aes(x = method, y = mdf_value, fill = instrument)) +
  geom_violin(alpha = 0.4, position = position_dodge(width = 0.8),
              draw_quantiles = c(0.5), scale = "width") +
  geom_boxplot(width = 0.15, position = position_dodge(width = 0.8),
               outlier.size = 0.5, alpha = 0.6) +
  scale_y_log10(
    labels = scales::label_number(drop0trailing = TRUE),
    breaks = c(0.001, 0.01, 0.1, 1, 10)
  ) +
  scale_fill_manual(
    values = c("Analyzer A" = "#E69F00", "Analyzer B" = "#56B4E9"),
    name = "Analyzer"
  ) +
  labs(
    x = "MDF method",
    y = expression(MDF~(nmol~m^{-2}~s^{-1})),
    title = expression(Distribution~of~per-measurement~CH[4]~MDF~thresholds),
    subtitle = paste0("Stem dataset (n=", nrow(df),
                      "). Violin width proportional to count; ",
                      "Christiansen requires Allan deviation (n=",
                      sum(!is.na(df$allan_sd_CH4)), ")")
  ) +
  annotation_logticks(sides = "l", size = 0.3) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(size = 9, lineheight = 0.9),
    legend.position = "bottom",
    plot.title = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 9, color = "grey40")
  )

ggsave(file.path(fig_dir, "mdf_threshold_distributions.png"),
       p_mdf_dist, width = 10, height = 6, dpi = 300)
ggsave(file.path(fig_dir, "mdf_threshold_distributions.pdf"),
       p_mdf_dist, width = 10, height = 6)
message("Saved: mdf_threshold_distributions")


# ============================================================================
# FIGURE C: MDF VS MEASUREMENT DURATION
# ============================================================================

message("\n=== Figure C: MDF vs measurement duration ===")

# Use median flux_term from stem data as representative chamber geometry
median_flux_term <- median(df$flux_term, na.rm = TRUE)
message("  Median flux_term: ", round(median_flux_term, 3))

# Precision values to compare
precision_values <- data.frame(
  label = c(
    paste0("Analyzer A manufacturer (", PREC_CH4_LGR, " ppb)"),
    paste0("Analyzer A empirical median (", round(median(df$allan_sd_CH4[df$year < 2025], na.rm = TRUE), 1), " ppb)"),
    paste0("Analyzer B manufacturer (", PREC_CH4_7810, " ppb)"),
    paste0("Analyzer B empirical median (", round(median(df$allan_sd_CH4[df$year == 2025], na.rm = TRUE), 3), " ppb)")
  ),
  precision = c(
    PREC_CH4_LGR,
    median(df$allan_sd_CH4[df$year < 2025], na.rm = TRUE),
    PREC_CH4_7810,
    median(df$allan_sd_CH4[df$year == 2025], na.rm = TRUE)
  ),
  linetype = c("dashed", "solid", "dashed", "solid"),
  instrument = c("Analyzer A", "Analyzer A", "Analyzer B", "Analyzer B"),
  stringsAsFactors = FALSE
)

# Generate curves: MDF = precision / t × flux_term
t_seq <- seq(60, 600, by = 5)

mdf_curves <- lapply(seq_len(nrow(precision_values)), function(i) {
  data.frame(
    t_sec     = t_seq,
    mdf       = precision_values$precision[i] / t_seq * median_flux_term,
    label     = precision_values$label[i],
    precision = precision_values$precision[i],
    linetype_grp = precision_values$linetype[i],
    instrument   = precision_values$instrument[i],
    stringsAsFactors = FALSE
  )
})
mdf_curve_df <- do.call(rbind, mdf_curves)

# Color by instrument, linetype by manufacturer vs empirical
mdf_curve_df$source <- ifelse(grepl("manufacturer", mdf_curve_df$label),
                               "Manufacturer spec", "Empirical (Allan)")

p_mdf_duration <- ggplot(mdf_curve_df,
                          aes(x = t_sec, y = mdf, color = instrument,
                              linetype = source)) +
  geom_line(linewidth = 0.9) +
  geom_vline(xintercept = 300, linetype = "dotted", color = "grey50", linewidth = 0.4) +
  annotate("text", x = 310, y = max(mdf_curve_df$mdf) * 0.9,
           label = "5 min", hjust = 0, size = 3, color = "grey50") +
  scale_color_manual(
    values = c("Analyzer A" = "#E69F00", "Analyzer B" = "#56B4E9"),
    name = "Analyzer"
  ) +
  scale_linetype_manual(
    values = c("Manufacturer spec" = "dashed", "Empirical (Allan)" = "solid"),
    name = "Precision source"
  ) +
  scale_y_continuous(
    labels = scales::label_number(drop0trailing = TRUE)
  ) +
  scale_x_continuous(breaks = seq(60, 600, by = 60),
                     labels = function(x) paste0(x / 60, " min")) +
  labs(
    x = "Measurement duration",
    y = expression(MDF~(nmol~m^{-2}~s^{-1})),
    title = expression(Minimal~detectable~CH[4]~flux~vs.~measurement~duration),
    subtitle = paste0("Stem collar geometry (flux term = ",
                      round(median_flux_term, 2),
                      "). MDF = precision / t \u00D7 flux term.")
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    plot.title = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 10, color = "grey40")
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
cat("  Wassmann 90-99%:        ",
    round(100 * mean(hf_out$CH4_below_MDF_wass90, na.rm = TRUE), 1), "-",
    round(100 * mean(hf_out$CH4_below_MDF_wass99, na.rm = TRUE), 1), "%\n")
cat("  Christiansen 90-99%:    ",
    round(100 * mean(hf_out$CH4_below_MDF_chr90, na.rm = TRUE), 1), "-",
    round(100 * mean(hf_out$CH4_below_MDF_chr99, na.rm = TRUE), 1), "%\n")

cat("\nStem (n=", nrow(df), "):\n")
for (inst in c("Analyzer A", "Analyzer B")) {
  mask <- if (inst == "Analyzer A") df$year < 2025 else df$year == 2025
  n_inst <- sum(mask)
  cat("  ", inst, " (n=", n_inst, "):\n", sep = "")
  cat("    Manufacturer:         ", round(100 * mean(df$CH4_below_MDF_manuf[mask], na.rm = TRUE), 1), "%\n")
  cat("    Wassmann 90-99%:      ",
      round(100 * mean(df$CH4_below_MDF_wass90[mask], na.rm = TRUE), 1), "-",
      round(100 * mean(df$CH4_below_MDF_wass99[mask], na.rm = TRUE), 1), "%\n")
  cat("    Christiansen 90-99%:  ",
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
cat("MDF at t=300s for manufacturer Analyzer A:",
    round(PREC_CH4_LGR / 300 * median_flux_term, 4), "nmol/m2/s\n")
cat("MDF at t=300s for empirical Analyzer A:",
    round(median(df$allan_sd_CH4[df$year < 2025], na.rm = TRUE) / 300 * median_flux_term, 4),
    "nmol/m2/s\n")
cat("MDF at t=300s for empirical Analyzer B:",
    round(median(df$allan_sd_CH4[df$year == 2025], na.rm = TRUE) / 300 * median_flux_term, 6),
    "nmol/m2/s\n")

message("\n=== MDF figures complete ===")
