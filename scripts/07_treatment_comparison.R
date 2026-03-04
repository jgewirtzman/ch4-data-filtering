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
# Plots:
#   1. Treatment ridges (overlaid density for 3 treatments per threshold)
#   2. Bias dot-line plot (mean/median shift under each treatment)
#   3. Treatment by instrument (faceted density showing instrument effects)
#   4. % change bars (percent change in mean/median relative to unfiltered)
#
# Loads: outputs/03_mdf_results.RData
# Saves: figures to fig_dir (PNG + PDF)
# ============================================================================

source("scripts/00_setup.R")

message("\n=== Loading MDF results from 03 ===")
load(file.path(output_dir, "03_mdf_results.RData"))

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ── Define asinh transform for ggplot scales ─────────────────────────────────
asinh_trans <- trans_new(
  name = "asinh",
  transform = asinh,
  inverse = sinh,
  breaks = function(x) {
    rng <- sinh(x)
    pretty(rng, n = 8)
  }
)

# ── Ensure inst_label exists ─────────────────────────────────────────────────
if (!"inst_label" %in% names(df)) {
  df <- df %>%
    mutate(
      inst_label = case_when(
        year == 2025 ~ "LI-7810 (2025)",
        year < 2025  ~ "LGR/UGGA (2023-24)"
      )
    )
}

# ============================================================
# PART 12: BELOW-MDF TREATMENT COMPARISON
# ============================================================

message("\n=== Below-MDF treatment comparison ===")

# Focus on MDF-based filters only (these have clear thresholds)
mdf_filters <- list(
  "Manufacturer MDF" = "CH4_below_MDF_manuf",
  "Wassmann 90%"     = "CH4_below_MDF_wass90",
  "Wassmann 95%"     = "CH4_below_MDF_wass95",
  "Wassmann 99%"     = "CH4_below_MDF_wass99",
  "Christiansen 90%" = "CH4_below_MDF_chr90",
  "Christiansen 95%" = "CH4_below_MDF_chr95",
  "Christiansen 99%" = "CH4_below_MDF_chr99"
)

# --- Build long data frame with three treatments per filter ---
treat_rows <- list()
k <- 0

for (filt_name in names(mdf_filters)) {
  flag_col <- mdf_filters[[filt_name]]
  below <- df[[flag_col]]
  below[is.na(below)] <- FALSE  # if no Allan data, keep measurement

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
filter_order_mdf <- c("Manufacturer MDF",
                       "Wassmann 90%", "Wassmann 95%", "Wassmann 99%",
                       "Christiansen 90%", "Christiansen 95%", "Christiansen 99%")
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

# --- PLOT 1: Overlaid densities, faceted by filter ---
# Shows how each treatment changes the distribution shape

p_treat_ridges <- ggplot(treat_df,
                          aes(x = CH4_flux, y = filter,
                              fill = treatment, color = treatment)) +
  geom_density_ridges(
    alpha = 0.35, scale = 0.85,
    bandwidth = 0.25,
    from = asinh(-3), to = asinh(20),
    rel_min_height = 0.005
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red", linewidth = 0.4) +
  scale_x_continuous(
    trans = asinh_trans,
    breaks = c(-2, -1, 0, 0.5, 1, 2, 5, 10),
    labels = c("-2", "-1", "0", "0.5", "1", "2", "5", "10")
  ) +
  scale_fill_manual(
    values = c("Keep all (unfiltered)" = "grey70",
               "Set below MDF to zero" = "#2166AC",
               "Remove below MDF"      = "#B2182B"),
    name = "Below-MDF treatment"
  ) +
  scale_color_manual(
    values = c("Keep all (unfiltered)" = "grey40",
               "Set below MDF to zero" = "#2166AC",
               "Remove below MDF"      = "#B2182B"),
    name = "Below-MDF treatment"
  ) +
  labs(
    x = expression(CH[4]~flux~(nmol~m^{-2}~s^{-1})),
    y = NULL,
    title = expression(Effect~of~below-MDF~treatment~on~CH[4]~flux~distribution),
    subtitle = paste0("n = ", n_total, " measurements | ",
                      "Remove: excludes below-MDF | ",
                      "Set to zero: replaces with 0 | ",
                      "Unfiltered: keeps original values")
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.y  = element_text(size = 9),
    legend.position = "bottom",
    plot.title    = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 8.5, color = "grey40")
  ) +
  coord_cartesian(xlim = c(asinh(-3), asinh(15)))

ggsave(file.path(fig_dir, "ch4_mdf_treatment_ridges.png"),
       p_treat_ridges, width = 11, height = 8, dpi = 300)
ggsave(file.path(fig_dir, "ch4_mdf_treatment_ridges.pdf"),
       p_treat_ridges, width = 11, height = 8)
message("Saved MDF treatment ridges to: ", fig_dir)

# --- PLOT 2: Mean / median shift by treatment (dot-and-line) ---
# Directly shows the positive bias from removal

treat_summary_long <- treat_stats %>%
  select(filter, treatment, mean_flux, median_flux) %>%
  pivot_longer(cols = c(mean_flux, median_flux),
               names_to = "statistic", values_to = "value") %>%
  mutate(statistic = recode(statistic,
                             "mean_flux" = "Mean",
                             "median_flux" = "Median"))

p_bias <- ggplot(treat_summary_long,
                  aes(x = filter, y = value, color = treatment, group = treatment)) +
  geom_line(linewidth = 0.7, alpha = 0.8) +
  geom_point(size = 2.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  facet_wrap(~ statistic, ncol = 1, scales = "free_y") +
  scale_color_manual(
    values = c("Keep all (unfiltered)" = "grey50",
               "Set below MDF to zero" = "#2166AC",
               "Remove below MDF"      = "#B2182B"),
    name = "Below-MDF treatment"
  ) +
  scale_y_continuous(
    trans = asinh_trans,
    breaks = c(0, 0.05, 0.1, 0.2, 0.5, 1, 2, 5),
    labels = c("0", "0.05", "0.1", "0.2", "0.5", "1", "2", "5")
  ) +
  labs(
    x = NULL,
    y = expression(CH[4]~flux~(nmol~m^{-2}~s^{-1})),
    title = expression(Positive~bias~from~removing~below-MDF~measurements),
    subtitle = "Removing below-MDF shifts mean/median upward; setting to zero preserves sample size with minimal bias"
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x  = element_text(angle = 30, hjust = 1, size = 8),
    legend.position = "bottom",
    plot.title    = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 9, color = "grey40"),
    strip.text    = element_text(size = 11, face = "bold")
  )

ggsave(file.path(fig_dir, "ch4_mdf_treatment_bias.png"),
       p_bias, width = 10, height = 7, dpi = 300)
ggsave(file.path(fig_dir, "ch4_mdf_treatment_bias.pdf"),
       p_bias, width = 10, height = 7)
message("Saved MDF treatment bias plot to: ", fig_dir)

# --- PLOT 3: Same as PLOT 1 but split by instrument ---
# Shows the differential impact: removal biases LGR far more than LI-7810

p_treat_inst <- ggplot(treat_df,
                        aes(x = CH4_flux, fill = treatment, color = treatment)) +
  geom_density(alpha = 0.3, linewidth = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red", linewidth = 0.3) +
  facet_grid(filter ~ instrument, scales = "free_y") +
  scale_x_continuous(
    trans = asinh_trans,
    breaks = c(-2, -1, 0, 0.5, 1, 2, 5, 10),
    labels = c("-2", "-1", "0", "0.5", "1", "2", "5", "10")
  ) +
  scale_fill_manual(
    values = c("Keep all (unfiltered)" = "grey70",
               "Set below MDF to zero" = "#2166AC",
               "Remove below MDF"      = "#B2182B"),
    name = "Treatment"
  ) +
  scale_color_manual(
    values = c("Keep all (unfiltered)" = "grey40",
               "Set below MDF to zero" = "#2166AC",
               "Remove below MDF"      = "#B2182B"),
    name = "Treatment"
  ) +
  labs(
    x = expression(CH[4]~flux~(nmol~m^{-2}~s^{-1})),
    y = "Density",
    title = expression(Below-MDF~treatment~effect~by~instrument),
    subtitle = "Removal causes rightward shift primarily for LGR; LI-7810 is minimally affected"
  ) +
  theme_classic(base_size = 9) +
  theme(
    legend.position = "bottom",
    plot.title    = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 9, color = "grey40"),
    strip.text.y  = element_text(size = 7, angle = 0),
    strip.text.x  = element_text(size = 9, face = "bold")
  ) +
  coord_cartesian(xlim = c(asinh(-3), asinh(15)))

ggsave(file.path(fig_dir, "ch4_mdf_treatment_by_instrument.png"),
       p_treat_inst, width = 12, height = 14, dpi = 300)
ggsave(file.path(fig_dir, "ch4_mdf_treatment_by_instrument.pdf"),
       p_treat_inst, width = 12, height = 14)
message("Saved MDF treatment x instrument plot to: ", fig_dir)

# --- PLOT 4: Percent change in mean/median from removal vs set-to-zero ---
# relative to the unfiltered baseline

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

pct_long <- pct_change %>%
  select(filter, treatment, pct_change_mean, pct_change_median) %>%
  pivot_longer(cols = c(pct_change_mean, pct_change_median),
               names_to = "statistic", values_to = "pct_change") %>%
  mutate(statistic = recode(statistic,
                             "pct_change_mean"   = "Mean",
                             "pct_change_median" = "Median"))

p_pct <- ggplot(pct_long,
                 aes(x = filter, y = pct_change, fill = treatment)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.85) +
  geom_hline(yintercept = 0, color = "grey30") +
  facet_wrap(~ statistic, ncol = 1) +
  scale_fill_manual(
    values = c("Set below MDF to zero" = "#2166AC",
               "Remove below MDF"      = "#B2182B"),
    name = "Treatment"
  ) +
  labs(
    x = NULL,
    y = "% change relative to unfiltered",
    title = expression(Bias~introduced~by~below-MDF~treatment),
    subtitle = "Removing below-MDF inflates mean/median; setting to zero introduces much less bias"
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, size = 8),
    legend.position = "bottom",
    plot.title    = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 9, color = "grey40"),
    strip.text    = element_text(size = 11, face = "bold")
  )

ggsave(file.path(fig_dir, "ch4_mdf_treatment_pct_change.png"),
       p_pct, width = 10, height = 7, dpi = 300)
ggsave(file.path(fig_dir, "ch4_mdf_treatment_pct_change.pdf"),
       p_pct, width = 10, height = 7)
message("Saved MDF treatment % change to: ", fig_dir)

message("\n=== Treatment comparison complete ===")
