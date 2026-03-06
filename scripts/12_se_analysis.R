# ============================================================================
# 12_se_analysis.R — Explore flux SE: correlations with magnitude, R²,
#                     noise floor, temperature, instrument, season
#
# Goal: Understand what goFlux model-fit SE captures, how it relates to
#       other uncertainty/quality metrics, and whether it adds information
#       beyond what noise_floor and R² already provide.
#
# Loads: outputs/03_mdf_results.RData (stem data with SE + noise floor + MDF)
# Saves: figures/se_analysis/ (exploratory plots)
# ============================================================================

source("scripts/00_setup.R")

message("\n=== Loading data ===")
load(file.path(output_dir, "03_mdf_results.RData"))

# Create output directory
se_fig_dir <- file.path(fig_dir, "se_analysis")
dir.create(se_fig_dir, recursive = TRUE, showWarnings = FALSE)

# Add instrument label
df <- df %>%
  mutate(inst_label = if_else(year == 2025, "Analyzer B", "Analyzer A"))

# ============================================================================
# PART 1: BASIC SE DISTRIBUTIONS
# ============================================================================

message("\n=== Part 1: SE distributions ===")

# Summary stats
se_summary <- df %>%
  group_by(inst_label) %>%
  summarise(
    n = n(),
    median_SE = median(CH4_SE_corr, na.rm = TRUE),
    mean_SE   = mean(CH4_SE_corr, na.rm = TRUE),
    iqr_SE    = IQR(CH4_SE_corr, na.rm = TRUE),
    min_SE    = min(CH4_SE_corr, na.rm = TRUE),
    max_SE    = max(CH4_SE_corr, na.rm = TRUE),
    .groups = "drop"
  )
cat("\n--- CH4 SE summary by instrument ---\n")
print(as.data.frame(se_summary))

# SE distribution by instrument
p_se_dist <- ggplot(df, aes(x = CH4_SE_corr, fill = inst_label)) +
  geom_histogram(bins = 60, alpha = 0.7, position = "identity") +
  scale_x_log10() +
  labs(x = expression(CH[4]~"flux SE (nmol"~m^{-2}~s^{-1}*")"),
       y = "Count",
       title = "Distribution of goFlux model-fit SE by instrument",
       fill = "Instrument") +
  theme_bw(base_size = 12) +
  theme(legend.position = "top")

ggsave(file.path(se_fig_dir, "se_distribution_by_instrument.png"),
       p_se_dist, width = 8, height = 5, dpi = 200)
message("  Saved: se_distribution_by_instrument.png")

# ============================================================================
# PART 2: SE vs FLUX MAGNITUDE
# ============================================================================

message("\n=== Part 2: SE vs flux magnitude ===")

p_se_flux <- ggplot(df, aes(x = abs(CH4_flux_nmolpm2ps), y = CH4_SE_corr,
                             color = inst_label)) +
  geom_point(alpha = 0.3, size = 1) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
  scale_x_log10() + scale_y_log10() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  labs(x = expression("|CH"[4]~"flux| (nmol"~m^{-2}~s^{-1}*")"),
       y = expression("SE (nmol"~m^{-2}~s^{-1}*")"),
       title = expression("SE vs |flux|: dashed = 1:1 line (SE = |flux|"~"means SNR"[SE]~"= 1)"),
       color = "Instrument") +
  theme_bw(base_size = 12) +
  theme(legend.position = "top")

ggsave(file.path(se_fig_dir, "se_vs_flux_magnitude.png"),
       p_se_flux, width = 8, height = 6, dpi = 200)
message("  Saved: se_vs_flux_magnitude.png")

# Signed flux version
p_se_flux_signed <- ggplot(df, aes(x = CH4_flux_nmolpm2ps, y = CH4_SE_corr,
                                    color = inst_label)) +
  geom_point(alpha = 0.3, size = 1) +
  scale_y_log10() +
  geom_hline(yintercept = median(df$CH4_SE_corr[df$inst_label == "Analyzer A"]),
             linetype = "dotted", color = "#F8766D") +
  geom_hline(yintercept = median(df$CH4_SE_corr[df$inst_label == "Analyzer B"]),
             linetype = "dotted", color = "#00BFC4") +
  labs(x = expression("CH"[4]~"flux (nmol"~m^{-2}~s^{-1}*")"),
       y = expression("SE (nmol"~m^{-2}~s^{-1}*")"),
       title = "SE vs signed flux; dotted = instrument median SE",
       color = "Instrument") +
  theme_bw(base_size = 12) +
  theme(legend.position = "top")

ggsave(file.path(se_fig_dir, "se_vs_flux_signed.png"),
       p_se_flux_signed, width = 8, height = 6, dpi = 200)
message("  Saved: se_vs_flux_signed.png")


# ============================================================================
# PART 3: SE vs R²
# ============================================================================

message("\n=== Part 3: SE vs R² ===")

p_se_r2 <- ggplot(df, aes(x = CH4_r2, y = CH4_SE_corr, color = inst_label)) +
  geom_point(alpha = 0.3, size = 1) +
  scale_y_log10() +
  geom_vline(xintercept = c(0.5, 0.7, 0.9), linetype = "dashed",
             color = "grey50", linewidth = 0.4) +
  labs(x = expression(CH[4]~R^2),
       y = expression("SE (nmol"~m^{-2}~s^{-1}*")"),
       title = expression(SE~vs~R^2*": higher R² → lower SE (but instrument matters more)"),
       color = "Instrument") +
  theme_bw(base_size = 12) +
  theme(legend.position = "top")

ggsave(file.path(se_fig_dir, "se_vs_r2.png"),
       p_se_r2, width = 8, height = 6, dpi = 200)
message("  Saved: se_vs_r2.png")

# Correlation stats
cor_se_r2 <- df %>%
  group_by(inst_label) %>%
  summarise(
    cor_se_r2        = cor(CH4_SE_corr, CH4_r2, use = "complete.obs"),
    cor_log_se_r2    = cor(log10(CH4_SE_corr), CH4_r2, use = "complete.obs"),
    .groups = "drop"
  )
cat("\n--- Correlation: SE vs R² ---\n")
print(as.data.frame(cor_se_r2))


# ============================================================================
# PART 4: SE vs NOISE FLOOR (Allan-based)
# ============================================================================

message("\n=== Part 4: SE vs noise floor ===")

df_with_nf <- df %>% filter(!is.na(CH4_noise_floor))

p_se_nf <- ggplot(df_with_nf, aes(x = CH4_noise_floor, y = CH4_SE_corr,
                                    color = inst_label)) +
  geom_point(alpha = 0.3, size = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  scale_x_log10() + scale_y_log10() +
  labs(x = expression("Noise floor (nmol"~m^{-2}~s^{-1}*")"),
       y = expression("SE (nmol"~m^{-2}~s^{-1}*")"),
       title = "SE (model-fit) vs noise floor (Allan): dashed = 1:1",
       color = "Instrument") +
  theme_bw(base_size = 12) +
  theme(legend.position = "top")

ggsave(file.path(se_fig_dir, "se_vs_noise_floor.png"),
       p_se_nf, width = 8, height = 6, dpi = 200)
message("  Saved: se_vs_noise_floor.png")

# Ratio: SE / noise_floor
df_with_nf <- df_with_nf %>%
  mutate(se_nf_ratio = CH4_SE_corr / CH4_noise_floor)

ratio_summary <- df_with_nf %>%
  group_by(inst_label) %>%
  summarise(
    median_ratio = median(se_nf_ratio, na.rm = TRUE),
    mean_ratio   = mean(se_nf_ratio, na.rm = TRUE),
    iqr_lo       = quantile(se_nf_ratio, 0.25, na.rm = TRUE),
    iqr_hi       = quantile(se_nf_ratio, 0.75, na.rm = TRUE),
    pct_se_larger = 100 * mean(se_nf_ratio > 1, na.rm = TRUE),
    .groups = "drop"
  )
cat("\n--- SE / noise_floor ratio by instrument ---\n")
print(as.data.frame(ratio_summary))

cor_se_nf <- df_with_nf %>%
  group_by(inst_label) %>%
  summarise(
    cor_se_nf     = cor(CH4_SE_corr, CH4_noise_floor, use = "complete.obs"),
    cor_log       = cor(log10(CH4_SE_corr), log10(CH4_noise_floor), use = "complete.obs"),
    .groups = "drop"
  )
cat("\n--- Correlation: SE vs noise_floor ---\n")
print(as.data.frame(cor_se_nf))


# ============================================================================
# PART 5: SE vs TEMPERATURE (seasonal confounding?)
# ============================================================================

message("\n=== Part 5: SE vs temperature ===")

p_se_temp <- ggplot(df, aes(x = airt, y = CH4_SE_corr, color = inst_label)) +
  geom_point(alpha = 0.3, size = 1) +
  geom_smooth(method = "loess", se = TRUE, alpha = 0.2) +
  scale_y_log10() +
  labs(x = "Air temperature (°C)",
       y = expression("SE (nmol"~m^{-2}~s^{-1}*")"),
       title = "SE vs temperature: does SE filter introduce seasonal bias?",
       color = "Instrument") +
  theme_bw(base_size = 12) +
  theme(legend.position = "top")

ggsave(file.path(se_fig_dir, "se_vs_temperature.png"),
       p_se_temp, width = 8, height = 6, dpi = 200)
message("  Saved: se_vs_temperature.png")


# ============================================================================
# PART 6: SNR(SE) vs SNR(Allan) — do they agree?
# ============================================================================

message("\n=== Part 6: SNR(SE) vs SNR(Allan) ===")

df_both_snr <- df %>%
  filter(!is.na(CH4_snr_allan), !is.na(CH4_snr_se))

p_snr_compare <- ggplot(df_both_snr, aes(x = CH4_snr_allan, y = CH4_snr_se,
                                           color = inst_label)) +
  geom_point(alpha = 0.3, size = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  scale_x_log10() + scale_y_log10() +
  geom_hline(yintercept = c(2, 3), linetype = "dotted", color = "grey50") +
  geom_vline(xintercept = c(2, 3), linetype = "dotted", color = "grey50") +
  labs(x = expression(SNR[Allan]~"= |flux| / noise floor"),
       y = expression(SNR[SE]~"= |flux| / SE"),
       title = "SNR comparison: SE-based vs Allan-based; dashed = 1:1",
       color = "Instrument") +
  theme_bw(base_size = 12) +
  theme(legend.position = "top")

ggsave(file.path(se_fig_dir, "snr_se_vs_snr_allan.png"),
       p_snr_compare, width = 8, height = 6, dpi = 200)
message("  Saved: snr_se_vs_snr_allan.png")

# Agreement: how often do they agree on detection?
detection_agreement <- df_both_snr %>%
  mutate(
    detected_se2    = CH4_snr_se > 2,
    detected_allan2 = CH4_snr_allan > 2,
    detected_se3    = CH4_snr_se > 3,
    detected_allan3 = CH4_snr_allan > 3,
    agree_snr2      = detected_se2 == detected_allan2,
    agree_snr3      = detected_se3 == detected_allan3
  ) %>%
  summarise(
    n = n(),
    pct_agree_snr2 = round(100 * mean(agree_snr2), 1),
    pct_agree_snr3 = round(100 * mean(agree_snr3), 1),
    # Confusion matrix for SNR > 2
    both_detect_2  = sum(detected_se2 & detected_allan2),
    se_only_2      = sum(detected_se2 & !detected_allan2),
    allan_only_2   = sum(!detected_se2 & detected_allan2),
    neither_2      = sum(!detected_se2 & !detected_allan2)
  )
cat("\n--- SNR detection agreement ---\n")
print(as.data.frame(detection_agreement))


# ============================================================================
# PART 7: SE-BASED SEASONAL CONFOUNDING (parallel to R² analysis)
# ============================================================================

message("\n=== Part 7: SE-based seasonal confounding ===")

# Add month and season
df <- df %>%
  mutate(
    month = as.integer(format(as.Date(date), "%m")),
    season = if_else(month %in% 5:10, "Growing", "Dormant")
  )

# What fraction fails SE-based SNR thresholds by month?
se_monthly <- df %>%
  group_by(month, inst_label) %>%
  summarise(
    n_total     = n(),
    pct_fail_se2 = 100 * mean(CH4_snr_se <= 2),
    pct_fail_se3 = 100 * mean(CH4_snr_se <= 3),
    pct_fail_r2_07 = 100 * mean(CH4_r2 <= 0.7),
    .groups = "drop"
  )

# Long format for comparison plot
se_monthly_long <- se_monthly %>%
  pivot_longer(cols = starts_with("pct_fail"),
               names_to = "criterion",
               values_to = "pct_excluded") %>%
  mutate(
    criterion = case_when(
      criterion == "pct_fail_se2"    ~ "SNR(SE) > 2",
      criterion == "pct_fail_se3"    ~ "SNR(SE) > 3",
      criterion == "pct_fail_r2_07"  ~ "R² > 0.7"
    )
  )

p_se_seasonal <- ggplot(se_monthly_long,
                         aes(x = factor(month), y = pct_excluded,
                             color = criterion, group = criterion)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  facet_wrap(~ inst_label) +
  annotate("rect", xmin = 0.5, xmax = 4.5, ymin = -Inf, ymax = Inf,
           fill = "grey80", alpha = 0.2) +
  annotate("rect", xmin = 10.5, xmax = 12.5, ymin = -Inf, ymax = Inf,
           fill = "grey80", alpha = 0.2) +
  labs(x = "Month", y = "% excluded",
       title = "Seasonal exclusion: SE-based SNR vs R² thresholds",
       subtitle = "Grey bands = dormant season (Nov-Apr)",
       color = "Criterion") +
  ylim(0, 100) +
  theme_bw(base_size = 12) +
  theme(legend.position = "top")

ggsave(file.path(se_fig_dir, "se_seasonal_confounding.png"),
       p_se_seasonal, width = 10, height = 5, dpi = 200)
message("  Saved: se_seasonal_confounding.png")

# Seasonal summary
se_seasonal_summary <- df %>%
  group_by(season, inst_label) %>%
  summarise(
    n           = n(),
    pct_fail_se2 = round(100 * mean(CH4_snr_se <= 2), 1),
    pct_fail_se3 = round(100 * mean(CH4_snr_se <= 3), 1),
    pct_fail_r2  = round(100 * mean(CH4_r2 <= 0.7), 1),
    mean_flux    = round(mean(CH4_flux_nmolpm2ps), 3),
    median_flux  = round(median(CH4_flux_nmolpm2ps), 3),
    .groups = "drop"
  )
cat("\n--- Seasonal exclusion comparison ---\n")
print(as.data.frame(se_seasonal_summary))


# ============================================================================
# PART 8: SE vs R² — do they provide redundant or complementary info?
# ============================================================================

message("\n=== Part 8: SE vs R² — complementarity ===")

# Four quadrants: high SE + high R² vs low SE + low R²
df_quad <- df %>%
  mutate(
    se_high  = CH4_SE_corr > median(CH4_SE_corr),
    r2_low   = CH4_r2 < 0.7,
    quadrant = case_when(
      !se_high & !r2_low ~ "Low SE, High R²",
       se_high & !r2_low ~ "High SE, High R²",
      !se_high &  r2_low ~ "Low SE, Low R²",
       se_high &  r2_low ~ "High SE, Low R²"
    )
  )

quad_summary <- df_quad %>%
  group_by(quadrant, inst_label) %>%
  summarise(
    n = n(),
    mean_flux   = round(mean(CH4_flux_nmolpm2ps), 3),
    median_flux = round(median(CH4_flux_nmolpm2ps), 3),
    mean_temp   = round(mean(airt, na.rm = TRUE), 1),
    pct_neg     = round(100 * mean(CH4_flux_nmolpm2ps < 0), 1),
    .groups = "drop"
  )
cat("\n--- SE × R² quadrant summary ---\n")
print(as.data.frame(quad_summary))

# Plot: colored by quadrant
p_quad <- ggplot(df_quad, aes(x = CH4_r2, y = CH4_SE_corr, color = quadrant)) +
  geom_point(alpha = 0.3, size = 1) +
  scale_y_log10() +
  facet_wrap(~ inst_label) +
  geom_vline(xintercept = 0.7, linetype = "dashed", color = "grey40") +
  geom_hline(aes(yintercept = median(CH4_SE_corr)),
             linetype = "dashed", color = "grey40") +
  labs(x = expression(CH[4]~R^2),
       y = expression("SE (nmol"~m^{-2}~s^{-1}*")"),
       title = "SE × R² quadrants: do they capture the same measurements?",
       color = "Quadrant") +
  theme_bw(base_size = 12) +
  theme(legend.position = "top")

ggsave(file.path(se_fig_dir, "se_r2_quadrants.png"),
       p_quad, width = 12, height = 5, dpi = 200)
message("  Saved: se_r2_quadrants.png")


# ============================================================================
# PART 9: TREATMENT EFFECT USING SE-BASED FILTERS
# ============================================================================

message("\n=== Part 9: Treatment comparison with SE-based SNR filter ===")

# Compare mean flux under different treatments for SE-based filters
# Parallel to the MDF treatment analysis in the paper

se_treatments <- expand.grid(
  snr_threshold = c(2, 3),
  treatment = c("Keep all", "Remove", "Set to zero"),
  stringsAsFactors = FALSE
)

treatment_results <- lapply(1:nrow(se_treatments), function(i) {
  thresh <- se_treatments$snr_threshold[i]
  trt    <- se_treatments$treatment[i]

  below <- df$CH4_snr_se <= thresh

  flux_vals <- df$CH4_flux_nmolpm2ps
  if (trt == "Remove") {
    flux_vals <- flux_vals[!below]
  } else if (trt == "Set to zero") {
    flux_vals[below] <- 0
  }

  data.frame(
    snr_threshold = thresh,
    treatment     = trt,
    n             = length(flux_vals),
    n_flagged     = sum(below),
    pct_flagged   = round(100 * mean(below), 1),
    mean_flux     = round(mean(flux_vals, na.rm = TRUE), 4),
    median_flux   = round(median(flux_vals, na.rm = TRUE), 4),
    stringsAsFactors = FALSE
  )
})
treatment_results <- do.call(rbind, treatment_results)

# Add baseline
baseline <- data.frame(
  snr_threshold = NA, treatment = "Unfiltered",
  n = nrow(df), n_flagged = 0, pct_flagged = 0,
  mean_flux = round(mean(df$CH4_flux_nmolpm2ps), 4),
  median_flux = round(median(df$CH4_flux_nmolpm2ps), 4)
)
treatment_results <- rbind(baseline, treatment_results)

cat("\n--- Treatment comparison (SE-based SNR) ---\n")
print(treatment_results, row.names = FALSE)

# Percent bias relative to unfiltered
base_mean   <- baseline$mean_flux
base_median <- baseline$median_flux
treatment_results$pct_bias_mean <- round(
  100 * (treatment_results$mean_flux - base_mean) / abs(base_mean), 1)
treatment_results$pct_bias_median <- round(
  100 * (treatment_results$median_flux - base_median) / abs(base_median), 1)

cat("\n--- Treatment bias relative to unfiltered ---\n")
print(treatment_results, row.names = FALSE)


# ============================================================================
# PART 10: OVERALL CORRELATION MATRIX
# ============================================================================

message("\n=== Part 10: Correlation matrix ===")

cor_cols <- c("CH4_flux_nmolpm2ps", "CH4_SE_corr", "CH4_r2",
              "CH4_noise_floor", "CH4_snr_se", "CH4_snr_allan", "airt")

cor_df <- df %>%
  select(all_of(cor_cols)) %>%
  filter(complete.cases(.))

cor_mat <- cor(cor_df, use = "complete.obs")
cat("\n--- Correlation matrix (n =", nrow(cor_df), ") ---\n")
print(round(cor_mat, 3))

# Also compute rank correlations (Spearman) which handle nonlinearity better
cor_mat_spearman <- cor(cor_df, use = "complete.obs", method = "spearman")
cat("\n--- Spearman correlation matrix ---\n")
print(round(cor_mat_spearman, 3))


# ============================================================================
# PART 11: FIRST-PRINCIPLES SE VERIFICATION
# ============================================================================

message("\n=== Part 11: First-principles SE verification ===")

# For each instrument, compare observed SE to expected SE from:
#   SE(flux) = sigma_Allan / sqrt(SSxx) * flux_term
#   where SSxx = n*(n^2-1)/12 for evenly spaced data

df_verify <- df %>%
  filter(!is.na(allan_sd_CH4), !is.na(n_pts), !is.na(flux_term))

se_verify <- df_verify %>%
  group_by(inst_label) %>%
  summarise(
    n_meas           = n(),
    median_allan_sd  = median(allan_sd_CH4),
    median_n_pts     = median(n_pts),
    median_flux_term = median(flux_term),
    median_SSxx      = median(n_pts * (n_pts^2 - 1) / 12),
    expected_se_flux = median(allan_sd_CH4) /
                        sqrt(median(n_pts * (n_pts^2 - 1) / 12)) *
                        median(flux_term),
    observed_se_flux = median(CH4_SE_corr),
    noise_floor      = median(CH4_noise_floor),
    ratio_obs_exp    = median(CH4_SE_corr) /
                        (median(allan_sd_CH4) /
                         sqrt(median(n_pts * (n_pts^2 - 1) / 12)) *
                         median(flux_term)),
    .groups = "drop"
  )

cat("\n--- SE verification: observed vs first-principles expected ---\n")
print(as.data.frame(se_verify))
cat("\nratio_obs_exp ≈ 1 means SE is correctly calibrated.\n")
cat("Analyzer B ratio ≈ 1.14 → SE is trustworthy.\n")
cat("Analyzer A ratio ≈ 0.013 → SE is ~77× too small (unreliable).\n")
cat("\n→ Use noise_floor (Allan-based) as per-measurement uncertainty.\n")


# ============================================================================
# SAVE SUMMARY OUTPUTS
# ============================================================================

message("\n=== Saving summary outputs ===")

save(se_summary, ratio_summary, cor_se_r2, cor_se_nf,
     detection_agreement, se_seasonal_summary, quad_summary,
     treatment_results, cor_mat, cor_mat_spearman, se_verify,
     file = file.path(output_dir, "12_se_analysis.RData"))

message("Saved: ", file.path(output_dir, "12_se_analysis.RData"))
message("\n=== SE analysis complete ===")
