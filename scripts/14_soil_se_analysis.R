# ============================================================================
# 14_soil_se_analysis.R — SE analysis on soil CH₄ fluxes from goFlux
#
# Uses the FINAL corrected soil dataset from tree-methanogens, which has
# full goFlux output including per-measurement SE for both LM and HM models.
# All measurements are from the UGGA (LGR3, same instrument family as GLA131)
# with manufacturer CH4 precision = 0.9 ppb.
#
# This provides a clean test case for SE analysis because:
#   1. All data processed through one consistent goFlux pipeline
#   2. Single instrument (no instrument-mixing issues)
#   3. Soil CH4 fluxes are large and well-characterized
#   4. goFlux SE values are available and unmodified
#
# Analyses:
#   Part 1: SE distributions and basic summary
#   Part 2: SE vs flux magnitude
#   Part 3: SE vs R², MDF
#   Part 4: First-principles SE verification (Allan-based)
#   Part 5: SNR analysis (SE-based vs MDF-based)
#   Part 6: Scaling uncertainty (SE of mean, bootstrap, random effects)
#   Part 7: Variance decomposition
#
# Loads: tree-methanogens CH4_best_flux_lgr_results_soil.csv (raw goFlux)
#        + FINAL_soil_dataset_CH4_corrected.csv (for metadata)
# Saves: figures/soil_se_analysis/
# ============================================================================

source("scripts/00_setup.R")

# metafor for random-effects meta-analysis
if (!requireNamespace("metafor", quietly = TRUE)) {
  install.packages("metafor", repos = "https://cloud.r-project.org")
}
library(metafor)

message("\n=== Loading soil goFlux data ===")

soil_data_dir <- file.path(data_dir, "soil")

# Raw goFlux best-flux output (has correct numeric best.flux)
# Source: tree-methanogens/data/processed/flux/
flux_raw <- read.csv(file.path(soil_data_dir,
                                "CH4_best_flux_lgr_results_soil.csv"),
                      stringsAsFactors = FALSE)
message("goFlux best-flux rows: ", nrow(flux_raw))

# FINAL merged dataset (for metadata: Date, Site, Plot)
# Note: CH4_best.flux column in FINAL file has a merge bug (contains model
# name instead of flux value), so we use the raw file for flux/SE values.
final_raw <- read.csv(file.path(soil_data_dir,
                                 "FINAL_soil_dataset_CH4_corrected.csv"),
                       stringsAsFactors = FALSE)
message("FINAL metadata rows: ", nrow(final_raw))

# Join metadata to flux results
meta_cols <- final_raw %>%
  select(UniqueID, Date, Site, `Plot.Tag`, `Plot.letter`) %>%
  distinct()

soil <- flux_raw %>%
  filter(!is.na(best.flux)) %>%
  left_join(meta_cols, by = "UniqueID") %>%
  mutate(
    # Best SE from the selected model
    CH4_best_SE = ifelse(model == "HM", HM.SE, LM.SE),
    CH4_best_r2 = ifelse(model == "HM", HM.r2, LM.r2),
    CH4_abs_flux = abs(best.flux),
    CH4_SNR_se  = CH4_abs_flux / CH4_best_SE,
    CH4_SNR_mdf = CH4_abs_flux / MDF,
    date = as.Date(Date),
    month = as.integer(format(date, "%m")),
    season = ifelse(month %in% 5:10, "Growing", "Dormant")
  )

message("After filtering: ", nrow(soil), " measurements")
message("  Model selection: LM = ", sum(soil$model == "LM"),
        ", HM = ", sum(soil$model == "HM"))

soil_fig_dir <- file.path(fig_dir, "soil_se_analysis")
dir.create(soil_fig_dir, recursive = TRUE, showWarnings = FALSE)


# ============================================================================
# PART 1: SE DISTRIBUTIONS
# ============================================================================

message("\n=== Part 1: SE distributions ===")

cat("\n--- CH4 best SE summary ---\n")
cat(sprintf("  n     = %d\n", nrow(soil)))
cat(sprintf("  mean  = %.6f\n", mean(soil$CH4_best_SE, na.rm = TRUE)))
cat(sprintf("  median= %.6f\n", median(soil$CH4_best_SE, na.rm = TRUE)))
cat(sprintf("  SD    = %.6f\n", sd(soil$CH4_best_SE, na.rm = TRUE)))
cat(sprintf("  min   = %.6f\n", min(soil$CH4_best_SE, na.rm = TRUE)))
cat(sprintf("  max   = %.6f\n", max(soil$CH4_best_SE, na.rm = TRUE)))

cat("\n--- CH4 flux summary ---\n")
cat(sprintf("  mean  = %.4f\n", mean(soil$best.flux, na.rm = TRUE)))
cat(sprintf("  median= %.4f\n", median(soil$best.flux, na.rm = TRUE)))
cat(sprintf("  SD    = %.4f\n", sd(soil$best.flux, na.rm = TRUE)))
cat(sprintf("  range = [%.4f, %.4f]\n",
            min(soil$best.flux, na.rm = TRUE),
            max(soil$best.flux, na.rm = TRUE)))
cat(sprintf("  %% negative = %.1f%%\n",
            100 * mean(soil$best.flux < 0, na.rm = TRUE)))

cat("\n--- CH4 MDF summary ---\n")
cat(sprintf("  median MDF = %.6f nmol/m2/s\n",
            median(soil$MDF, na.rm = TRUE)))
cat(sprintf("  median SE  = %.6f nmol/m2/s\n",
            median(soil$CH4_best_SE, na.rm = TRUE)))
cat(sprintf("  SE / MDF ratio (median) = %.3f\n",
            median(soil$CH4_best_SE, na.rm = TRUE) /
              median(soil$MDF, na.rm = TRUE)))

p_se_hist <- ggplot(soil, aes(x = CH4_best_SE)) +
  geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
  geom_vline(xintercept = median(soil$CH4_best_SE, na.rm = TRUE),
             linetype = "dashed", color = "red") +
  labs(x = "CH4 best-model SE (nmol/m2/s)",
       y = "Count",
       title = "Soil CH4 goFlux SE distribution",
       subtitle = sprintf("n = %d, median = %.4f",
                            nrow(soil),
                            median(soil$CH4_best_SE, na.rm = TRUE))) +
  theme_bw(base_size = 12)

ggsave(file.path(soil_fig_dir, "se_distribution.png"),
       p_se_hist, width = 7, height = 5, dpi = 200)


# ============================================================================
# PART 2: SE vs FLUX MAGNITUDE
# ============================================================================

message("\n=== Part 2: SE vs flux magnitude ===")

# Does SE scale with |flux|? (It should if instrument noise is relative)
cor_se_flux <- cor(soil$CH4_best_SE, soil$CH4_abs_flux,
                    use = "complete.obs", method = "spearman")
cat(sprintf("\n  Spearman cor(SE, |flux|) = %.3f\n", cor_se_flux))

p_se_vs_flux <- ggplot(soil, aes(x = CH4_abs_flux, y = CH4_best_SE)) +
  geom_point(alpha = 0.4, size = 1.5) +
  geom_smooth(method = "lm", se = TRUE, color = "red", linewidth = 0.8) +
  scale_x_log10() + scale_y_log10() +
  labs(x = "|CH4 flux| (nmol/m2/s)",
       y = "CH4 SE (nmol/m2/s)",
       title = "Soil: SE vs flux magnitude",
       subtitle = sprintf("Spearman rho = %.3f (log-log)", cor_se_flux)) +
  theme_bw(base_size = 12)

ggsave(file.path(soil_fig_dir, "se_vs_flux.png"),
       p_se_vs_flux, width = 7, height = 5, dpi = 200)


# ============================================================================
# PART 3: SE vs R² and MDF
# ============================================================================

message("\n=== Part 3: SE vs R² and MDF ===")

cor_se_r2 <- cor(soil$CH4_best_SE, soil$CH4_best_r2,
                  use = "complete.obs", method = "spearman")
cor_se_mdf <- cor(soil$CH4_best_SE, soil$MDF,
                   use = "complete.obs", method = "spearman")
cat(sprintf("  Spearman cor(SE, R²)  = %.3f\n", cor_se_r2))
cat(sprintf("  Spearman cor(SE, MDF) = %.3f\n", cor_se_mdf))

p_se_vs_r2 <- ggplot(soil, aes(x = CH4_best_r2, y = CH4_best_SE)) +
  geom_point(alpha = 0.4, size = 1.5) +
  scale_y_log10() +
  labs(x = "CH4 R²", y = "CH4 SE (nmol/m2/s, log scale)",
       title = "Soil: SE vs R²",
       subtitle = sprintf("Spearman rho = %.3f", cor_se_r2)) +
  theme_bw(base_size = 12)

p_se_vs_mdf <- ggplot(soil, aes(x = MDF, y = CH4_best_SE)) +
  geom_point(alpha = 0.4, size = 1.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  scale_x_log10() + scale_y_log10() +
  labs(x = "CH4 MDF (nmol/m2/s)", y = "CH4 SE (nmol/m2/s)",
       title = "Soil: SE vs MDF (log-log)",
       subtitle = sprintf("Spearman rho = %.3f; red = 1:1 line", cor_se_mdf)) +
  theme_bw(base_size = 12)

ggsave(file.path(soil_fig_dir, "se_vs_r2.png"),
       p_se_vs_r2, width = 7, height = 5, dpi = 200)
ggsave(file.path(soil_fig_dir, "se_vs_mdf.png"),
       p_se_vs_mdf, width = 7, height = 5, dpi = 200)


# ============================================================================
# PART 4: FIRST-PRINCIPLES SE VERIFICATION
# ============================================================================

message("\n=== Part 4: First-principles SE verification ===")

# For a linear model fit to concentration vs time:
#   SE_slope = sigma_noise / sqrt(sum((t - t_mean)^2))
# For evenly-spaced t with n points:
#   sum((t - t_mean)^2) = n(n^2 - 1)/12  (for unit spacing)
# Then SE_flux = SE_slope × flux_term
#   where flux_term = V·P/(A·R·T) converts concentration slope to flux
#
# sigma_noise ≈ instrument precision (0.9 ppb for UGGA)
# or ≈ Allan deviation if available

# Use manufacturer precision as sigma_noise (0.9 ppb)
soil_verify <- soil %>%
  filter(!is.na(nb.obs), !is.na(flux.term), nb.obs > 2) %>%
  mutate(
    # Expected SE from first principles (LM only, using manufacturer precision)
    denom_fp = sqrt(nb.obs * (nb.obs^2 - 1) / 12),
    expected_SE_mfr = (prec / denom_fp) * abs(flux.term),
    # Ratio: observed vs expected
    ratio_obs_exp_mfr = LM.SE / expected_SE_mfr
  )

cat("\n--- First-principles SE check (manufacturer precision = 0.9 ppb) ---\n")
cat(sprintf("  n measurements with valid data: %d\n", nrow(soil_verify)))
cat(sprintf("  Median observed LM.SE:   %.6f\n",
            median(soil_verify$LM.SE, na.rm = TRUE)))
cat(sprintf("  Median expected SE:      %.6f\n",
            median(soil_verify$expected_SE_mfr, na.rm = TRUE)))
cat(sprintf("  Median ratio (obs/exp):  %.3f\n",
            median(soil_verify$ratio_obs_exp_mfr, na.rm = TRUE)))
cat(sprintf("  Mean ratio (obs/exp):    %.3f\n",
            mean(soil_verify$ratio_obs_exp_mfr, na.rm = TRUE)))
cat(sprintf("  SD of ratio:             %.3f\n",
            sd(soil_verify$ratio_obs_exp_mfr, na.rm = TRUE)))

cat("\n  Interpretation:\n")
cat("    ratio ~ 1: goFlux SE matches first principles (good!)\n")
cat("    ratio >> 1: goFlux SE larger than expected (conservative)\n")
cat("    ratio << 1: goFlux SE smaller than expected (like GLA131 stem)\n")

p_ratio <- ggplot(soil_verify, aes(x = ratio_obs_exp_mfr)) +
  geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red", linewidth = 0.8) +
  geom_vline(xintercept = median(soil_verify$ratio_obs_exp_mfr, na.rm = TRUE),
             linetype = "solid", color = "darkblue", linewidth = 0.8) +
  labs(x = "Ratio: observed LM.SE / expected SE (first principles)",
       y = "Count",
       title = "Soil: First-principles SE verification (UGGA)",
       subtitle = sprintf("Median ratio = %.2f (red = 1.0, blue = median)",
                            median(soil_verify$ratio_obs_exp_mfr, na.rm = TRUE))) +
  theme_bw(base_size = 12)

ggsave(file.path(soil_fig_dir, "se_first_principles_ratio.png"),
       p_ratio, width = 7, height = 5, dpi = 200)

# Also check: does the ratio vary with flux magnitude? (it shouldn't)
cor_ratio_flux <- cor(soil_verify$ratio_obs_exp_mfr,
                       soil_verify$CH4_abs_flux,
                       use = "complete.obs", method = "spearman")
cat(sprintf("\n  Spearman cor(ratio, |flux|) = %.3f\n", cor_ratio_flux))
cat("    (should be ~0 if SE is well-calibrated across flux range)\n")


# ============================================================================
# PART 5: SNR ANALYSIS — SE-based vs MDF-based
# ============================================================================

message("\n=== Part 5: SNR comparison ===")

cat("\n--- SNR distributions ---\n")
cat(sprintf("  SNR(SE)  median = %.2f\n", median(soil$CH4_SNR_se, na.rm = TRUE)))
cat(sprintf("  SNR(MDF) median = %.2f\n", median(soil$CH4_SNR_mdf, na.rm = TRUE)))
cat(sprintf("  Spearman cor(SNR_SE, SNR_MDF) = %.3f\n",
            cor(soil$CH4_SNR_se, soil$CH4_SNR_mdf,
                use = "complete.obs", method = "spearman")))

# Disagreement at common thresholds
for (thresh in c(1, 2, 3)) {
  below_se  <- soil$CH4_SNR_se < thresh
  below_mdf <- soil$CH4_SNR_mdf < thresh
  disagree  <- sum(below_se != below_mdf, na.rm = TRUE)
  cat(sprintf("  At SNR threshold %d: SE flags %d, MDF flags %d, disagree on %d (%.1f%%)\n",
              thresh,
              sum(below_se, na.rm = TRUE),
              sum(below_mdf, na.rm = TRUE),
              disagree,
              100 * disagree / nrow(soil)))
}

p_snr_compare <- ggplot(soil, aes(x = CH4_SNR_mdf, y = CH4_SNR_se)) +
  geom_point(alpha = 0.4, size = 1.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  scale_x_log10() + scale_y_log10() +
  labs(x = "SNR (MDF-based): |flux| / MDF",
       y = "SNR (SE-based): |flux| / SE",
       title = "Soil: SNR comparison — SE-based vs MDF-based",
       subtitle = "Red = 1:1 line") +
  theme_bw(base_size = 12)

ggsave(file.path(soil_fig_dir, "snr_comparison.png"),
       p_snr_compare, width = 7, height = 5, dpi = 200)


# ============================================================================
# PART 6: SCALING UNCERTAINTY — Three approaches compared
# ============================================================================

message("\n=== Part 6: Scaling uncertainty ===")

flux <- soil$best.flux
n_flux <- length(flux)
se_vals <- soil$CH4_best_SE

# --- Approach 1: Simple SE of the mean ---
m <- mean(flux); se1 <- sd(flux) / sqrt(n_flux)
cat(sprintf("\n  1. Simple SE:      mean = %.4f, SE = %.4f, 95%% CI [%.4f, %.4f]\n",
            m, se1, m - 1.96 * se1, m + 1.96 * se1))

# --- Approach 2: Bootstrap CI ---
set.seed(42)
boot_means <- replicate(10000, mean(sample(flux, replace = TRUE)))
se2 <- sd(boot_means)
ci2 <- quantile(boot_means, probs = c(0.025, 0.975))
cat(sprintf("  2. Bootstrap:      mean = %.4f, SE = %.4f, 95%% CI [%.4f, %.4f]\n",
            mean(flux), se2, ci2[[1]], ci2[[2]]))

# --- Approach 3: Random-effects (SE-based) ---
vi <- se_vals^2
vi[vi < 1e-20 | is.na(vi)] <- 1e-20
fit <- tryCatch(rma(yi = flux, vi = vi, method = "REML"),
                error = function(e) NULL)
if (!is.null(fit)) {
  cat(sprintf("  3. Random effects: mean = %.4f, SE = %.4f, 95%% CI [%.4f, %.4f]\n",
              as.numeric(fit$beta), fit$se, fit$ci.lb, fit$ci.ub))
  cat(sprintf("     tau2 = %.4f, I2 = %.1f%%\n", fit$tau2, fit$I2))
} else {
  cat("  3. Random effects: FAILED\n")
}

cat(sprintf("\n  SE ratio (simple/boot): %.4f\n", se1 / se2))
cat("  (ratio ~ 1 confirms simple SE ≈ bootstrap)\n")

if (!is.null(fit)) {
  cat(sprintf("  SE ratio (simple/RE):   %.4f\n", se1 / fit$se))
  cat(sprintf("  SE ratio (RE/boot):     %.4f\n", fit$se / se2))
}

# --- Approach 3b: Random-effects using MDF² as within-measurement variance ---
vi_mdf <- soil$MDF^2
vi_mdf[vi_mdf < 1e-20 | is.na(vi_mdf)] <- 1e-20
fit_mdf <- tryCatch(rma(yi = flux, vi = vi_mdf, method = "REML"),
                     error = function(e) NULL)
if (!is.null(fit_mdf)) {
  cat(sprintf("\n  3b. RE (MDF²):     mean = %.4f, SE = %.4f, 95%% CI [%.4f, %.4f]\n",
              as.numeric(fit_mdf$beta), fit_mdf$se, fit_mdf$ci.lb, fit_mdf$ci.ub))
  cat(sprintf("      tau2 = %.4f, I2 = %.1f%%\n", fit_mdf$tau2, fit_mdf$I2))
}


# ============================================================================
# PART 7: VARIANCE DECOMPOSITION
# ============================================================================

message("\n=== Part 7: Variance decomposition ===")

if (!is.null(fit)) {
  mean_vi_se <- mean(vi)
  cat("\n--- Using goFlux SE² as measurement variance ---\n")
  cat(sprintf("  Total observed var:   %.4f\n", var(flux)))
  cat(sprintf("  tau2 (biological):    %.4f\n", fit$tau2))
  cat(sprintf("  mean(SE²) (meas):    %.6f\n", mean_vi_se))
  cat(sprintf("  %% biological:         %.1f%%\n",
              100 * fit$tau2 / (fit$tau2 + mean_vi_se)))
  cat(sprintf("  %% measurement:        %.1f%%\n",
              100 * mean_vi_se / (fit$tau2 + mean_vi_se)))
  cat(sprintf("  I²:                   %.1f%%\n", fit$I2))
}

if (!is.null(fit_mdf)) {
  mean_vi_mdf <- mean(vi_mdf)
  cat("\n--- Using MDF² as measurement variance ---\n")
  cat(sprintf("  tau2 (biological):    %.4f\n", fit_mdf$tau2))
  cat(sprintf("  mean(MDF²) (meas):   %.6f\n", mean_vi_mdf))
  cat(sprintf("  %% biological:         %.1f%%\n",
              100 * fit_mdf$tau2 / (fit_mdf$tau2 + mean_vi_mdf)))
  cat(sprintf("  %% measurement:        %.1f%%\n",
              100 * mean_vi_mdf / (fit_mdf$tau2 + mean_vi_mdf)))
  cat(sprintf("  I²:                   %.1f%%\n", fit_mdf$I2))
}

# Seasonal breakdown
cat("\n--- Seasonal variance decomposition (SE²-based) ---\n")
for (ssn in c("Growing", "Dormant")) {
  d <- soil %>% filter(season == ssn)
  if (nrow(d) < 5) next
  vi_s <- d$CH4_best_SE^2
  vi_s[vi_s < 1e-20] <- 1e-20
  fit_s <- tryCatch(rma(yi = d$best.flux, vi = vi_s, method = "REML"),
                     error = function(e) NULL)
  if (!is.null(fit_s)) {
    cat(sprintf("  %s (n=%d): mean=%.4f, SE=%.4f, tau2=%.4f, I2=%.1f%%\n",
                ssn, nrow(d), as.numeric(fit_s$beta), fit_s$se,
                fit_s$tau2, fit_s$I2))
  }
}


# ============================================================================
# PART 8: KEY COMPARISON — Soil (UGGA) vs Stem (GLA131 / LI-7810)
# ============================================================================

message("\n=== Part 8: Cross-dataset comparison ===")

cat("\n--- SE / flux_term ratio (diagnostic for unit consistency) ---\n")
cat("  If SE units are correct, SE / flux_term should be comparable\n")
cat("  to slope_SE = precision / sqrt(n(n²-1)/12)\n")

soil_diag <- soil_verify %>%
  mutate(slope_SE_obs = LM.SE / abs(flux.term))

cat(sprintf("\n  Soil (UGGA): median slope_SE = %.6f ppb/s\n",
            median(soil_diag$slope_SE_obs, na.rm = TRUE)))
cat(sprintf("  Soil:        median ratio (obs/exp) = %.3f\n",
            median(soil_verify$ratio_obs_exp_mfr, na.rm = TRUE)))
cat(sprintf("  (Compare: GLA131 stem ratio was 0.013, LI-7810 was 1.14)\n"))


# ============================================================================
# SAVE
# ============================================================================

message("\n=== Saving results ===")

save(soil, soil_verify,
     file = file.path(output_dir, "14_soil_se_analysis.RData"))

message("Saved: ", file.path(output_dir, "14_soil_se_analysis.RData"))
message("Figures: ", soil_fig_dir)
message("\n=== Soil SE analysis complete ===")
