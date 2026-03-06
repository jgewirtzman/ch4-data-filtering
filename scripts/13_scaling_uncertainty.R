# ============================================================================
# 13_scaling_uncertainty.R — Compare approaches for estimating uncertainty
#                             on mean component flux for scaling
#
# Three approaches:
#   1. Simple SE of the mean = SD(fluxes) / sqrt(n)
#   2. Bootstrap CI (resample measured values, compute mean per replicate)
#   3. Random-effects meta-analytic framework (rma from metafor)
#      — uses noise_floor² as within-measurement variance (vi)
#      — decomposes total variance into measurement (vi) + biological (τ²)
#
# Key findings from SE verification (script 12):
#   - GLA131 model-fit SE is unreliable (77× too small vs first principles)
#   - LI-7810 model-fit SE matches first principles (ratio = 1.14)
#   - Allan noise floor is the trustworthy per-measurement uncertainty metric
#   → Use noise_floor² (not SE²) as within-measurement variance in approach 3
#
# Loads: outputs/03_mdf_results.RData
# Saves: figures/scaling_uncertainty/
# ============================================================================

source("scripts/00_setup.R")

# metafor for random-effects meta-analysis
if (!requireNamespace("metafor", quietly = TRUE)) {
  message("Installing metafor package...")
  install.packages("metafor", repos = "https://cloud.r-project.org")
}
library(metafor)

message("\n=== Loading data ===")
load(file.path(output_dir, "03_mdf_results.RData"))

scale_fig_dir <- file.path(fig_dir, "scaling_uncertainty")
dir.create(scale_fig_dir, recursive = TRUE, showWarnings = FALSE)

# Add instrument label and season
df <- df %>%
  mutate(
    inst_label = if_else(year == 2025, "Analyzer B", "Analyzer A"),
    month  = as.integer(format(as.Date(date), "%m")),
    season = if_else(month %in% 5:10, "Growing", "Dormant")
  )

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

boot_ci <- function(x, n_boot = 10000, conf = 0.95) {
  set.seed(42)
  boot_means <- replicate(n_boot, mean(sample(x, replace = TRUE)))
  alpha <- (1 - conf) / 2
  ci <- quantile(boot_means, probs = c(alpha, 1 - alpha))
  list(mean = mean(x), se_boot = sd(boot_means),
       ci_lo = ci[[1]], ci_hi = ci[[2]])
}

# --- Approach 1: Simple SE of the mean ---
approach1 <- function(flux) {
  n <- length(flux); m <- mean(flux); se <- sd(flux) / sqrt(n)
  data.frame(approach = "1. Simple SE", mean = m, se = se,
             ci_lo = m - 1.96 * se, ci_hi = m + 1.96 * se,
             ci_width = 2 * 1.96 * se, tau2 = NA_real_, I2 = NA_real_)
}

# --- Approach 2: Bootstrap CI ---
approach2 <- function(flux) {
  b <- boot_ci(flux)
  data.frame(approach = "2. Bootstrap", mean = b$mean, se = b$se_boot,
             ci_lo = b$ci_lo, ci_hi = b$ci_hi,
             ci_width = b$ci_hi - b$ci_lo, tau2 = NA_real_, I2 = NA_real_)
}

# --- Approach 3: Random-effects (noise-floor-based) ---
approach3 <- function(flux, noise_floor_vals) {
  vi <- noise_floor_vals^2
  vi[vi < 1e-20 | is.na(vi)] <- 1e-20

  fit <- tryCatch(rma(yi = flux, vi = vi, method = "REML"),
                  error = function(e) NULL)
  if (is.null(fit)) {
    return(data.frame(approach = "3. Random effects (NF)", mean = NA, se = NA,
                      ci_lo = NA, ci_hi = NA, ci_width = NA,
                      tau2 = NA, I2 = NA))
  }
  data.frame(approach = "3. Random effects (NF)",
             mean = as.numeric(fit$beta), se = fit$se,
             ci_lo = fit$ci.lb, ci_hi = fit$ci.ub,
             ci_width = fit$ci.ub - fit$ci.lb,
             tau2 = fit$tau2, I2 = fit$I2)
}


# ============================================================================
# PART 1: MAIN COMPARISON — ALL SUBSETS
# ============================================================================

message("\n=== Part 1: Three-approach comparison across subsets ===")

# Restrict to measurements with noise_floor for fair comparison
df_nf <- df %>% filter(!is.na(CH4_noise_floor))
message("Using ", nrow(df_nf), " measurements with noise floor estimates")

run_all_three <- function(d, label) {
  rbind(
    approach1(d$CH4_flux_nmolpm2ps),
    approach2(d$CH4_flux_nmolpm2ps),
    approach3(d$CH4_flux_nmolpm2ps, d$CH4_noise_floor)
  ) %>% mutate(subset = label)
}

# Full dataset
results_all <- run_all_three(df_nf, "All instruments")

# By instrument
results_inst <- lapply(c("Analyzer A", "Analyzer B"), function(inst) {
  run_all_three(df_nf %>% filter(inst_label == inst), inst)
}) %>% do.call(rbind, .)

# By season × instrument
results_season <- lapply(c("Analyzer A", "Analyzer B"), function(inst) {
  lapply(c("Growing", "Dormant"), function(ssn) {
    d <- df_nf %>% filter(inst_label == inst, season == ssn)
    if (nrow(d) < 5) return(NULL)
    run_all_three(d, paste(inst, ssn))
  }) %>% do.call(rbind, .)
}) %>% do.call(rbind, .)

# Removal bias
results_removal <- lapply(c(FALSE, TRUE), function(remove) {
  if (remove) {
    d <- df_nf %>% filter(!is.na(CH4_below_MDF_chr95), !CH4_below_MDF_chr95)
    label <- "After removal (Chr 95%)"
  } else {
    d <- df_nf %>% filter(!is.na(CH4_below_MDF_chr95))
    label <- "All retained"
  }
  run_all_three(d, label)
}) %>% do.call(rbind, .)

all_results <- rbind(results_all, results_inst, results_season, results_removal)

cat("\n--- All results ---\n")
print(all_results, row.names = FALSE)


# ============================================================================
# PART 2: DIRECT COMPARISON — Simple SE vs Bootstrap
# ============================================================================

message("\n=== Part 2: Simple SE vs Bootstrap — are they equivalent? ===")

# Compare SE and CI width for every subset
se_vs_boot <- all_results %>%
  filter(approach %in% c("1. Simple SE", "2. Bootstrap")) %>%
  select(subset, approach, se, ci_width) %>%
  pivot_wider(names_from = approach, values_from = c(se, ci_width),
              names_sep = "__")

# Clean column names
names(se_vs_boot) <- gsub("1\\. Simple SE", "simple", names(se_vs_boot))
names(se_vs_boot) <- gsub("2\\. Bootstrap", "boot", names(se_vs_boot))

se_vs_boot <- se_vs_boot %>%
  mutate(
    se_ratio     = round(se__simple / se__boot, 4),
    ci_ratio     = round(ci_width__simple / ci_width__boot, 4)
  )

cat("\n--- Simple SE vs Bootstrap comparison ---\n")
print(as.data.frame(se_vs_boot), row.names = FALSE)
cat("\nRatio ≈ 1.0 means the approaches give equivalent uncertainty estimates.\n")
cat("Simple SE is therefore the preferred default (no resampling needed).\n")


# ============================================================================
# PART 3: VARIANCE DECOMPOSITION (τ² vs noise floor²)
# ============================================================================

message("\n=== Part 3: Variance decomposition ===")

var_decomp <- lapply(c("Analyzer A", "Analyzer B"), function(inst) {
  lapply(c("All", "Growing", "Dormant"), function(ssn) {
    if (ssn == "All") {
      d <- df_nf %>% filter(inst_label == inst)
    } else {
      d <- df_nf %>% filter(inst_label == inst, season == ssn)
    }
    if (nrow(d) < 5) return(NULL)

    vi <- d$CH4_noise_floor^2
    vi[vi < 1e-20] <- 1e-20

    fit <- tryCatch(rma(yi = d$CH4_flux_nmolpm2ps, vi = vi, method = "REML"),
                    error = function(e) NULL)
    if (is.null(fit)) return(NULL)

    mean_vi <- mean(vi)
    data.frame(
      instrument      = inst,
      season          = ssn,
      n               = nrow(d),
      total_var       = round(var(d$CH4_flux_nmolpm2ps), 4),
      tau2_biological = round(fit$tau2, 4),
      mean_nf2_meas   = round(mean_vi, 4),
      pct_biological  = round(100 * fit$tau2 / (fit$tau2 + mean_vi), 1),
      pct_measurement = round(100 * mean_vi / (fit$tau2 + mean_vi), 1),
      I2              = round(fit$I2, 1)
    )
  }) %>% do.call(rbind, .)
}) %>% do.call(rbind, .)

cat("\n--- Variance decomposition: biological (τ²) vs measurement (noise_floor²) ---\n")
print(var_decomp, row.names = FALSE)
cat("\nI² = % of observed heterogeneity beyond measurement noise\n")
cat("I² ≈ 100: biological variation dominates → simple SE is sufficient\n")
cat("I² << 100: measurement noise matters → per-measurement weighting may help\n")


# ============================================================================
# PART 4: REMOVAL BIAS ON UNCERTAINTY — DETAILED
# ============================================================================

message("\n=== Part 4: Removal bias on mean AND variance ===")

removal_detail <- lapply(c("No filter", "Chr 90%", "Chr 95%", "Chr 99%",
                            "R² > 0.5", "R² > 0.7"), function(filt_name) {
  if (filt_name == "No filter") {
    d <- df_nf
  } else if (grepl("Chr", filt_name)) {
    conf <- sub("Chr (\\d+)%", "\\1", filt_name)
    col <- paste0("CH4_below_MDF_chr", conf)
    d <- df_nf %>% filter(!is.na(.data[[col]]), !.data[[col]])
  } else {
    thresh <- as.numeric(sub("R² > ", "", filt_name))
    d <- df_nf %>% filter(CH4_r2 > thresh)
  }

  if (nrow(d) < 5) return(NULL)

  n <- nrow(d)
  m <- mean(d$CH4_flux_nmolpm2ps)
  se <- sd(d$CH4_flux_nmolpm2ps) / sqrt(n)

  data.frame(
    filter    = filt_name,
    n         = n,
    mean      = round(m, 4),
    se_mean   = round(se, 4),
    ci_lo     = round(m - 1.96 * se, 4),
    ci_hi     = round(m + 1.96 * se, 4),
    ci_width  = round(2 * 1.96 * se, 4),
    sd        = round(sd(d$CH4_flux_nmolpm2ps), 4)
  )
}) %>% do.call(rbind, .)

# Add percent bias columns
base <- removal_detail[removal_detail$filter == "No filter",]
removal_detail$mean_bias_pct <- round(
  100 * (removal_detail$mean - base$mean) / abs(base$mean), 1)
removal_detail$ci_width_ratio <- round(
  removal_detail$ci_width / base$ci_width, 2)

cat("\n--- Filter effect on mean AND CI width ---\n")
print(removal_detail, row.names = FALSE)
cat("\nci_width_ratio > 1: CI wider after removal (more variance in retained subset)\n")
cat("ci_width_ratio < 1: CI narrower (fewer measurements → less power)\n")


# ============================================================================
# PART 5: FOREST PLOT VISUALIZATION
# ============================================================================

message("\n=== Part 5: Visualization ===")

all_results$subset <- factor(all_results$subset,
                              levels = rev(unique(all_results$subset)))

p_forest <- ggplot(all_results, aes(x = mean, y = subset,
                                     color = approach, shape = approach)) +
  geom_pointrange(aes(xmin = ci_lo, xmax = ci_hi),
                  position = position_dodge(width = 0.5), size = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  labs(x = expression("Mean CH"[4]~"flux (nmol"~m^{-2}~s^{-1}*") with 95% CI"),
       y = NULL,
       title = "Scaling uncertainty: three approaches compared",
       subtitle = "Random effects uses noise floor² as within-measurement variance",
       color = "Approach", shape = "Approach") +
  theme_bw(base_size = 12) +
  theme(legend.position = "top")

ggsave(file.path(scale_fig_dir, "forest_plot_comparison.png"),
       p_forest, width = 10, height = 8, dpi = 200)
message("  Saved: forest_plot_comparison.png")

# Removal bias visualization
removal_detail$filter <- factor(removal_detail$filter,
                                 levels = rev(removal_detail$filter))
p_removal <- ggplot(removal_detail, aes(x = mean, y = filter)) +
  geom_pointrange(aes(xmin = ci_lo, xmax = ci_hi), size = 0.5) +
  geom_vline(xintercept = base$mean, linetype = "dashed", color = "blue") +
  geom_text(aes(label = paste0(ifelse(mean_bias_pct > 0, "+", ""),
                                mean_bias_pct, "%")),
            hjust = -0.3, size = 3.5, color = "red") +
  labs(x = expression("Mean CH"[4]~"flux ± 95% CI (nmol"~m^{-2}~s^{-1}*")"),
       y = NULL,
       title = "Filter removal shifts both the mean and the CI",
       subtitle = "Blue dashed = unfiltered mean; red = % bias") +
  theme_bw(base_size = 12)

ggsave(file.path(scale_fig_dir, "removal_bias_with_ci.png"),
       p_removal, width = 9, height = 5, dpi = 200)
message("  Saved: removal_bias_with_ci.png")


# ============================================================================
# SAVE ALL RESULTS
# ============================================================================

message("\n=== Saving results ===")

save(all_results, se_vs_boot, var_decomp, removal_detail,
     file = file.path(output_dir, "13_scaling_uncertainty.RData"))

write.csv(all_results,
          file.path(output_dir, "scaling_uncertainty_comparison.csv"),
          row.names = FALSE)

message("Saved: ", file.path(output_dir, "13_scaling_uncertainty.RData"))
message("\n=== Scaling uncertainty analysis complete ===")
