# ============================================================================
# 08_transformation_comparison.R — Compare raw, log, and arcsinh transforms
#
# Demonstrates that log transforms silently drop negative values (true CH4
# uptake), while arcsinh handles the full range including negative fluxes.
# Uses three datasets: wetland soil, tree stem (YMF), tree canopy (HF).
#
# Builds the plot data (transformed values, per-panel counts, arcsinh curve)
# drawn in 19_main_figures.R, and runs Shapiro-Wilk normality tests and a
# variance stabilization comparison (console).
#
# Loads: outputs/01_loaded_data.RData, outputs/03_mdf_results.RData
# Saves: outputs/08_transformation.RData
# ============================================================================

source("scripts/00_setup.R")

message("\n=== Loading data from 01 and 03 ===")
load(file.path(output_dir, "01_loaded_data.RData"))
load(file.path(output_dir, "03_mdf_results.RData"))

# ── Load the three CH4 flux datasets ────────────────────────────────────────

# 1. Soils along an upland-wetland gradient
soil <- soil_flux %>%
  filter(measurement_type == "soil", !is.na(CH4_best_flux_nmol_m2_s)) %>%
  mutate(flux = CH4_best_flux_nmol_m2_s, dataset = "Upland–wetland soils")

# 2. Forest stem fluxes (both analyzers)
tree_stem <- df %>%
  filter(is.finite(CH4_flux_nmolpm2ps)) %>%
  mutate(flux = CH4_flux_nmolpm2ps, dataset = "Forest stems")

# 3. Harvard Forest canopy fluxes
hf_flux <- hf_out %>%
  filter(!is.na(CH4_best.flux)) %>%
  mutate(flux = CH4_best.flux, dataset = "Tree canopy (HF)")

# Combine (two datasets only; canopy dropped for clarity)
all_data <- bind_rows(
  soil %>% select(flux, dataset),
  tree_stem %>% select(flux, dataset)
)

# ── Summary statistics ───────────────────────────────────────────────────────

cat("\n=== Data summary ===\n")
all_data %>%
  group_by(dataset) %>%
  summarise(
    n = n(),
    n_negative = sum(flux < 0),
    pct_negative = round(100 * mean(flux < 0), 1),
    min = round(min(flux), 3),
    median = round(median(flux), 4),
    mean = round(mean(flux), 4),
    max = round(max(flux), 3),
    .groups = "drop"
  ) %>%
  print(width = 120)

# ── Apply transformations ────────────────────────────────────────────────────

transform_data <- all_data %>%
  mutate(
    raw = flux,
    log10 = ifelse(flux > 0, log10(flux), NA),  # NA for non-positive
    asinh = asinh(flux)
  )

# Count what log loses
log_loss <- transform_data %>%
  group_by(dataset) %>%
  summarise(
    n_total = n(),
    n_lost_log = sum(is.na(log10)),
    pct_lost = round(100 * mean(is.na(log10)), 1),
    .groups = "drop"
  )

cat("\n=== Values lost by log10 transform ===\n")
print(log_loss)

# ── Transformed values by dataset (raw, log10, arcsinh) ─────────────────────

transform_long <- transform_data %>%
  pivot_longer(
    cols = c(raw, log10, asinh),
    names_to = "transform",
    values_to = "value"
  ) %>%
  mutate(
    transform = factor(transform,
                       levels = c("raw", "log10", "asinh"),
                       labels = c("Raw (untransformed)",
                                  "log10 (drops negatives)",
                                  "arcsinh (preserves all)")),
    dataset = factor(dataset,
                     levels = c("Upland–wetland soils", "Forest stems"))
  )

# Count retained vs dropped per panel
panel_counts <- transform_long %>%
  group_by(dataset, transform) %>%
  summarise(
    n_total = n(),
    n_valid = sum(!is.na(value)),
    n_dropped = sum(is.na(value)),
    .groups = "drop"
  ) %>%
  mutate(
    label = ifelse(n_dropped > 0,
                   paste0("n = ", n_valid, " / ", n_total,
                          "\n(", n_dropped, " dropped, ",
                          round(100 * n_dropped / n_total), "%)"),
                   paste0("n = ", n_valid))
  )

# Create combined facet label (dataset x transform)
transform_long <- transform_long %>%
  mutate(
    panel_label = paste0(dataset, "\n", transform),
    panel_label = factor(panel_label,
      levels = as.vector(t(outer(
        levels(dataset),
        levels(transform),
        function(d, t) paste0(d, "\n", t)
      )))
    )
  )

# Update panel_counts to match
panel_counts <- panel_counts %>%
  mutate(
    panel_label = paste0(dataset, "\n", transform),
    panel_label = factor(panel_label, levels = levels(transform_long$panel_label))
  )

# Compute back-transformed mean and median for each panel
# Back-transform each value individually, then take mean/median
bt_stats <- transform_long %>%
  filter(!is.na(value)) %>%
  mutate(
    bt_value = case_when(
      grepl("Raw", transform)     ~ value,
      grepl("log10", transform)   ~ 10^value,
      grepl("arcsinh", transform) ~ sinh(value)
    )
  ) %>%
  group_by(dataset, transform, panel_label) %>%
  summarise(
    bt_mean   = mean(bt_value),
    bt_median = median(bt_value),
    .groups = "drop"
  )

cat("\n=== Back-transformed summary statistics ===\n")
bt_stats %>% print(width = 100)

# Merge back-transformed stats into panel_counts for a single top-right annotation
panel_counts <- panel_counts %>%
  left_join(bt_stats %>% select(dataset, transform, bt_mean, bt_median),
            by = c("dataset", "transform")) %>%
  mutate(
    label = ifelse(n_dropped > 0,
                   paste0("n = ", n_valid, " / ", n_total,
                          " (", n_dropped, " dropped, ",
                          round(100 * n_dropped / n_total), "%)",
                          "\nmean = ", round(bt_mean, 3),
                          "  med = ", round(bt_median, 3)),
                   paste0("n = ", n_valid,
                          "\nmean = ", round(bt_mean, 3),
                          "  med = ", round(bt_median, 3)))
  )

# ── Arcsinh behavior curve (plot data for 19_main_figures.R) ───────────────
# arcsinh(x) ~ x near zero and ~ sign(x)*log(2|x|) for large |x|

x <- seq(-200, 200, length.out = 1000)
demo <- data.frame(
  x = x,
  asinh = asinh(x),
  log_approx = sign(x) * log(2 * abs(x)),
  linear = x
)

# ── Normality tests ──────────────────────────────────────────────────────────

cat("\n=== Shapiro-Wilk normality tests (wetland soil CH4) ===\n")
soil_flux_vals <- soil$flux

# Raw
sw_raw <- shapiro.test(sample(soil_flux_vals, min(5000, length(soil_flux_vals))))
cat(sprintf("Raw:     W = %.4f, p = %.2e\n", sw_raw$statistic, sw_raw$p.value))

# Log (positive only)
pos <- soil_flux_vals[soil_flux_vals > 0]
if (length(pos) > 3) {
  sw_log <- shapiro.test(log10(pos))
  cat(sprintf("log10:   W = %.4f, p = %.2e (n = %d of %d; %d%% of data dropped)\n",
              sw_log$statistic, sw_log$p.value, length(pos), length(soil_flux_vals),
              round(100 * (1 - length(pos)/length(soil_flux_vals)))))
}

# Arcsinh
sw_asinh <- shapiro.test(asinh(soil_flux_vals))
cat(sprintf("arcsinh: W = %.4f, p = %.2e (n = %d, all values retained)\n",
            sw_asinh$statistic, sw_asinh$p.value, length(soil_flux_vals)))

# ── Variance homogeneity demonstration ───────────────────────────────────────

cat("\n=== Variance stabilization (all datasets) ===\n")
all_data %>%
  group_by(dataset) %>%
  summarise(
    raw_var = var(flux),
    raw_sd = sd(flux),
    asinh_var = var(asinh(flux)),
    asinh_sd = sd(asinh(flux)),
    var_ratio_raw = var(flux) / min(var(flux)),  # within-group; cross-group below
    .groups = "drop"
  ) %>%
  print()

# Cross-dataset variance ratio
vars_raw <- all_data %>% group_by(dataset) %>% summarise(v = var(flux), .groups = "drop")
vars_asinh <- all_data %>% group_by(dataset) %>% summarise(v = var(asinh(flux)), .groups = "drop")
cat(sprintf("\nRaw variance ratio (max/min across datasets): %.1f\n",
            max(vars_raw$v) / min(vars_raw$v)))
cat(sprintf("arcsinh variance ratio (max/min across datasets): %.1f\n",
            max(vars_asinh$v) / min(vars_asinh$v)))

message("\n=== Transformation comparison complete ===")

# Plot-ready data for the assembled main figures (19_main_figures.R)
save(transform_long, panel_counts, demo, file = file.path(output_dir, "08_transformation.RData"))
