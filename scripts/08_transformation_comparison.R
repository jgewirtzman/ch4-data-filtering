# ============================================================================
# 08_transformation_comparison.R — Compare raw, log, and arcsinh transforms
#
# Demonstrates that log transforms silently drop negative values (true CH4
# uptake), while arcsinh handles the full range including negative fluxes.
# Uses three datasets: wetland soil, tree stem (YMF), tree canopy (HF).
#
# Figures:
#   1. 3x3 faceted histograms (datasets x transforms)
#   2. Q-Q plots for wetland soil under all 3 transforms
#   3. Arcsinh behavior demo (full range + zoom near zero)
#
# Also runs Shapiro-Wilk normality tests and variance stabilization comparison.
#
# Loads: outputs/01_loaded_data.RData, outputs/03_mdf_results.RData
# Saves: figures to fig_dir (PNG + PDF)
# ============================================================================

source("scripts/00_setup.R")

message("\n=== Loading data from 01 and 03 ===")
load(file.path(output_dir, "01_loaded_data.RData"))
load(file.path(output_dir, "03_mdf_results.RData"))

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ── Load the three CH4 flux datasets ────────────────────────────────────────

# 1. Wetland soil CH4 fluxes (semirigid chambers -- soil measurements)
soil <- soil_flux %>%
  filter(measurement_type == "soil", !is.na(CH4_best_flux_nmol_m2_s)) %>%
  mutate(flux = CH4_best_flux_nmol_m2_s, dataset = "Wetland soil")

# 2. Tree stem fluxes from same dataset
tree_stem <- soil_flux %>%
  filter(measurement_type == "tree_stem", !is.na(CH4_best_flux_nmol_m2_s)) %>%
  mutate(flux = CH4_best_flux_nmol_m2_s, dataset = "Tree stem (YMF)")

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

# ── FIGURE 1: Three-panel transformation comparison ─────────────────────────
# Each dataset gets a row; columns are raw, log10, arcsinh

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
                     levels = c("Wetland soil", "Tree stem (YMF)"))
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

# Main figure: faceted histograms
# Use facet_wrap (not facet_grid) so each panel gets independent y-axes,
# since sample sizes differ dramatically across dataset-transform combos
# (e.g. log10 drops 83% of wetland soil but only 10% of canopy)

# Create combined facet label
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

p_hist <- ggplot(transform_long, aes(x = value)) +
  geom_histogram(aes(fill = dataset), bins = 40, alpha = 0.8, color = "grey30", linewidth = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.4) +
  facet_wrap(~ panel_label, ncol = 3, scales = "free") +
  geom_text(
    data = panel_counts,
    aes(label = label),
    x = Inf, y = Inf, hjust = 1.05, vjust = 1.1,
    size = 2.6, color = "grey20", inherit.aes = FALSE
  ) +
  scale_fill_manual(values = c(
    "Wetland soil" = "#2166AC",
    "Tree stem (YMF)" = "#4DAF4A"
  )) +
  scale_y_sqrt() +
  labs(
    x = "CH4 flux (transformed value)",
    y = "Count (sqrt scale)",
    title = "Effect of data transformation on CH4 flux distributions",
    subtitle = "log10 silently drops negative values; arcsinh preserves the full distribution\nmean/median back-transformed to original units (nmol m\u207b\u00b2 s\u207b\u00b9)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 9, face = "bold"),
    plot.title = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 10, color = "grey30"),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(fig_dir, "ch4_transform_comparison.png"), p_hist,
       width = 12, height = 6, dpi = 300)
ggsave(file.path(fig_dir, "ch4_transform_comparison.pdf"), p_hist,
       width = 12, height = 6)
cat("\nSaved: ch4_transform_comparison.png/pdf\n")

# ── FIGURE 2: QQ-plot comparison ────────────────────────────────────────────
# Show how arcsinh normalizes better than raw or log

qq_data <- transform_data %>%
  filter(dataset == "Wetland soil") %>%
  select(flux, raw, log10, asinh) %>%
  pivot_longer(cols = c(raw, log10, asinh), names_to = "transform", values_to = "value") %>%
  filter(!is.na(value)) %>%
  mutate(
    transform = factor(transform,
                       levels = c("raw", "log10", "asinh"),
                       labels = c("Raw", "log10", "arcsinh"))
  )

# Compute theoretical quantiles per transform
qq_data <- qq_data %>%
  group_by(transform) %>%
  arrange(value) %>%
  mutate(
    theoretical = qnorm(ppoints(n()))
  ) %>%
  ungroup()

p_qq <- ggplot(qq_data, aes(x = theoretical, y = value)) +
  geom_point(alpha = 0.5, size = 1, color = "#2166AC") +
  geom_smooth(method = "lm", se = FALSE, color = "red", linewidth = 0.7) +
  facet_wrap(~ transform, scales = "free_y", ncol = 3) +
  labs(
    x = "Theoretical quantiles (normal)",
    y = "Sample quantiles",
    title = "Q-Q plots: wetland soil CH4 flux under three transformations",
    subtitle = "arcsinh produces a more symmetric, approximately normal distribution"
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.text = element_text(size = 11, face = "bold"),
    plot.title = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 10, color = "grey30")
  )

ggsave(file.path(fig_dir, "ch4_transform_qq.png"), p_qq,
       width = 10, height = 4, dpi = 300)
ggsave(file.path(fig_dir, "ch4_transform_qq.pdf"), p_qq,
       width = 10, height = 4)
cat("Saved: ch4_transform_qq.png/pdf\n")

# ── FIGURE 3: Arcsinh behavior demonstration ───────────────────────────────
# Show that arcsinh(x) ~ x near zero and ~ sign(x)*log(2|x|) for large |x|

x <- seq(-200, 200, length.out = 1000)
demo <- data.frame(
  x = x,
  asinh = asinh(x),
  log_approx = sign(x) * log(2 * abs(x)),
  linear = x
)

# Zoomed near zero
demo_zoom <- demo %>% filter(abs(x) <= 5)

p_behavior <- ggplot(demo, aes(x = x)) +
  geom_line(aes(y = asinh, color = "arcsinh(x)"), linewidth = 1) +
  geom_line(aes(y = log_approx, color = "sign(x) * log(2|x|)"),
            linewidth = 0.7, linetype = "dashed") +
  scale_color_manual(values = c("arcsinh(x)" = "#2166AC",
                                "sign(x) * log(2|x|)" = "#D6604D")) +
  labs(
    x = "x", y = "f(x)", color = NULL,
    title = "arcsinh: linear near zero, logarithmic for large values"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = c(0.2, 0.85),
    legend.background = element_rect(fill = "white", color = "grey70"),
    plot.title = element_text(size = 12, face = "bold")
  )

p_zoom <- ggplot(demo_zoom, aes(x = x)) +
  geom_line(aes(y = asinh, color = "arcsinh(x)"), linewidth = 1) +
  geom_line(aes(y = x, color = "y = x (linear)"), linewidth = 0.7, linetype = "dashed") +
  geom_vline(xintercept = 0, color = "grey60", linewidth = 0.3) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.3) +
  scale_color_manual(values = c("arcsinh(x)" = "#2166AC",
                                "y = x (linear)" = "#4DAF4A")) +
  labs(
    x = "x", y = "f(x)", color = NULL,
    title = "Near zero: arcsinh is approximately linear"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = c(0.25, 0.85),
    legend.background = element_rect(fill = "white", color = "grey70"),
    plot.title = element_text(size = 12, face = "bold")
  )

p_combined <- p_behavior + p_zoom + plot_layout(ncol = 2)

ggsave(file.path(fig_dir, "arcsinh_behavior.png"), p_combined,
       width = 10, height = 4.5, dpi = 300)
ggsave(file.path(fig_dir, "arcsinh_behavior.pdf"), p_combined,
       width = 10, height = 4.5)
cat("Saved: arcsinh_behavior.png/pdf\n")

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
