# ============================================================================
# 05_instrument_comparison.R — Stem flux ridgeline plot (all quality filters)
# and instrument-specific analyses (Analyzer A vs Analyzer B)
#
# PART 9: Main ggridges rainfall plot showing CH4 flux distribution under
#         ALL 13+ quality filter criteria, rows ordered by % retained,
#         with annotations (n and %)
# PART 10: Instrument-specific analyses:
#          - Faceted ridgeline by instrument (Analyzer A vs Analyzer B)
#          - Violin plot of Allan deviation by instrument
#          - Density comparison overlay
#
# Loads: outputs/03_mdf_results.RData
# Saves: figures to fig_dir (PNG + PDF)
# ============================================================================

source("scripts/00_setup.R")

message("\n=== Loading MDF results from 03 ===")
# Load 02 first so that 03's df (with MDF columns) takes precedence
load(file.path(output_dir, "02_allan_deviation.RData"))
load(file.path(output_dir, "03_mdf_results.RData"))

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# PART 9: GGRIDGES RAINFALL PLOT
# ============================================================

message("\n=== Building main ridgeline plot (all filters) ===")

asinh_trans <- trans_new(
  name = "asinh",
  transform = asinh,
  inverse = sinh,
  breaks = function(x) {
    rng <- sinh(x)
    pretty(rng, n = 8)
  }
)

p_ridges <- ggplot(ridge_df, aes(x = CH4_flux, y = filter)) +
  geom_density_ridges(
    aes(fill = pct_retained, point_color = ifelse(passes, "Retained", "Excluded")),
    jittered_points = TRUE,
    point_size = 0.6, point_alpha = 0.4,
    scale = 0.9, alpha = 0.7,
    position = position_raincloud(width = 0.05, ygap = 0.05),
    bandwidth = 0.3,
    from = asinh(-5), to = asinh(30)
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red", linewidth = 0.5) +
  geom_point(data = filter_stats, aes(x = mean_flux, y = filter),
             shape = "|", size = 4, color = "red", inherit.aes = FALSE) +
  scale_x_continuous(trans = asinh_trans,
                     breaks = c(-5, -2, -1, 0, 0.5, 1, 2, 5, 10, 20),
                     labels = c("-5", "-2", "-1", "0", "0.5", "1", "2",
                                "5", "10", "20")) +
  scale_fill_gradient2(low = "#B2182B", mid = "#F7F7F7", high = "#2166AC",
                       midpoint = 70, limits = c(20, 100),
                       name = "% retained") +
  scale_discrete_manual("point_color",
                        values = c("Retained" = "grey40", "Excluded" = "red"),
                        name = "Measurement") +
  labs(
    x = expression(CH[4]~flux~(nmol~m^{-2}~s^{-1})),
    y = NULL,
    title = expression(CH[4]~flux~distribution~under~quality~filters),
    subtitle = paste0("n = ", n_total,
                      " measurements (", sum(df$year < 2025), " Analyzer A + ",
                      sum(df$year == 2025), " Analyzer B)",
                      "\nAllan deviation coverage: ", n_with_allan, "/", n_total,
                      " | Labels: (n retained, % | % negative)")
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.y = element_text(size = 7, lineheight = 0.85),
    legend.position = "bottom",
    plot.title = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 9, color = "grey40")
  ) +
  coord_cartesian(xlim = c(asinh(-5), asinh(25)))

# ============================================================
# SAVE RIDGELINE PLOT
# ============================================================

ggsave(file.path(fig_dir, "ch4_flux_ridges_by_filter.png"),
       p_ridges, width = 13, height = 12, dpi = 300)
ggsave(file.path(fig_dir, "ch4_flux_ridges_by_filter.pdf"),
       p_ridges, width = 13, height = 12)

message("Saved ridges plot to: ", fig_dir)

# ============================================================
# SUMMARY TABLE
# ============================================================

cat("\n=== % Negative CH4 fluxes by filter ===\n")
neg_summary <- data.frame(
  filter = filter_order,
  stringsAsFactors = FALSE
) %>%
  rowwise() %>%
  mutate(
    n_pass = {
      fails <- quality_filters[[filter]](df)
      sum(!fails)
    },
    n_neg = {
      fails <- quality_filters[[filter]](df)
      sum(df$CH4_flux_nmolpm2ps[!fails] < 0)
    },
    pct_neg = round(100 * n_neg / n_pass, 1)
  ) %>%
  ungroup()

print(as.data.frame(neg_summary), row.names = FALSE)

# ============================================================
# ALLAN DEVIATION COMPARISON
# ============================================================

cat("\n=== Allan deviation comparison by instrument ===\n")
if (nrow(allan_df_lgr) > 0) {
  cat("Analyzer A:\n")
  cat("  CH4: median =", round(median(allan_df_lgr$allan_sd_CH4, na.rm = TRUE), 4),
      "ppb, IQR =", round(quantile(allan_df_lgr$allan_sd_CH4, 0.25, na.rm = TRUE), 4),
      "-", round(quantile(allan_df_lgr$allan_sd_CH4, 0.75, na.rm = TRUE), 4), "\n")
  cat("  CO2: median =", round(median(allan_df_lgr$allan_sd_CO2, na.rm = TRUE), 3),
      "ppm, IQR =", round(quantile(allan_df_lgr$allan_sd_CO2, 0.25, na.rm = TRUE), 3),
      "-", round(quantile(allan_df_lgr$allan_sd_CO2, 0.75, na.rm = TRUE), 3), "\n")
  cat("  Manufacturer spec: CH4 =", PREC_CH4_LGR, "ppb, CO2 =", PREC_CO2_LGR, "ppm\n")
}
if (nrow(allan_df_7810) > 0) {
  cat("Analyzer B:\n")
  cat("  CH4: median =", round(median(allan_df_7810$allan_sd_CH4, na.rm = TRUE), 4),
      "ppb, IQR =", round(quantile(allan_df_7810$allan_sd_CH4, 0.25, na.rm = TRUE), 4),
      "-", round(quantile(allan_df_7810$allan_sd_CH4, 0.75, na.rm = TRUE), 4), "\n")
  cat("  CO2: median =", round(median(allan_df_7810$allan_sd_CO2, na.rm = TRUE), 3),
      "ppm, IQR =", round(quantile(allan_df_7810$allan_sd_CO2, 0.25, na.rm = TRUE), 3),
      "-", round(quantile(allan_df_7810$allan_sd_CO2, 0.75, na.rm = TRUE), 3), "\n")
  cat("  Manufacturer spec: CH4 =", PREC_CH4_7810, "ppb, CO2 =", PREC_CO2_7810, "ppm\n")
}

# ============================================================
# PART 10: INSTRUMENT COMPARISON
# ============================================================

message("\n=== Instrument comparison ===")

# Label each measurement by instrument
df <- df %>%
  mutate(
    inst_label = case_when(
      year == 2025 ~ "Analyzer B",
      year < 2025  ~ "Analyzer A"
    )
  )

# --- Negative flux breakdown by instrument ---
cat("\n=== Negative CH4 fluxes by instrument ===\n")
inst_neg <- df %>%
  group_by(inst_label) %>%
  summarise(
    n_total   = n(),
    n_neg     = sum(CH4_flux_nmolpm2ps < 0),
    pct_neg   = round(100 * n_neg / n_total, 1),
    mean_flux = round(mean(CH4_flux_nmolpm2ps), 4),
    median_flux = round(median(CH4_flux_nmolpm2ps), 4),
    min_flux  = round(min(CH4_flux_nmolpm2ps), 4),
    max_flux  = round(max(CH4_flux_nmolpm2ps), 4),
    .groups   = "drop"
  )
print(as.data.frame(inst_neg), row.names = FALSE)

# --- Per-instrument filter sensitivity ---
cat("\n=== % Negative by filter, split by instrument ===\n")
inst_filter_summary <- lapply(filter_order, function(fn) {
  fails <- quality_filters[[fn]](df)
  df_pass <- df[!fails, ]
  df_pass %>%
    group_by(inst_label) %>%
    summarise(
      filter  = fn,
      n_pass  = n(),
      n_neg   = sum(CH4_flux_nmolpm2ps < 0),
      pct_neg = round(100 * n_neg / n_pass, 1),
      .groups = "drop"
    )
})
inst_filter_df <- do.call(rbind, inst_filter_summary)
# Pivot wide for readability
inst_wide <- inst_filter_df %>%
  select(filter, inst_label, n_pass, pct_neg) %>%
  pivot_wider(
    names_from  = inst_label,
    values_from = c(n_pass, pct_neg),
    names_glue  = "{inst_label}_{.value}"
  )
print(as.data.frame(inst_wide), row.names = FALSE)

# --- PLOT: Filter sensitivity by instrument (faceted ridges) ---

# Build long data with instrument facet
rows_inst <- lapply(seq_along(filter_order), function(i) {
  fn <- filter_order[i]
  fails <- quality_filters[[fn]](df)
  data.frame(
    filter     = fn,
    CH4_flux   = df$CH4_flux_nmolpm2ps,
    passes     = !fails,
    instrument = df$inst_label,
    stringsAsFactors = FALSE
  )
})
ridge_inst <- do.call(rbind, rows_inst)

# Compute per-filter per-instrument labels
ridge_inst <- ridge_inst %>%
  group_by(filter, instrument) %>%
  mutate(
    n_pass = sum(passes),
    n_inst = n(),
    pct_neg = round(100 * sum(CH4_flux[passes] < 0) / sum(passes), 1)
  ) %>%
  ungroup()

# Simplify filter labels for faceted plot
ridge_inst$filter_label <- sapply(ridge_inst$filter, function(fn) {
  paste0(fn, " (", ridge_inst$n_pass[ridge_inst$filter == fn][1], ")")
})
# Actually, make clean per-instrument per-filter labels
inst_stats <- ridge_inst %>%
  filter(passes) %>%
  group_by(filter, instrument) %>%
  summarise(
    n_pass  = n(),
    n_neg   = sum(CH4_flux < 0),
    pct_neg = round(100 * n_neg / n_pass, 1),
    mean_flux = mean(CH4_flux),
    .groups = "drop"
  )

# Order filters by total retention
filter_order_fct <- factor(ridge_inst$filter,
                            levels = rev(filter_order))
ridge_inst$filter <- filter_order_fct

p_inst_ridges <- ggplot(
  ridge_inst %>% filter(passes),
  aes(x = CH4_flux, y = filter, fill = instrument)
) +
  geom_density_ridges(
    alpha = 0.5, scale = 0.9,
    bandwidth = 0.3,
    from = asinh(-5), to = asinh(30)
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red",
             linewidth = 0.5) +
  scale_x_continuous(
    trans = asinh_trans,
    breaks = c(-5, -2, -1, 0, 0.5, 1, 2, 5, 10, 20),
    labels = c("-5", "-2", "-1", "0", "0.5", "1", "2", "5", "10", "20")
  ) +
  scale_fill_manual(
    values = c("Analyzer A" = "#E69F00", "Analyzer B" = "#56B4E9"),
    name = "Instrument"
  ) +
  labs(
    x = expression(CH[4]~flux~(nmol~m^{-2}~s^{-1})),
    y = NULL,
    title = expression(CH[4]~flux~distribution~by~instrument~and~quality~filter),
    subtitle = paste0("Overlaid densities: Analyzer A (n=", sum(df$year < 2025),
                      ") vs Analyzer B (n=", sum(df$year == 2025), ")")
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.y = element_text(size = 7),
    legend.position = "bottom",
    plot.title = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 9, color = "grey40")
  ) +
  coord_cartesian(xlim = c(asinh(-5), asinh(25)))

ggsave(file.path(fig_dir, "ch4_flux_ridges_by_instrument.png"),
       p_inst_ridges, width = 13, height = 12, dpi = 300)
ggsave(file.path(fig_dir, "ch4_flux_ridges_by_instrument.pdf"),
       p_inst_ridges, width = 13, height = 12)
message("Saved instrument comparison ridges to: ", fig_dir)

# --- PLOT: Allan deviation by instrument (violin + boxplot) ---

allan_combined <- bind_rows(
  allan_df_lgr %>%
    transmute(instrument = "Analyzer A",
              CH4_allan_sd = allan_sd_CH4,
              CO2_allan_sd = allan_sd_CO2),
  allan_df_7810 %>%
    transmute(instrument = "Analyzer B",
              CH4_allan_sd = allan_sd_CH4,
              CO2_allan_sd = allan_sd_CO2)
)

p_allan <- ggplot(allan_combined, aes(x = instrument, y = CH4_allan_sd,
                                       fill = instrument)) +
  geom_violin(alpha = 0.4, draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_jitter(width = 0.15, alpha = 0.15, size = 0.8) +
  geom_hline(yintercept = PREC_CH4_LGR,  linetype = "dashed", color = "#E69F00") +
  geom_hline(yintercept = PREC_CH4_7810, linetype = "dashed", color = "#56B4E9") +
  annotate("text", x = 2.4, y = PREC_CH4_LGR,
           label = paste0("Analyzer A spec: ", PREC_CH4_LGR, " ppb"),
           hjust = 1, size = 3, color = "#E69F00") +
  annotate("text", x = 2.4, y = PREC_CH4_7810,
           label = paste0("Analyzer B spec: ", PREC_CH4_7810, " ppb"),
           hjust = 1, size = 3, color = "#56B4E9") +
  scale_y_log10() +
  scale_fill_manual(
    values = c("Analyzer A" = "#E69F00", "Analyzer B" = "#56B4E9"),
    guide = "none"
  ) +
  labs(
    x = NULL,
    y = expression(CH[4]~Allan~deviation~(ppb)),
    title = expression(Per-measurement~CH[4]~Allan~deviation~by~instrument),
    subtitle = "Dashed lines = manufacturer specs | Violin shows full distribution"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 9, color = "grey40")
  )

ggsave(file.path(fig_dir, "ch4_allan_deviation_by_instrument.png"),
       p_allan, width = 7, height = 6, dpi = 300)
ggsave(file.path(fig_dir, "ch4_allan_deviation_by_instrument.pdf"),
       p_allan, width = 7, height = 6)
message("Saved Allan deviation comparison to: ", fig_dir)

# --- PLOT: Flux distributions by instrument (simple density comparison) ---

p_flux_inst <- ggplot(df, aes(x = CH4_flux_nmolpm2ps, fill = inst_label,
                                color = inst_label)) +
  geom_density(alpha = 0.3, linewidth = 0.8) +
  geom_rug(aes(color = inst_label), alpha = 0.1, sides = "b") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  scale_x_continuous(
    trans = asinh_trans,
    breaks = c(-5, -2, -1, 0, 0.5, 1, 2, 5, 10, 20),
    labels = c("-5", "-2", "-1", "0", "0.5", "1", "2", "5", "10", "20")
  ) +
  scale_fill_manual(
    values = c("Analyzer A" = "#E69F00", "Analyzer B" = "#56B4E9"),
    name = "Instrument"
  ) +
  scale_color_manual(
    values = c("Analyzer A" = "#E69F00", "Analyzer B" = "#56B4E9"),
    name = "Instrument"
  ) +
  labs(
    x = expression(CH[4]~flux~(nmol~m^{-2}~s^{-1})),
    y = "Density",
    title = expression(CH[4]~flux~distribution~by~instrument),
    subtitle = paste0(
      "Analyzer A: n=", sum(df$year < 2025),
      " (", round(100 * mean(df$CH4_flux_nmolpm2ps[df$year < 2025] < 0), 1), "% neg)",
      " | Analyzer B: n=", sum(df$year == 2025),
      " (", round(100 * mean(df$CH4_flux_nmolpm2ps[df$year == 2025] < 0), 1), "% neg)"
    )
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 10, color = "grey40")
  ) +
  coord_cartesian(xlim = c(asinh(-5), asinh(25)))

ggsave(file.path(fig_dir, "ch4_flux_density_by_instrument.png"),
       p_flux_inst, width = 9, height = 5, dpi = 300)
ggsave(file.path(fig_dir, "ch4_flux_density_by_instrument.pdf"),
       p_flux_inst, width = 9, height = 5)
message("Saved flux density comparison to: ", fig_dir)

message("\n=== Instrument comparison complete ===")
