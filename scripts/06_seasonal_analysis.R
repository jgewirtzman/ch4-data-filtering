# ============================================================================
# 06_seasonal_analysis.R — Seasonal / monthly breakdown of CH4 fluxes
#
# Disentangles biological signal (dormant-season uptake) from instrument
# noise artifacts by comparing LGR (all seasons) vs LI-7810 (growing only).
#
# Analyses:
#   1. Monthly negative flux table by instrument
#   2. Seasonal summary (growing vs dormant)
#   3. Overlapping months comparison (Jun-Oct) between LGR and LI-7810
#   4. LGR dormant vs growing season breakdown
#   5. Negative flux magnitude vs noise floor analysis (% within 1/2/3 sigma)
#
# Plots:
#   1. Monthly histograms faceted by month
#   2. Bar chart of % negative by month and instrument
#   3. Scatter of flux vs noise floor with noise wedge
#   4. Monthly time series with noise envelope ribbons
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
# PART 11: SEASONAL / MONTHLY BREAKDOWN
# ============================================================
#
# LGR covered dormant + growing season (Jun 2023-Dec 2024)
# LI-7810 covered growing season only (Apr-Oct 2025)
# Need to disentangle: which negative fluxes are real biology
# (uptake in dormant season) vs noise artifacts from noisier LGR?
# ============================================================

message("\n=== Seasonal / monthly breakdown ===")

df <- df %>%
  mutate(
    date_parsed = as.Date(date, format = "%Y-%m-%d"),
    month       = as.integer(format(date_parsed, "%m")),
    month_name  = factor(format(date_parsed, "%b"),
                          levels = c("Jan","Feb","Mar","Apr","May","Jun",
                                     "Jul","Aug","Sep","Oct","Nov","Dec")),
    year_month  = format(date_parsed, "%Y-%m"),
    # Season labels
    season = case_when(
      month %in% c(12, 1, 2) ~ "Winter (Dec-Feb)",
      month %in% c(3, 4, 5)  ~ "Spring (Mar-May)",
      month %in% c(6, 7, 8)  ~ "Summer (Jun-Aug)",
      month %in% c(9, 10, 11) ~ "Fall (Sep-Nov)"
    ),
    season = factor(season, levels = c("Spring (Mar-May)", "Summer (Jun-Aug)",
                                        "Fall (Sep-Nov)", "Winter (Dec-Feb)"))
  )

# --- TABLE 1: Monthly breakdown by instrument ---
cat("\n=== Monthly negative flux breakdown by instrument ===\n")
monthly_inst <- df %>%
  group_by(month_name, inst_label) %>%
  summarise(
    n        = n(),
    n_neg    = sum(CH4_flux_nmolpm2ps < 0),
    pct_neg  = round(100 * n_neg / n, 1),
    mean_flux  = round(mean(CH4_flux_nmolpm2ps), 4),
    median_flux = round(median(CH4_flux_nmolpm2ps), 4),
    .groups  = "drop"
  ) %>%
  arrange(month_name, inst_label)
print(as.data.frame(monthly_inst), row.names = FALSE)

# --- TABLE 2: Seasonal summary ---
cat("\n=== Seasonal summary by instrument ===\n")
seasonal_inst <- df %>%
  group_by(season, inst_label) %>%
  summarise(
    n          = n(),
    n_neg      = sum(CH4_flux_nmolpm2ps < 0),
    pct_neg    = round(100 * n_neg / n, 1),
    mean_flux  = round(mean(CH4_flux_nmolpm2ps), 4),
    median_flux = round(median(CH4_flux_nmolpm2ps), 4),
    .groups    = "drop"
  )
print(as.data.frame(seasonal_inst), row.names = FALSE)

# --- TABLE 3: Overlapping months only (Jun-Oct, when both instruments measured) ---
# This is the fairest apples-to-apples comparison
overlap_months <- c(6, 7, 8, 9, 10)  # months present in both instrument periods
cat("\n=== Overlapping months (Jun-Oct) — apples-to-apples comparison ===\n")
df_overlap <- df %>% filter(month %in% overlap_months)

overlap_summary <- df_overlap %>%
  group_by(inst_label) %>%
  summarise(
    n          = n(),
    n_neg      = sum(CH4_flux_nmolpm2ps < 0),
    pct_neg    = round(100 * n_neg / n, 1),
    mean_flux  = round(mean(CH4_flux_nmolpm2ps), 4),
    median_flux = round(median(CH4_flux_nmolpm2ps), 4),
    sd_flux    = round(sd(CH4_flux_nmolpm2ps), 4),
    min_flux   = round(min(CH4_flux_nmolpm2ps), 4),
    max_flux   = round(max(CH4_flux_nmolpm2ps), 4),
    .groups    = "drop"
  )
print(as.data.frame(overlap_summary), row.names = FALSE)

# How many LGR negatives are in dormant vs growing months?
cat("\n=== LGR negative fluxes: dormant vs growing season ===\n")
lgr_only <- df %>% filter(year < 2025)
dormant_months  <- c(11, 12, 1, 2, 3, 4)
growing_months  <- c(5, 6, 7, 8, 9, 10)
lgr_season <- lgr_only %>%
  mutate(period = if_else(month %in% dormant_months,
                           "Dormant (Nov-Apr)", "Growing (May-Oct)")) %>%
  group_by(period) %>%
  summarise(
    n          = n(),
    n_neg      = sum(CH4_flux_nmolpm2ps < 0),
    pct_neg    = round(100 * n_neg / n, 1),
    mean_flux  = round(mean(CH4_flux_nmolpm2ps), 4),
    median_flux = round(median(CH4_flux_nmolpm2ps), 4),
    .groups    = "drop"
  )
print(as.data.frame(lgr_season), row.names = FALSE)

# --- TABLE 4: Negative fluxes vs noise floor ---
# For measurements with Allan deviation: are negatives within the noise?
cat("\n=== Negative flux magnitude vs noise floor ===\n")
df_neg <- df %>%
  filter(CH4_flux_nmolpm2ps < 0, !is.na(CH4_noise_floor))
cat("Negative fluxes with Allan-derived noise floor:", nrow(df_neg), "\n")
if (nrow(df_neg) > 0) {
  df_neg <- df_neg %>%
    mutate(
      abs_flux    = abs(CH4_flux_nmolpm2ps),
      within_1sd  = abs_flux < CH4_noise_floor,
      within_2sd  = abs_flux < 2 * CH4_noise_floor,
      within_3sd  = abs_flux < 3 * CH4_noise_floor
    )
  noise_check <- df_neg %>%
    group_by(inst_label) %>%
    summarise(
      n_neg        = n(),
      within_1_noise = sum(within_1sd),
      within_2_noise = sum(within_2sd),
      within_3_noise = sum(within_3sd),
      pct_within_1   = round(100 * within_1_noise / n_neg, 1),
      pct_within_2   = round(100 * within_2_noise / n_neg, 1),
      pct_within_3   = round(100 * within_3_noise / n_neg, 1),
      median_abs_flux = round(median(abs_flux), 4),
      median_noise    = round(median(CH4_noise_floor), 4),
      .groups = "drop"
    )
  print(as.data.frame(noise_check), row.names = FALSE)
}

# --- PLOT 1: Monthly flux distributions by instrument (faceted) ---

p_monthly <- ggplot(df, aes(x = CH4_flux_nmolpm2ps, fill = inst_label)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 50, alpha = 0.6, position = "identity") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red", linewidth = 0.4) +
  facet_wrap(~ month_name, ncol = 4, scales = "free_y") +
  scale_x_continuous(
    trans = asinh_trans,
    breaks = c(-2, -1, 0, 0.5, 1, 2, 5, 10),
    labels = c("-2", "-1", "0", "0.5", "1", "2", "5", "10")
  ) +
  scale_fill_manual(
    values = c("LGR/UGGA (2023-24)" = "#E69F00", "LI-7810 (2025)" = "#56B4E9"),
    name = "Instrument"
  ) +
  labs(
    x = expression(CH[4]~flux~(nmol~m^{-2}~s^{-1})),
    y = "Density",
    title = expression(Monthly~CH[4]~flux~distributions~by~instrument),
    subtitle = "LGR covers dormant + growing season; LI-7810 covers growing season only (Apr-Oct 2025)"
  ) +
  theme_classic(base_size = 10) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 9, color = "grey40"),
    strip.text = element_text(size = 9, face = "bold")
  ) +
  coord_cartesian(xlim = c(asinh(-3), asinh(15)))

ggsave(file.path(fig_dir, "ch4_flux_monthly_by_instrument.png"),
       p_monthly, width = 12, height = 8, dpi = 300)
ggsave(file.path(fig_dir, "ch4_flux_monthly_by_instrument.pdf"),
       p_monthly, width = 12, height = 8)
message("Saved monthly breakdown to: ", fig_dir)

# --- PLOT 2: Negative flux fraction by month + instrument ---

monthly_neg_rate <- df %>%
  group_by(month_name, inst_label) %>%
  summarise(
    n       = n(),
    pct_neg = 100 * mean(CH4_flux_nmolpm2ps < 0),
    .groups = "drop"
  )

p_neg_month <- ggplot(monthly_neg_rate,
                       aes(x = month_name, y = pct_neg,
                           fill = inst_label, group = inst_label)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.8) +
  geom_text(aes(label = paste0(n)),
            position = position_dodge(width = 0.7),
            vjust = -0.4, size = 2.5, color = "grey40") +
  scale_fill_manual(
    values = c("LGR/UGGA (2023-24)" = "#E69F00", "LI-7810 (2025)" = "#56B4E9"),
    name = "Instrument"
  ) +
  labs(
    x = "Month",
    y = "% negative fluxes",
    title = expression(Fraction~of~negative~CH[4]~fluxes~by~month~and~instrument),
    subtitle = "Numbers above bars = sample size"
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 9, color = "grey40")
  ) +
  coord_cartesian(ylim = c(0, 55))

ggsave(file.path(fig_dir, "ch4_pct_negative_by_month_instrument.png"),
       p_neg_month, width = 9, height = 5, dpi = 300)
ggsave(file.path(fig_dir, "ch4_pct_negative_by_month_instrument.pdf"),
       p_neg_month, width = 9, height = 5)
message("Saved negative % by month to: ", fig_dir)

# --- PLOT 3: Flux vs noise floor scatter, colored by sign ---
# Shows whether negative fluxes sit within the instrument noise

if (sum(!is.na(df$CH4_noise_floor)) > 50) {
  p_noise_scatter <- ggplot(
    df %>% filter(!is.na(CH4_noise_floor)),
    aes(x = CH4_noise_floor, y = CH4_flux_nmolpm2ps, color = inst_label)
  ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.4) +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey50") +
    geom_abline(slope = -1, intercept = 0, linetype = "dotted", color = "grey50") +
    geom_point(alpha = 0.25, size = 1) +
    scale_x_log10(labels = scales::label_number()) +
    scale_y_continuous(
      trans = asinh_trans,
      breaks = c(-5, -2, -1, 0, 0.5, 1, 2, 5, 10, 20),
      labels = c("-5", "-2", "-1", "0", "0.5", "1", "2", "5", "10", "20")
    ) +
    scale_color_manual(
      values = c("LGR/UGGA (2023-24)" = "#E69F00", "LI-7810 (2025)" = "#56B4E9"),
      name = "Instrument"
    ) +
    annotate("text", x = 0.3, y = 0.4, label = "flux = +noise floor",
             angle = 18, size = 3, color = "grey50") +
    annotate("text", x = 0.3, y = -0.4, label = "flux = \u2013noise floor",
             angle = -18, size = 3, color = "grey50") +
    labs(
      x = expression(CH[4]~noise~floor~(nmol~m^{-2}~s^{-1})),
      y = expression(CH[4]~flux~(nmol~m^{-2}~s^{-1})),
      title = expression(CH[4]~flux~vs.~per-measurement~noise~floor),
      subtitle = "Dotted lines = \u00b1 noise floor (1\u03c3 Allan). Points within dotted wedge are indistinguishable from zero."
    ) +
    theme_classic(base_size = 11) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(size = 13, face = "bold"),
      plot.subtitle = element_text(size = 9, color = "grey40")
    )

  ggsave(file.path(fig_dir, "ch4_flux_vs_noise_floor.png"),
         p_noise_scatter, width = 8, height = 6, dpi = 300)
  ggsave(file.path(fig_dir, "ch4_flux_vs_noise_floor.pdf"),
         p_noise_scatter, width = 8, height = 6)
  message("Saved flux vs noise floor scatter to: ", fig_dir)
}

# --- PLOT 4: Monthly time series with noise envelope ---
# Shows median flux +/- noise floor by month and instrument

monthly_ts <- df %>%
  filter(!is.na(CH4_noise_floor)) %>%
  group_by(year_month, inst_label) %>%
  summarise(
    date_mid    = mean(date_parsed),
    median_flux = median(CH4_flux_nmolpm2ps),
    q25_flux    = quantile(CH4_flux_nmolpm2ps, 0.25),
    q75_flux    = quantile(CH4_flux_nmolpm2ps, 0.75),
    median_noise = median(CH4_noise_floor),
    n            = n(),
    .groups = "drop"
  )

p_ts <- ggplot(monthly_ts, aes(x = date_mid, color = inst_label, fill = inst_label)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.3) +
  # Noise envelope around zero
  geom_ribbon(aes(ymin = -median_noise, ymax = median_noise), alpha = 0.15,
              color = NA) +
  # IQR ribbon
  geom_ribbon(aes(ymin = q25_flux, ymax = q75_flux), alpha = 0.25, color = NA) +
  # Median line
  geom_line(aes(y = median_flux), linewidth = 0.8) +
  geom_point(aes(y = median_flux, size = n), alpha = 0.7) +
  scale_size_continuous(range = c(1.5, 4), name = "n measurements") +
  scale_color_manual(
    values = c("LGR/UGGA (2023-24)" = "#E69F00", "LI-7810 (2025)" = "#56B4E9"),
    name = "Instrument"
  ) +
  scale_fill_manual(
    values = c("LGR/UGGA (2023-24)" = "#E69F00", "LI-7810 (2025)" = "#56B4E9"),
    name = "Instrument"
  ) +
  scale_y_continuous(
    trans = asinh_trans,
    breaks = c(-2, -1, -0.5, 0, 0.5, 1, 2, 5, 10),
    labels = c("-2", "-1", "-0.5", "0", "0.5", "1", "2", "5", "10")
  ) +
  labs(
    x = "Date",
    y = expression(CH[4]~flux~(nmol~m^{-2}~s^{-1})),
    title = expression(Monthly~median~CH[4]~flux~with~noise~envelope),
    subtitle = "Line = median flux | Dark ribbon = IQR | Light ribbon = \u00b1 median noise floor (1\u03c3 Allan)"
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 9, color = "grey40")
  )

ggsave(file.path(fig_dir, "ch4_flux_monthly_timeseries_noise.png"),
       p_ts, width = 11, height = 5, dpi = 300)
ggsave(file.path(fig_dir, "ch4_flux_monthly_timeseries_noise.pdf"),
       p_ts, width = 11, height = 5)
message("Saved monthly time series with noise to: ", fig_dir)

message("\n=== Seasonal analysis complete ===")
