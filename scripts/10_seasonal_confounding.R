# ============================================================================
# 10_seasonal_confounding.R — Seasonal structure in quality filter exclusions
#
# Demonstrates that quality filters (R², MDF, SNR) do not exclude measurements
# randomly but preferentially exclude cold-season, low-flux conditions ---
# introducing structured, ecologically correlated bias.
#
# Figures:
#   A. Seasonal fraction flagged by month per filter criterion
#   B. Flux magnitude vs R² colored by air temperature
#   C. Annual time-series colored by filter status
#   D. Flux vs air temperature colored by filter status
#   E. Cumulative annual budget under different filter approaches
#
# Loads: outputs/03_mdf_results.RData
# Saves: figures to fig_dir (PNG + PDF)
# ============================================================================

source("scripts/00_setup.R")

message("\n=== Loading MDF results ===")
load(file.path(output_dir, "03_mdf_results.RData"))

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# --- Prepare date columns ---
df$date_parsed <- as.Date(df$date)
df$month_num   <- as.numeric(format(df$date_parsed, "%m"))
df$month_lab   <- factor(month.abb[df$month_num], levels = month.abb)
df$year_month  <- format(df$date_parsed, "%Y-%m")

# Instrument label
df$inst_label <- ifelse(df$year == 2025, "Analyzer B",
                         "Analyzer A")

# Season
df$season <- ifelse(df$month_num %in% 5:10, "Growing (May\u2013Oct)",
                     "Dormant (Nov\u2013Apr)")

# Four-season classification
df$season4 <- factor(
  case_when(
    df$month_num %in% c(12, 1, 2) ~ "Winter (DJF)",
    df$month_num %in% 3:5          ~ "Spring (MAM)",
    df$month_num %in% 6:8          ~ "Summer (JJA)",
    df$month_num %in% 9:11         ~ "Fall (SON)"
  ),
  levels = c("Winter (DJF)", "Spring (MAM)", "Summer (JJA)", "Fall (SON)")
)


# ============================================================================
# FIGURE A: SEASONAL FRACTION FLAGGED BY SEASON PER FILTER CRITERION
# ============================================================================

message("\n=== Figure A: Seasonal fraction flagged ===")

# Select key filter criteria to show
key_filters <- list(
  "CH4 R\u00b2 > 0.5"     = function(d) d$CH4_r2 <= 0.5,
  "CH4 R\u00b2 > 0.7"     = function(d) d$CH4_r2 <= 0.7,
  "Manufacturer MDF"       = function(d) ifelse(!is.na(d$CH4_below_MDF_manuf),
                                                  d$CH4_below_MDF_manuf, FALSE),
  "Wassmann 95%"           = function(d) ifelse(!is.na(d$CH4_below_MDF_wass95),
                                                  d$CH4_below_MDF_wass95, FALSE),
  "Christiansen 95%"       = function(d) ifelse(!is.na(d$CH4_below_MDF_chr95),
                                                  d$CH4_below_MDF_chr95, FALSE),
  "Allan SNR > 2"          = function(d) ifelse(!is.na(d$CH4_snr_allan),
                                                  d$CH4_snr_allan <= 2, FALSE)
)

# Compute per-SEASON fraction flagged for each criterion
seasonal_flagged <- lapply(names(key_filters), function(fn) {
  fails <- key_filters[[fn]](df)
  out <- df %>%
    mutate(fails_filter = fails) %>%
    group_by(season4) %>%
    summarise(
      n_total     = n(),
      n_flagged   = sum(fails_filter, na.rm = TRUE),
      pct_flagged = 100 * n_flagged / n_total,
      mean_airt   = mean(airt, na.rm = TRUE),
      .groups = "drop"
    )
  out$filter <- fn
  out
})
seasonal_flagged_df <- do.call(rbind, seasonal_flagged)

# Order filters for legend
seasonal_flagged_df$filter <- factor(seasonal_flagged_df$filter,
                                      levels = names(key_filters))

# Keep dormant_rects for cumulative budget figure later
dormant_rects <- data.frame(
  xmin = c(0.5, 10.5),
  xmax = c(4.5, 12.5),
  ymin = -Inf, ymax = Inf
)

p_seasonal_frac <- ggplot(seasonal_flagged_df,
                           aes(x = season4, y = pct_flagged,
                               color = filter, group = filter)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  scale_color_brewer(palette = "Dark2", name = "Filter criterion") +
  labs(
    x = NULL,
    y = "% of measurements excluded",
    title = "Quality filter exclusion rates by season",
    subtitle = paste0("Stem flux dataset (n=", nrow(df),
                      "). All filter criteria show elevated exclusion\n",
                      "in cold seasons (winter/spring) when true fluxes approach zero.")
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 9),
    plot.title = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 9, color = "grey40")
  ) +
  guides(color = guide_legend(nrow = 2))

ggsave(file.path(fig_dir, "seasonal_fraction_flagged.png"),
       p_seasonal_frac, width = 10, height = 6, dpi = 300)
ggsave(file.path(fig_dir, "seasonal_fraction_flagged.pdf"),
       p_seasonal_frac, width = 10, height = 6)
message("Saved: seasonal_fraction_flagged")

# Print summary: seasonal exclusion rates
cat("\n=== Exclusion rates by season ===\n")
season_summary <- seasonal_flagged_df %>%
  group_by(filter, season4) %>%
  summarise(
    pct_flagged = round(mean(pct_flagged), 1),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = season4, values_from = pct_flagged)
print(as.data.frame(season_summary), row.names = FALSE)


# ============================================================================
# FIGURE A2: FRACTION FLAGGED BY 5°C TEMPERATURE BIN
# ============================================================================

message("\n=== Figure A2: Fraction flagged by temperature bin ===")

# Create 10°C bins
df$temp_bin <- cut(df$airt,
                   breaks = seq(-10, 30, by = 10),
                   include.lowest = TRUE,
                   right = FALSE,
                   labels = paste0(seq(-10, 20, by = 10), " to ", seq(0, 30, by = 10)))

# Drop bins with no data
df$temp_bin <- droplevels(df$temp_bin)

# Compute per-temp-bin fraction flagged for each criterion
temp_flagged <- lapply(names(key_filters), function(fn) {
  fails <- key_filters[[fn]](df)
  out <- df %>%
    mutate(fails_filter = fails) %>%
    filter(!is.na(temp_bin)) %>%
    group_by(temp_bin) %>%
    summarise(
      n_total     = n(),
      n_flagged   = sum(fails_filter, na.rm = TRUE),
      pct_flagged = 100 * n_flagged / n_total,
      .groups = "drop"
    )
  out$filter <- fn
  out
})
temp_flagged_df <- do.call(rbind, temp_flagged)
temp_flagged_df$filter <- factor(temp_flagged_df$filter, levels = names(key_filters))

p_temp_frac <- ggplot(temp_flagged_df,
                       aes(x = temp_bin, y = pct_flagged,
                           color = filter, group = filter)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  scale_color_brewer(palette = "Dark2", name = "Filter criterion") +
  labs(
    x = "Air temperature bin (C)",
    y = "% of measurements excluded",
    title = "Quality filter exclusion rates by temperature",
    subtitle = paste0("Stem flux dataset (n=", nrow(df),
                      "). Exclusion rates decrease monotonically\n",
                      "with temperature, confirming temperature-driven bias.")
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 9),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    plot.title = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 9, color = "grey40")
  ) +
  guides(color = guide_legend(nrow = 2))

ggsave(file.path(fig_dir, "temp_bin_fraction_flagged.png"),
       p_temp_frac, width = 10, height = 6, dpi = 300)
ggsave(file.path(fig_dir, "temp_bin_fraction_flagged.pdf"),
       p_temp_frac, width = 10, height = 6)
message("Saved: temp_bin_fraction_flagged")

cat("\n=== Exclusion rates by temperature bin ===\n")
temp_summary <- temp_flagged_df %>%
  group_by(filter, temp_bin) %>%
  summarise(pct = round(mean(pct_flagged), 1), .groups = "drop") %>%
  pivot_wider(names_from = temp_bin, values_from = pct)
print(as.data.frame(temp_summary), row.names = FALSE)


# ============================================================================
# FIGURE B: FLUX MAGNITUDE VS R² COLORED BY TEMPERATURE
# ============================================================================

message("\n=== Figure B: |Flux| vs R\u00b2 by temperature ===")

# Use Analyzer A only (Analyzer B covers only May-Oct, biasing the seasonal pattern)
df_gla <- df %>% filter(inst_label != "Analyzer B")
message("  Analyzer A measurements: ", nrow(df_gla))

# Compute mean temperature per R\u00b2 bin x sign for annotations
r2_sign_annot <- df_gla %>%
  mutate(
    r2_bin = case_when(
      CH4_r2 > 0.7 ~ "R2 > 0.7",
      CH4_r2 > 0.5 ~ "0.5 < R2 <= 0.7",
      TRUE          ~ "R2 <= 0.5"
    ),
    sign = ifelse(CH4_flux_nmolpm2ps >= 0, "positive", "negative")
  ) %>%
  group_by(r2_bin, sign) %>%
  summarise(
    mean_temp = round(mean(airt, na.rm = TRUE), 1),
    n = n(),
    .groups = "drop"
  )

cat("\nMean air temp by R2 bin x sign:\n")
print(as.data.frame(r2_sign_annot), row.names = FALSE)

# Build annotation table for the plot
annot_tbl <- r2_sign_annot %>%
  mutate(
    label = paste0(r2_bin, " (", sign, "): T=", mean_temp, "C, n=", n)
  )

# Position annotations: right side for positive flux, left side for negative
annot_pos <- annot_tbl %>% filter(sign == "positive")
annot_neg <- annot_tbl %>% filter(sign == "negative")

p_r2_flux <- ggplot(df_gla, aes(x = CH4_flux_nmolpm2ps, y = CH4_r2)) +
  geom_point(aes(color = airt), alpha = 0.5, size = 1.2) +
  geom_hline(yintercept = 0.5, linetype = "dashed",
             color = "red", linewidth = 0.4) +
  geom_hline(yintercept = 0.7, linetype = "dashed",
             color = "darkred", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey50", linewidth = 0.3) +
  # Annotate: positive-flux side (right)
  annotate("label", x = 10, y = 0.87,
           label = paste0("R2>0.7, (+): T=", annot_pos$mean_temp[annot_pos$r2_bin=="R2 > 0.7"],
                          "C  n=", annot_pos$n[annot_pos$r2_bin=="R2 > 0.7"]),
           size = 2.5, fill = "white", alpha = 0.8, hjust = 0, color = "darkred") +
  annotate("label", x = 10, y = 0.6,
           label = paste0("0.5<R2<0.7, (+): T=", annot_pos$mean_temp[annot_pos$r2_bin=="0.5 < R2 <= 0.7"],
                          "C  n=", annot_pos$n[annot_pos$r2_bin=="0.5 < R2 <= 0.7"]),
           size = 2.5, fill = "white", alpha = 0.8, hjust = 0, color = "red") +
  annotate("label", x = 10, y = 0.25,
           label = paste0("R2<0.5, (+): T=", annot_pos$mean_temp[annot_pos$r2_bin=="R2 <= 0.5"],
                          "C  n=", annot_pos$n[annot_pos$r2_bin=="R2 <= 0.5"]),
           size = 2.5, fill = "white", alpha = 0.8, hjust = 0, color = "#B2182B") +
  # Annotate: negative-flux side (left)
  annotate("label", x = -10, y = 0.87,
           label = paste0("R2>0.7, (-): T=",
                          ifelse(any(annot_neg$r2_bin=="R2 > 0.7"),
                                 annot_neg$mean_temp[annot_neg$r2_bin=="R2 > 0.7"], "NA"),
                          "C  n=",
                          ifelse(any(annot_neg$r2_bin=="R2 > 0.7"),
                                 annot_neg$n[annot_neg$r2_bin=="R2 > 0.7"], 0)),
           size = 2.5, fill = "white", alpha = 0.8, hjust = 1, color = "darkred") +
  annotate("label", x = -10, y = 0.6,
           label = paste0("0.5<R2<0.7, (-): T=",
                          ifelse(any(annot_neg$r2_bin=="0.5 < R2 <= 0.7"),
                                 annot_neg$mean_temp[annot_neg$r2_bin=="0.5 < R2 <= 0.7"], "NA"),
                          "C  n=",
                          ifelse(any(annot_neg$r2_bin=="0.5 < R2 <= 0.7"),
                                 annot_neg$n[annot_neg$r2_bin=="0.5 < R2 <= 0.7"], 0)),
           size = 2.5, fill = "white", alpha = 0.8, hjust = 1, color = "red") +
  annotate("label", x = -10, y = 0.25,
           label = paste0("R2<0.5, (-): T=", annot_neg$mean_temp[annot_neg$r2_bin=="R2 <= 0.5"],
                          "C  n=", annot_neg$n[annot_neg$r2_bin=="R2 <= 0.5"]),
           size = 2.5, fill = "white", alpha = 0.8, hjust = 1, color = "#B2182B") +
  scale_x_continuous(
    trans = scales::pseudo_log_trans(sigma = 0.01),
    limits = c(-50, 50),
    breaks = c(-50, -10, -1, -0.1, 0, 0.1, 1, 10, 50),
    labels = c("-50", "-10", "-1", "-0.1", "0", "0.1", "1", "10", "50")
  ) +
  scale_color_viridis_c(option = "plasma", name = "Air temp (C)") +
  labs(
    x = expression(CH[4]~flux~(nmol~m^{-2}~s^{-1})),
    y = expression(CH[4]~R^2),
    title = expression(CH[4]~R^2~vs~signed~flux~and~temperature),
    subtitle = paste0("Analyzer A stem data (n=", nrow(df_gla), "). ",
                      "Negative fluxes cluster below R2 thresholds and at cold temperatures.")
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 9, color = "grey40")
  )

ggsave(file.path(fig_dir, "r2_vs_flux_by_temperature.png"),
       p_r2_flux, width = 11, height = 6, dpi = 300)
ggsave(file.path(fig_dir, "r2_vs_flux_by_temperature.pdf"),
       p_r2_flux, width = 11, height = 6)
message("Saved: r2_vs_flux_by_temperature")

# --- SIGNED flux vs R² ---
message("\n=== Figure B2: Signed flux vs R\u00b2 by temperature ===")

p_r2_flux_signed <- ggplot(df, aes(x = CH4_flux_nmolpm2ps, y = CH4_r2)) +
  geom_point(aes(color = airt), alpha = 0.5, size = 1.2) +
  geom_hline(yintercept = 0.5, linetype = "dashed",
             color = "red", linewidth = 0.4) +
  geom_hline(yintercept = 0.7, linetype = "dashed",
             color = "darkred", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey50", linewidth = 0.3) +
  scale_x_continuous(
    trans = scales::pseudo_log_trans(sigma = 0.01),
    breaks = c(-1, -0.1, 0, 0.1, 1, 10, 50),
    labels = c("-1", "-0.1", "0", "0.1", "1", "10", "50")
  ) +
  scale_color_viridis_c(option = "plasma", name = "Air temp (\u00b0C)") +
  facet_wrap(~ inst_label, ncol = 2) +
  labs(
    x = expression(CH[4]~flux~(nmol~m^{-2}~s^{-1})),
    y = expression(CH[4]~R^2),
    title = expression(CH[4]~R^2~vs~signed~flux~and~temperature),
    subtitle = paste0("Negative fluxes cluster below R\u00b2 thresholds and at cold temperatures.\n",
                      "R\u00b2 filters remove near-zero and negative values non-randomly.")
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "right",
    strip.text = element_text(size = 11, face = "bold"),
    plot.title = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 9, color = "grey40")
  )

ggsave(file.path(fig_dir, "r2_vs_flux_signed_by_temperature.png"),
       p_r2_flux_signed, width = 12, height = 6, dpi = 300)
ggsave(file.path(fig_dir, "r2_vs_flux_signed_by_temperature.pdf"),
       p_r2_flux_signed, width = 12, height = 6)
message("Saved: r2_vs_flux_signed_by_temperature")

# --- SIGNED flux vs R², LINEAR x-axis ---
message("\n=== Figure B3: Signed flux vs R\u00b2 by temperature (linear scale) ===")

p_r2_flux_linear <- ggplot(df, aes(x = CH4_flux_nmolpm2ps, y = CH4_r2)) +
  geom_point(aes(color = airt), alpha = 0.5, size = 1.2) +
  geom_hline(yintercept = 0.5, linetype = "dashed",
             color = "red", linewidth = 0.4) +
  geom_hline(yintercept = 0.7, linetype = "dashed",
             color = "darkred", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey50", linewidth = 0.3) +
  scale_color_viridis_c(option = "plasma", name = "Air temp (\u00b0C)") +
  facet_wrap(~ inst_label, ncol = 2) +
  labs(
    x = expression(CH[4]~flux~(nmol~m^{-2}~s^{-1})),
    y = expression(CH[4]~R^2),
    title = expression(CH[4]~R^2~vs~signed~flux~and~temperature~(linear~scale)),
    subtitle = paste0("Linear x-axis emphasises the skewed flux distribution. ",
                      "Low-R\u00b2 measurements\n",
                      "are compressed near zero while high-flux outliers dominate the axis range.")
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "right",
    strip.text = element_text(size = 11, face = "bold"),
    plot.title = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 9, color = "grey40")
  )

ggsave(file.path(fig_dir, "r2_vs_flux_linear_by_temperature.png"),
       p_r2_flux_linear, width = 12, height = 6, dpi = 300)
ggsave(file.path(fig_dir, "r2_vs_flux_linear_by_temperature.pdf"),
       p_r2_flux_linear, width = 12, height = 6)
message("Saved: r2_vs_flux_linear_by_temperature")

# Print correlation stats
cat("\n=== R\u00b2 vs |flux| and temperature ===\n")
cat("Correlation |flux| vs R\u00b2:",
    round(cor(abs(df$CH4_flux_nmolpm2ps), df$CH4_r2, use = "complete"), 3), "\n")
cat("Correlation airt vs R\u00b2:",
    round(cor(df$airt, df$CH4_r2, use = "complete"), 3), "\n")

# R² by season
r2_season <- df %>%
  group_by(season) %>%
  summarise(
    median_r2 = round(median(CH4_r2), 3),
    pct_below_0.5 = round(100 * mean(CH4_r2 < 0.5), 1),
    pct_below_0.7 = round(100 * mean(CH4_r2 < 0.7), 1),
    mean_airt = round(mean(airt, na.rm = TRUE), 1),
    .groups = "drop"
  )
cat("\nR\u00b2 by season:\n")
print(as.data.frame(r2_season), row.names = FALSE)


# ============================================================================
# FIGURE C: ANNUAL TIME-SERIES COLORED BY FILTER STATUS
# ============================================================================

message("\n=== Figure C: Time-series by filter status ===")

# Classify each measurement by filter status
# Use two key filters: R² > 0.5 and Christiansen 95% MDF
df$fails_r2 <- df$CH4_r2 <= 0.5
df$fails_mdf <- ifelse(!is.na(df$CH4_below_MDF_chr95),
                        df$CH4_below_MDF_chr95, FALSE)

df$filter_status <- case_when(
  df$fails_r2 & df$fails_mdf  ~ "Fails both",
  df$fails_r2 & !df$fails_mdf ~ "Fails R\u00b2 only",
  !df$fails_r2 & df$fails_mdf ~ "Fails MDF only",
  TRUE                         ~ "Passes all"
)

df$filter_status <- factor(df$filter_status,
                            levels = c("Passes all", "Fails MDF only",
                                       "Fails R\u00b2 only", "Fails both"))

status_colors <- c("Passes all"      = "#2166AC",
                    "Fails MDF only"  = "#F4A582",
                    "Fails R\u00b2 only"  = "#D6604D",
                    "Fails both"      = "#B2182B")

# Monthly temperature ribbon for context
monthly_temp <- df %>%
  group_by(year_month) %>%
  summarise(
    date_mid = median(date_parsed),
    mean_airt = mean(airt, na.rm = TRUE),
    .groups = "drop"
  )

p_timeseries <- ggplot(df, aes(x = date_parsed, y = CH4_flux_nmolpm2ps)) +
  # Temperature ribbon on secondary axis (scaled to fit)
  geom_point(aes(color = filter_status), alpha = 0.5, size = 1.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  scale_color_manual(values = status_colors, name = "Filter status") +
  scale_x_date(date_breaks = "3 months", date_labels = "%b\n%Y") +
  scale_y_continuous(
    trans = scales::pseudo_log_trans(sigma = 0.01),
    breaks = c(-1, -0.1, 0, 0.1, 1, 10, 50),
    labels = c("-1", "-0.1", "0", "0.1", "1", "10", "50")
  ) +
  labs(
    x = NULL,
    y = expression(CH[4]~flux~(nmol~m^{-2}~s^{-1})),
    title = "Temporal distribution of measurements failing quality filters",
    subtitle = paste0("Filters: CH4 R\u00b2 > 0.5 and Christiansen MDF 95%. ",
                      "Excluded measurements cluster in dormant months,\n",
                      "not randomly distributed in time.")
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 9, color = "grey40")
  )

ggsave(file.path(fig_dir, "timeseries_filter_status.png"),
       p_timeseries, width = 12, height = 6, dpi = 300)
ggsave(file.path(fig_dir, "timeseries_filter_status.pdf"),
       p_timeseries, width = 12, height = 6)
message("Saved: timeseries_filter_status")

# Print filter status breakdown by season
cat("\n=== Filter status by season ===\n")
status_season <- df %>%
  group_by(season, filter_status) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(season) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  ungroup()
print(as.data.frame(status_season), row.names = FALSE)


# ============================================================================
# FIGURE D: FLUX VS AIR TEMPERATURE COLORED BY FILTER STATUS
# ============================================================================

message("\n=== Figure D: Flux vs temperature by filter status ===")

p_flux_temp <- ggplot(df, aes(x = airt, y = CH4_flux_nmolpm2ps)) +
  geom_point(aes(color = filter_status), alpha = 0.5, size = 1.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey70", linewidth = 0.3) +
  scale_color_manual(values = status_colors, name = "Filter status") +
  scale_y_continuous(
    trans = scales::pseudo_log_trans(sigma = 0.01),
    breaks = c(-1, -0.1, 0, 0.1, 1, 10, 50),
    labels = c("-1", "-0.1", "0", "0.1", "1", "10", "50")
  ) +
  facet_wrap(~ inst_label, ncol = 2) +
  labs(
    x = "Air temperature (\u00b0C)",
    y = expression(CH[4]~flux~(nmol~m^{-2}~s^{-1})),
    title = "Quality filter exclusions are structured by temperature",
    subtitle = paste0("Excluded measurements (red tones) cluster at low temperatures.\n",
                      "Filters: CH4 R\u00b2 > 0.5 and Christiansen MDF 95%.")
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 11, face = "bold"),
    plot.title = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 9, color = "grey40")
  )

ggsave(file.path(fig_dir, "flux_vs_temperature_filter_status.png"),
       p_flux_temp, width = 12, height = 6, dpi = 300)
ggsave(file.path(fig_dir, "flux_vs_temperature_filter_status.pdf"),
       p_flux_temp, width = 12, height = 6)
message("Saved: flux_vs_temperature_filter_status")

# Correlation of exclusion with temperature
cat("\n=== Temperature distribution of excluded vs retained ===\n")
temp_by_status <- df %>%
  mutate(excluded = filter_status != "Passes all") %>%
  group_by(excluded) %>%
  summarise(
    n = n(),
    mean_airt = round(mean(airt, na.rm = TRUE), 1),
    median_airt = round(median(airt, na.rm = TRUE), 1),
    .groups = "drop"
  )
print(as.data.frame(temp_by_status), row.names = FALSE)


# ============================================================================
# FIGURE E: CUMULATIVE ANNUAL BUDGET UNDER DIFFERENT FILTER APPROACHES
# ============================================================================

message("\n=== Figure E: Cumulative annual budget ===")

# Use calendar year 2024 (Analyzer A only, all 12 months represented)
df_2024 <- df %>%
  filter(year == 2024) %>%
  arrange(date_parsed)

n_2024 <- nrow(df_2024)
message("  2024 measurements: ", n_2024)
message("  Months covered: ", paste(sort(unique(df_2024$month_num)), collapse = ", "))

# Define filter treatments
treatments <- list(
  "No filter (keep all)"        = function(d) d$CH4_flux_nmolpm2ps,
  "R\u00b2 > 0.5 (remove)"     = function(d) {
    keep <- d$CH4_r2 > 0.5
    d$CH4_flux_nmolpm2ps[keep]
  },
  "Christiansen 95%\n(remove)"  = function(d) {
    keep <- !ifelse(!is.na(d$CH4_below_MDF_chr95), d$CH4_below_MDF_chr95, FALSE)
    d$CH4_flux_nmolpm2ps[keep]
  },
  "Christiansen 95%\n(keep all, flag)" = function(d) d$CH4_flux_nmolpm2ps
)

# Compute monthly mean flux under each treatment
monthly_budget <- lapply(names(treatments), function(trt) {
  vals <- treatments[[trt]](df_2024)

  if (trt %in% c("No filter (keep all)", "Christiansen 95%\n(keep all, flag)")) {
    # All measurements kept — compute monthly means from full data
    monthly <- df_2024 %>%
      group_by(month_num) %>%
      summarise(
        mean_flux = mean(CH4_flux_nmolpm2ps, na.rm = TRUE),
        n = n(),
        .groups = "drop"
      )
  } else {
    # Filtered — need to flag which rows are kept
    if (grepl("R\u00b2", trt)) {
      keep <- df_2024$CH4_r2 > 0.5
    } else {
      keep <- !ifelse(!is.na(df_2024$CH4_below_MDF_chr95),
                       df_2024$CH4_below_MDF_chr95, FALSE)
    }
    monthly <- df_2024 %>%
      mutate(keep = keep) %>%
      filter(keep) %>%
      group_by(month_num) %>%
      summarise(
        mean_flux = mean(CH4_flux_nmolpm2ps, na.rm = TRUE),
        n = n(),
        .groups = "drop"
      )
  }

  # Fill missing months with 0 for cumulative sum (or NA)
  all_months <- data.frame(month_num = 1:12)
  monthly <- merge(all_months, monthly, by = "month_num", all.x = TRUE)
  monthly$mean_flux[is.na(monthly$mean_flux)] <- 0
  monthly$n[is.na(monthly$n)] <- 0

  monthly$treatment <- trt
  monthly$cum_flux  <- cumsum(monthly$mean_flux)
  monthly
})
monthly_budget_df <- do.call(rbind, monthly_budget)

monthly_budget_df$treatment <- factor(monthly_budget_df$treatment,
                                       levels = names(treatments))

# "No filter" and "keep all, flag" should produce identical lines
# (they are the same data, just different conceptual framing)
# So merge them conceptually

treatment_colors <- c(
  "No filter (keep all)"               = "black",
  "R\u00b2 > 0.5 (remove)"            = "#D6604D",
  "Christiansen 95%\n(remove)"         = "#E69F00",
  "Christiansen 95%\n(keep all, flag)" = "#2166AC"
)

p_cumulative <- ggplot(monthly_budget_df,
                        aes(x = month_num, y = cum_flux,
                            color = treatment, linetype = treatment)) +
  geom_rect(data = dormant_rects,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE, fill = "grey90", alpha = 0.5) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +
  scale_color_manual(values = treatment_colors, name = "Treatment") +
  scale_linetype_manual(
    values = c("No filter (keep all)" = "solid",
               "R\u00b2 > 0.5 (remove)" = "dashed",
               "Christiansen 95%\n(remove)" = "dashed",
               "Christiansen 95%\n(keep all, flag)" = "solid"),
    name = "Treatment"
  ) +
  labs(
    x = NULL,
    y = expression(Cumulative~mean~CH[4]~flux~(nmol~m^{-2}~s^{-1})),
    title = "Cumulative annual CH\u2084 flux estimate diverges with filter treatment",
    subtitle = paste0("Analyzer A stem data, calendar year 2024 (n=", n_2024,
                      "). Grey bands = dormant season.\n",
                      "Removal treatments inflate the annual estimate; ",
                      "divergence grows through dormant months.")
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 9),
    plot.title = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 9, color = "grey40")
  ) +
  guides(color = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2))

ggsave(file.path(fig_dir, "cumulative_annual_budget.png"),
       p_cumulative, width = 10, height = 6, dpi = 300)
ggsave(file.path(fig_dir, "cumulative_annual_budget.pdf"),
       p_cumulative, width = 10, height = 6)
message("Saved: cumulative_annual_budget")

# Print final cumulative values
cat("\n=== Final annual cumulative flux (sum of monthly means) ===\n")
final_vals <- monthly_budget_df %>%
  filter(month_num == 12) %>%
  select(treatment, cum_flux)
final_vals$pct_change <- round(100 * (final_vals$cum_flux / final_vals$cum_flux[1] - 1), 1)
print(as.data.frame(final_vals), row.names = FALSE)


# ============================================================================
# SUMMARY STATISTICS
# ============================================================================

cat("\n=== Summary: Seasonal confounding evidence ===\n")

# Mean temperature of excluded vs retained measurements
cat("\nExcluded (by any of R\u00b2>0.5 or Chr95%) vs retained:\n")
excluded <- df$filter_status != "Passes all"
cat("  Retained: n=", sum(!excluded),
    ", mean airt=", round(mean(df$airt[!excluded], na.rm=TRUE), 1), "\u00b0C",
    ", median flux=", round(median(df$CH4_flux_nmolpm2ps[!excluded]), 4), "\n")
cat("  Excluded: n=", sum(excluded),
    ", mean airt=", round(mean(df$airt[excluded], na.rm=TRUE), 1), "\u00b0C",
    ", median flux=", round(median(df$CH4_flux_nmolpm2ps[excluded]), 4), "\n")

# Fraction excluded by season
cat("\nFraction excluded by season:\n")
cat("  Growing season:  ", round(100 * mean(excluded[df$season == "Growing (May\u2013Oct)"]), 1), "%\n")
cat("  Dormant season:  ", round(100 * mean(excluded[df$season == "Dormant (Nov\u2013Apr)"]), 1), "%\n")

# --- Save stem summary stats before soil section overwrites df ---
stem_n <- nrow(df)
stem_pct_dorm_excl <- round(100 * mean(
  excluded[df$season == "Dormant (Nov\u2013Apr)"]), 1)
stem_pct_grow_excl <- round(100 * mean(
  excluded[df$season == "Growing (May\u2013Oct)"]), 1)

message("\n=== Stem seasonal confounding analysis complete ===")


# ############################################################################
# ############################################################################
#
#  SOIL DATA: REPEAT SEASONAL CONFOUNDING ANALYSIS
#
#  The wetland soil dataset (YMF semirigid chamber, n~288) lacks raw
#  concentration timeseries, so MDF/SNR cannot be computed.  However the
#  R²-based confounding analysis translates directly.
#
#  Available columns
#    CH4_LM_r2             → R² of the linear model (used for ALL rows,
#                             including those where HM was selected as best)
#    chamber_temp_C        → chamber headspace temperature (°C)
#    CH4_best_flux_nmol_m2_s → best-fit CH4 flux
#    date                  → measurement date
#    measurement_type      → "soil" or "tree_stem"
#    landscape_position    → U / I / WD / WS
#
# ############################################################################
# ############################################################################

message("\n\n========================================================")
message("=== SOIL DATA: Seasonal confounding analysis ===")
message("========================================================\n")

# --- Load soil data from 01_loaded_data.RData ---
# NOTE: this overwrites df; stem stats saved above
load(file.path(output_dir, "01_loaded_data.RData"))

soil <- soil_flux %>%
  filter(measurement_type == "soil",
         !is.na(CH4_best_flux_nmol_m2_s),
         !is.na(CH4_LM_r2))

n_soil <- nrow(soil)
message("Soil measurements (after NA removal): ", n_soil)

# --- Standardise column names to match stem analysis ---
soil <- soil %>%
  rename(
    flux    = CH4_best_flux_nmol_m2_s,
    r2      = CH4_LM_r2,
    temp_C  = chamber_temp_C
  ) %>%
  mutate(
    date_parsed = as.Date(date),
    month_num   = as.numeric(format(date_parsed, "%m")),
    month_lab   = factor(month.abb[month_num], levels = month.abb),
    year        = as.numeric(format(date_parsed, "%Y")),
    year_month  = format(date_parsed, "%Y-%m"),
    season      = ifelse(month_num %in% 5:10, "Growing (May\u2013Oct)",
                          "Dormant (Nov\u2013Apr)"),
    season4     = factor(
      case_when(
        month_num %in% c(12, 1, 2) ~ "Winter (DJF)",
        month_num %in% 3:5          ~ "Spring (MAM)",
        month_num %in% 6:8          ~ "Summer (JJA)",
        month_num %in% 9:11         ~ "Fall (SON)"
      ),
      levels = c("Winter (DJF)", "Spring (MAM)", "Summer (JJA)", "Fall (SON)")
    )
  )

message("Date range: ", min(soil$date), " to ", max(soil$date))
message("Months:     ", paste(sort(unique(soil$month_num)), collapse = ", "))
message("Temp range: ", round(min(soil$temp_C, na.rm = TRUE), 1),
        " to ", round(max(soil$temp_C, na.rm = TRUE), 1), "\u00b0C")
pct_neg <- round(100 * mean(soil$flux < 0), 1)
message("% negative: ", pct_neg, "%")


# ============================================================================
# SOIL FIGURE A: SEASONAL FRACTION FLAGGED (R²-based filters only)
# ============================================================================

message("\n=== Soil Figure A: Seasonal fraction flagged ===")

soil_filters <- list(
  "CH4 R\u00b2 > 0.3" = function(d) d$r2 <= 0.3,
  "CH4 R\u00b2 > 0.5" = function(d) d$r2 <= 0.5,
  "CH4 R\u00b2 > 0.7" = function(d) d$r2 <= 0.7,
  "p-value < 0.05"     = function(d) ifelse(!is.na(d$CH4_LM_pval),
                                             d$CH4_LM_pval >= 0.05, FALSE)
)

soil_seasonal <- lapply(names(soil_filters), function(fn) {
  fails <- soil_filters[[fn]](soil)
  out <- soil %>%
    mutate(fails_filter = fails) %>%
    group_by(season4) %>%
    summarise(
      n_total     = n(),
      n_flagged   = sum(fails_filter, na.rm = TRUE),
      pct_flagged = 100 * n_flagged / n_total,
      mean_temp   = mean(temp_C, na.rm = TRUE),
      .groups = "drop"
    )
  out$filter <- fn
  out
})
soil_seasonal_df <- do.call(rbind, soil_seasonal)
soil_seasonal_df$filter <- factor(soil_seasonal_df$filter,
                                   levels = names(soil_filters))

p_soil_seasonal <- ggplot(soil_seasonal_df,
                           aes(x = season4, y = pct_flagged,
                               color = filter, group = filter)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  scale_color_brewer(palette = "Dark2", name = "Filter criterion") +
  labs(
    x = NULL,
    y = "% of measurements excluded",
    title = "Quality filter exclusion rates by season \u2014 wetland soil",
    subtitle = paste0("Soil flux dataset (n=", n_soil,
                      ", ", pct_neg, "% negative). ",
                      "R\u00b2 filters exclude a larger\nfraction of cold-season measurements.")
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 9),
    plot.title = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 9, color = "grey40")
  ) +
  guides(fill = guide_legend(nrow = 1))

ggsave(file.path(fig_dir, "soil_seasonal_fraction_flagged.png"),
       p_soil_seasonal, width = 10, height = 6, dpi = 300)
ggsave(file.path(fig_dir, "soil_seasonal_fraction_flagged.pdf"),
       p_soil_seasonal, width = 10, height = 6)
message("Saved: soil_seasonal_fraction_flagged")

# Seasonal exclusion summary
cat("\n=== Soil: Exclusion rates by season ===\n")
soil_season_sum <- soil_seasonal_df %>%
  group_by(filter, season4) %>%
  summarise(pct = round(mean(pct_flagged), 1), .groups = "drop") %>%
  pivot_wider(names_from = season4, values_from = pct)
print(as.data.frame(soil_season_sum), row.names = FALSE)


# ============================================================================
# SOIL FIGURE A2: FRACTION FLAGGED BY 10°C TEMPERATURE BIN
# ============================================================================

message("\n=== Soil Figure A2: Fraction flagged by temperature bin ===")

soil$temp_bin <- cut(soil$temp_C,
                     breaks = seq(0, 30, by = 10),
                     include.lowest = TRUE,
                     right = FALSE,
                     labels = paste0(seq(0, 20, by = 10), " to ", seq(10, 30, by = 10)))

soil$temp_bin <- droplevels(soil$temp_bin)

soil_temp_flagged <- lapply(names(soil_filters), function(fn) {
  fails <- soil_filters[[fn]](soil)
  out <- soil %>%
    mutate(fails_filter = fails) %>%
    filter(!is.na(temp_bin)) %>%
    group_by(temp_bin) %>%
    summarise(
      n_total     = n(),
      n_flagged   = sum(fails_filter, na.rm = TRUE),
      pct_flagged = 100 * n_flagged / n_total,
      .groups = "drop"
    )
  out$filter <- fn
  out
})
soil_temp_flagged_df <- do.call(rbind, soil_temp_flagged)
soil_temp_flagged_df$filter <- factor(soil_temp_flagged_df$filter,
                                       levels = names(soil_filters))

p_soil_temp_frac <- ggplot(soil_temp_flagged_df,
                            aes(x = temp_bin, y = pct_flagged,
                                color = filter, group = filter)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  scale_color_brewer(palette = "Dark2", name = "Filter criterion") +
  labs(
    x = "Chamber temperature bin (C)",
    y = "% of measurements excluded",
    title = "Quality filter exclusion rates by temperature -- wetland soil",
    subtitle = paste0("Soil flux dataset (n=", n_soil,
                      ", ", pct_neg, "% negative). ",
                      "R2 filters exclude more\nmeasurements at lower temperatures.")
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 9),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    plot.title = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 9, color = "grey40")
  ) +
  guides(color = guide_legend(nrow = 1))

ggsave(file.path(fig_dir, "soil_temp_bin_fraction_flagged.png"),
       p_soil_temp_frac, width = 10, height = 6, dpi = 300)
ggsave(file.path(fig_dir, "soil_temp_bin_fraction_flagged.pdf"),
       p_soil_temp_frac, width = 10, height = 6)
message("Saved: soil_temp_bin_fraction_flagged")

cat("\n=== Soil: Exclusion rates by temperature bin ===\n")
soil_temp_sum <- soil_temp_flagged_df %>%
  group_by(filter, temp_bin) %>%
  summarise(pct = round(mean(pct_flagged), 1), .groups = "drop") %>%
  pivot_wider(names_from = temp_bin, values_from = pct)
print(as.data.frame(soil_temp_sum), row.names = FALSE)


# ============================================================================
# SOIL FIGURE B: SIGNED FLUX VS R² COLORED BY TEMPERATURE
# ============================================================================

message("\n=== Soil Figure B: Signed flux vs R\u00b2 by temperature ===")

# Compute mean temperature per R\u00b2 bin x sign for annotations
soil_r2_sign_annot <- soil %>%
  mutate(
    r2_bin = case_when(
      r2 > 0.7 ~ "R2 > 0.7",
      r2 > 0.5 ~ "0.5 < R2 <= 0.7",
      TRUE      ~ "R2 <= 0.5"
    ),
    sign = ifelse(flux >= 0, "positive", "negative")
  ) %>%
  group_by(r2_bin, sign) %>%
  summarise(
    mean_temp = round(mean(temp_C, na.rm = TRUE), 1),
    n = n(),
    .groups = "drop"
  )

cat("\nSoil: Mean temp by R\u00b2 bin x sign:\n")
print(as.data.frame(soil_r2_sign_annot), row.names = FALSE)

soil_annot_pos <- soil_r2_sign_annot %>% filter(sign == "positive")
soil_annot_neg <- soil_r2_sign_annot %>% filter(sign == "negative")

# Helper to safely look up annotation values
safe_val <- function(tbl, bin, col) {
  row <- tbl[tbl$r2_bin == bin, ]
  if (nrow(row) == 0) return(if (col == "n") 0 else NA)
  row[[col]]
}

p_soil_r2 <- ggplot(soil, aes(x = flux, y = r2)) +
  geom_point(aes(color = temp_C), alpha = 0.6, size = 1.5) +
  geom_hline(yintercept = 0.5, linetype = "dashed",
             color = "red", linewidth = 0.4) +
  geom_hline(yintercept = 0.7, linetype = "dashed",
             color = "darkred", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey50", linewidth = 0.3) +
  # Annotate: positive-flux side (right)
  annotate("label", x = 1.5, y = 0.87,
           label = paste0("R2>0.7, (+): T=", safe_val(soil_annot_pos, "R2 > 0.7", "mean_temp"),
                          "C  n=", safe_val(soil_annot_pos, "R2 > 0.7", "n")),
           size = 2.5, fill = "white", alpha = 0.8, hjust = 0, color = "darkred") +
  annotate("label", x = 1.5, y = 0.6,
           label = paste0("0.5<R2<0.7, (+): T=", safe_val(soil_annot_pos, "0.5 < R2 <= 0.7", "mean_temp"),
                          "C  n=", safe_val(soil_annot_pos, "0.5 < R2 <= 0.7", "n")),
           size = 2.5, fill = "white", alpha = 0.8, hjust = 0, color = "red") +
  annotate("label", x = 1.5, y = 0.25,
           label = paste0("R2<0.5, (+): T=", safe_val(soil_annot_pos, "R2 <= 0.5", "mean_temp"),
                          "C  n=", safe_val(soil_annot_pos, "R2 <= 0.5", "n")),
           size = 2.5, fill = "white", alpha = 0.8, hjust = 0, color = "#B2182B") +
  # Annotate: negative-flux side (left)
  annotate("label", x = -1.5, y = 0.87,
           label = paste0("R2>0.7, (-): T=", safe_val(soil_annot_neg, "R2 > 0.7", "mean_temp"),
                          "C  n=", safe_val(soil_annot_neg, "R2 > 0.7", "n")),
           size = 2.5, fill = "white", alpha = 0.8, hjust = 1, color = "darkred") +
  annotate("label", x = -1.5, y = 0.6,
           label = paste0("0.5<R2<0.7, (-): T=", safe_val(soil_annot_neg, "0.5 < R2 <= 0.7", "mean_temp"),
                          "C  n=", safe_val(soil_annot_neg, "0.5 < R2 <= 0.7", "n")),
           size = 2.5, fill = "white", alpha = 0.8, hjust = 1, color = "red") +
  annotate("label", x = -1.5, y = 0.25,
           label = paste0("R2<0.5, (-): T=", safe_val(soil_annot_neg, "R2 <= 0.5", "mean_temp"),
                          "C  n=", safe_val(soil_annot_neg, "R2 <= 0.5", "n")),
           size = 2.5, fill = "white", alpha = 0.8, hjust = 1, color = "#B2182B") +
  scale_x_continuous(
    trans = scales::pseudo_log_trans(sigma = 0.001),
    limits = c(-5, 5),
    breaks = c(-5, -1, -0.1, -0.01, 0, 0.01, 0.1, 1, 5),
    labels = c("-5", "-1", "-0.1", "-0.01", "0", "0.01", "0.1", "1", "5")
  ) +
  scale_color_viridis_c(option = "plasma", name = "Chamber temp (C)") +
  labs(
    x = expression(CH[4]~flux~(nmol~m^{-2}~s^{-1})),
    y = expression(CH[4]~R^2~(linear~model)),
    title = expression(CH[4]~R^2~vs~signed~flux~and~temperature~"--"~wetland~soil),
    subtitle = paste0("n=", n_soil, " soil measurements (", pct_neg, "% negative). ",
                      "Negative fluxes cluster below R2 thresholds at cold temperatures.")
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 9, color = "grey40")
  )

ggsave(file.path(fig_dir, "soil_r2_vs_flux_by_temperature.png"),
       p_soil_r2, width = 11, height = 6, dpi = 300)
ggsave(file.path(fig_dir, "soil_r2_vs_flux_by_temperature.pdf"),
       p_soil_r2, width = 11, height = 6)
message("Saved: soil_r2_vs_flux_by_temperature")

# --- SIGNED flux vs R² for soil ---
message("\n=== Soil Figure B2: Signed flux vs R\u00b2 by temperature ===")

p_soil_r2_signed <- ggplot(soil, aes(x = flux, y = r2)) +
  geom_point(aes(color = temp_C), alpha = 0.6, size = 1.5) +
  geom_hline(yintercept = 0.5, linetype = "dashed",
             color = "red", linewidth = 0.4) +
  geom_hline(yintercept = 0.7, linetype = "dashed",
             color = "darkred", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey50", linewidth = 0.3) +
  annotate("text", x = max(soil$flux) * 0.9,
           y = 0.52, label = "R\u00b2 = 0.5", hjust = 1, size = 3, color = "red") +
  annotate("text", x = max(soil$flux) * 0.9,
           y = 0.72, label = "R\u00b2 = 0.7", hjust = 1, size = 3, color = "darkred") +
  scale_x_continuous(
    trans = scales::pseudo_log_trans(sigma = 0.001),
    breaks = c(-5, -1, -0.1, -0.01, 0, 0.01, 0.1, 1, 5),
    labels = c("-5", "-1", "-0.1", "-0.01", "0", "0.01", "0.1", "1", "5")
  ) +
  scale_color_viridis_c(option = "plasma", name = "Chamber temp (\u00b0C)") +
  labs(
    x = expression(CH[4]~flux~(nmol~m^{-2}~s^{-1})),
    y = expression(CH[4]~R^2~(linear~model)),
    title = expression(CH[4]~R^2~vs~signed~flux~and~temperature~"\u2014"~wetland~soil),
    subtitle = paste0("n=", n_soil, " soil measurements (", pct_neg, "% negative). ",
                      "High-R\u00b2 measurements span the full flux range;\n",
                      "low-R\u00b2 values cluster near zero where sign is ambiguous.")
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 9, color = "grey40")
  )

ggsave(file.path(fig_dir, "soil_r2_vs_flux_signed_by_temperature.png"),
       p_soil_r2_signed, width = 9, height = 6, dpi = 300)
ggsave(file.path(fig_dir, "soil_r2_vs_flux_signed_by_temperature.pdf"),
       p_soil_r2_signed, width = 9, height = 6)
message("Saved: soil_r2_vs_flux_signed_by_temperature")

# --- SIGNED flux vs R², LINEAR x-axis for soil ---
message("\n=== Soil Figure B3: Signed flux vs R\u00b2 by temperature (linear scale) ===")

p_soil_r2_linear <- ggplot(soil, aes(x = flux, y = r2)) +
  geom_point(aes(color = temp_C), alpha = 0.6, size = 1.5) +
  geom_hline(yintercept = 0.5, linetype = "dashed",
             color = "red", linewidth = 0.4) +
  geom_hline(yintercept = 0.7, linetype = "dashed",
             color = "darkred", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey50", linewidth = 0.3) +
  scale_color_viridis_c(option = "plasma", name = "Chamber temp (\u00b0C)") +
  labs(
    x = expression(CH[4]~flux~(nmol~m^{-2}~s^{-1})),
    y = expression(CH[4]~R^2~(linear~model)),
    title = expression(CH[4]~R^2~vs~signed~flux~and~temperature~(linear~scale)~"\u2014"~wetland~soil),
    subtitle = paste0("n=", n_soil, " soil measurements (", pct_neg, "% negative). ",
                      "Linear scale shows the dominant negative\n",
                      "flux cluster and the asymmetry between uptake and emission magnitudes.")
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 9, color = "grey40")
  )

ggsave(file.path(fig_dir, "soil_r2_vs_flux_linear_by_temperature.png"),
       p_soil_r2_linear, width = 9, height = 6, dpi = 300)
ggsave(file.path(fig_dir, "soil_r2_vs_flux_linear_by_temperature.pdf"),
       p_soil_r2_linear, width = 9, height = 6)
message("Saved: soil_r2_vs_flux_linear_by_temperature")

cat("\n=== Soil: R\u00b2 vs |flux| and temperature ===\n")
cat("Correlation |flux| vs R\u00b2:",
    round(cor(abs(soil$flux), soil$r2, use = "complete"), 3), "\n")
cat("Correlation flux (signed) vs R\u00b2:",
    round(cor(soil$flux, soil$r2, use = "complete"), 3), "\n")
cat("Correlation temp vs R\u00b2:",
    round(cor(soil$temp_C, soil$r2, use = "complete"), 3), "\n")


# ============================================================================
# SOIL FIGURE C: TIME-SERIES COLORED BY FILTER STATUS (R²-only)
# ============================================================================

message("\n=== Soil Figure C: Time-series by R\u00b2 filter status ===")

soil$fails_r2_05 <- soil$r2 <= 0.5
soil$fails_r2_07 <- soil$r2 <= 0.7

# Use ASCII-safe labels to avoid Unicode encoding mismatches in ggplot
soil_lbl_pass  <- "Passes R2 > 0.7"
soil_lbl_mid   <- "Fails R2 > 0.7 only"
soil_lbl_fail  <- "Fails R2 > 0.5"

soil$filter_status <- case_when(
  soil$fails_r2_05  ~ soil_lbl_fail,
  soil$fails_r2_07  ~ soil_lbl_mid,
  TRUE               ~ soil_lbl_pass
)

soil$filter_status <- factor(soil$filter_status,
                              levels = c(soil_lbl_pass,
                                         soil_lbl_mid,
                                         soil_lbl_fail))

soil_status_colors <- setNames(
  c("#2166AC", "#F4A582", "#B2182B"),
  c(soil_lbl_pass, soil_lbl_mid, soil_lbl_fail)
)

p_soil_ts <- ggplot(soil, aes(x = date_parsed, y = flux)) +
  geom_point(aes(color = filter_status), alpha = 0.6, size = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  scale_color_manual(values = soil_status_colors, name = "Filter status") +
  scale_x_date(date_breaks = "2 months", date_labels = "%b\n%Y") +
  scale_y_continuous(
    trans = scales::pseudo_log_trans(sigma = 0.001),
    breaks = c(-1, -0.1, -0.01, 0, 0.01, 0.1, 1, 5),
    labels = c("-1", "-0.1", "-0.01", "0", "0.01", "0.1", "1", "5")
  ) +
  labs(
    x = NULL,
    y = expression(CH[4]~flux~(nmol~m^{-2}~s^{-1})),
    title = "Temporal distribution of measurements failing R\u00b2 filters \u2014 wetland soil",
    subtitle = paste0("n=", n_soil, " soil measurements. ",
                      "Excluded measurements (red tones) cluster in cold months,\n",
                      "not randomly distributed in time.")
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 9, color = "grey40")
  )

ggsave(file.path(fig_dir, "soil_timeseries_filter_status.png"),
       p_soil_ts, width = 12, height = 6, dpi = 300)
ggsave(file.path(fig_dir, "soil_timeseries_filter_status.pdf"),
       p_soil_ts, width = 12, height = 6)
message("Saved: soil_timeseries_filter_status")

cat("\n=== Soil: Filter status by season ===\n")
soil_status_season <- soil %>%
  group_by(season, filter_status) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(season) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  ungroup()
print(as.data.frame(soil_status_season), row.names = FALSE)


# ============================================================================
# SOIL FIGURE D: FLUX VS TEMPERATURE COLORED BY FILTER STATUS
# ============================================================================

message("\n=== Soil Figure D: Flux vs temperature by filter status ===")

p_soil_flux_temp <- ggplot(soil, aes(x = temp_C, y = flux)) +
  geom_point(aes(color = filter_status), alpha = 0.6, size = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey70", linewidth = 0.3) +
  scale_color_manual(values = soil_status_colors, name = "Filter status") +
  scale_y_continuous(
    trans = scales::pseudo_log_trans(sigma = 0.001),
    breaks = c(-1, -0.1, -0.01, 0, 0.01, 0.1, 1, 5),
    labels = c("-1", "-0.1", "-0.01", "0", "0.01", "0.1", "1", "5")
  ) +
  labs(
    x = "Chamber temperature (\u00b0C)",
    y = expression(CH[4]~flux~(nmol~m^{-2}~s^{-1})),
    title = "Quality filter exclusions are structured by temperature \u2014 wetland soil",
    subtitle = paste0("n=", n_soil, " measurements. ",
                      "Excluded measurements (red tones) cluster at low temperatures.\n",
                      "Filters: R\u00b2 > 0.5 and R\u00b2 > 0.7 thresholds.")
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 9, color = "grey40")
  )

ggsave(file.path(fig_dir, "soil_flux_vs_temperature_filter_status.png"),
       p_soil_flux_temp, width = 9, height = 6, dpi = 300)
ggsave(file.path(fig_dir, "soil_flux_vs_temperature_filter_status.pdf"),
       p_soil_flux_temp, width = 9, height = 6)
message("Saved: soil_flux_vs_temperature_filter_status")

cat("\n=== Soil: Temperature of excluded vs retained ===\n")
soil_excl <- soil$fails_r2_05
cat("  Retained (R2>0.5): n=", sum(!soil_excl),
    ", mean temp=", round(mean(soil$temp_C[!soil_excl], na.rm=TRUE), 1), "\u00b0C",
    ", median flux=", round(median(soil$flux[!soil_excl]), 4), "\n")
cat("  Excluded (R2<=0.5): n=", sum(soil_excl),
    ", mean temp=", round(mean(soil$temp_C[soil_excl], na.rm=TRUE), 1), "\u00b0C",
    ", median flux=", round(median(soil$flux[soil_excl]), 4), "\n")


# ============================================================================
# SOIL FIGURE E: CUMULATIVE ANNUAL BUDGET UNDER R² FILTERS
# ============================================================================

message("\n=== Soil Figure E: Cumulative annual budget ===")

# The soil data spans ~June 2020 to ~May 2021 (one "field year")
# Use all available months for the cumulative budget
soil_sorted <- soil %>% arrange(date_parsed)

soil_treatments <- list(
  "No filter (keep all)"  = function(d) rep(TRUE, nrow(d)),
  "R\u00b2 > 0.3 (remove)"   = function(d) d$r2 > 0.3,
  "R\u00b2 > 0.5 (remove)"   = function(d) d$r2 > 0.5,
  "R\u00b2 > 0.7 (remove)"   = function(d) d$r2 > 0.7
)

# Use the year_month ordering so cumulation follows calendar time
unique_ym <- sort(unique(soil_sorted$year_month))

soil_budget <- lapply(names(soil_treatments), function(trt) {
  keep <- soil_treatments[[trt]](soil_sorted)
  monthly <- soil_sorted %>%
    mutate(keep = keep) %>%
    filter(keep) %>%
    group_by(year_month) %>%
    summarise(
      mean_flux = mean(flux, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    )

  # Ensure all months are present
  all_ym <- data.frame(year_month = unique_ym, stringsAsFactors = FALSE)
  monthly <- merge(all_ym, monthly, by = "year_month", all.x = TRUE)
  monthly$mean_flux[is.na(monthly$mean_flux)] <- 0
  monthly$n[is.na(monthly$n)] <- 0
  monthly <- monthly[order(monthly$year_month), ]

  monthly$treatment <- trt
  monthly$cum_flux  <- cumsum(monthly$mean_flux)
  monthly$month_idx <- seq_len(nrow(monthly))

  # Create a date midpoint for the x-axis
  monthly$date_mid <- as.Date(paste0(monthly$year_month, "-15"))
  monthly
})
soil_budget_df <- do.call(rbind, soil_budget)
soil_budget_df$treatment <- factor(soil_budget_df$treatment,
                                    levels = names(soil_treatments))

soil_trt_colors <- c(
  "No filter (keep all)"     = "black",
  "R\u00b2 > 0.3 (remove)"  = "#E69F00",
  "R\u00b2 > 0.5 (remove)"  = "#D6604D",
  "R\u00b2 > 0.7 (remove)"  = "#B2182B"
)

p_soil_cum <- ggplot(soil_budget_df,
                      aes(x = date_mid, y = cum_flux,
                          color = treatment, linetype = treatment)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  scale_x_date(date_breaks = "2 months", date_labels = "%b\n%Y") +
  scale_color_manual(values = soil_trt_colors, name = "Treatment") +
  scale_linetype_manual(
    values = c("No filter (keep all)"    = "solid",
               "R\u00b2 > 0.3 (remove)" = "dashed",
               "R\u00b2 > 0.5 (remove)" = "dashed",
               "R\u00b2 > 0.7 (remove)" = "dashed"),
    name = "Treatment"
  ) +
  labs(
    x = NULL,
    y = expression(Cumulative~mean~CH[4]~flux~(nmol~m^{-2}~s^{-1})),
    title = "Cumulative CH\u2084 flux estimate diverges with R\u00b2 filter \u2014 wetland soil",
    subtitle = paste0("Soil data, ", min(soil$date), " to ", max(soil$date),
                      " (n=", n_soil, "). ",
                      pct_neg, "% of fluxes are negative.\n",
                      "Removing low-R\u00b2 measurements inflates the annual estimate ",
                      "by preferentially excluding near-zero dormant fluxes.")
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 9),
    plot.title = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 9, color = "grey40")
  ) +
  guides(color = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2))

ggsave(file.path(fig_dir, "soil_cumulative_annual_budget.png"),
       p_soil_cum, width = 10, height = 6, dpi = 300)
ggsave(file.path(fig_dir, "soil_cumulative_annual_budget.pdf"),
       p_soil_cum, width = 10, height = 6)
message("Saved: soil_cumulative_annual_budget")

cat("\n=== Soil: Final cumulative flux by treatment ===\n")
soil_final <- soil_budget_df %>%
  filter(month_idx == max(month_idx)) %>%
  select(treatment, cum_flux, n)
soil_final$pct_change <- round(100 * (soil_final$cum_flux / soil_final$cum_flux[1] - 1), 1)
print(as.data.frame(soil_final), row.names = FALSE)


# ============================================================================
# SOIL: LANDSCAPE POSITION BREAKDOWN
# ============================================================================

message("\n=== Soil: Exclusion by landscape position ===")

if ("landscape_position" %in% names(soil) && length(unique(soil$landscape_position)) > 1) {
  cat("\n=== R\u00b2 > 0.5 exclusion by landscape position ===\n")
  lp_summary <- soil %>%
    mutate(excluded_r2 = r2 <= 0.5) %>%
    group_by(landscape_position) %>%
    summarise(
      n = n(),
      pct_excluded = round(100 * mean(excluded_r2), 1),
      pct_negative = round(100 * mean(flux < 0), 1),
      median_flux  = round(median(flux), 4),
      mean_temp    = round(mean(temp_C, na.rm = TRUE), 1),
      .groups = "drop"
    )
  print(as.data.frame(lp_summary), row.names = FALSE)
}


# ============================================================================
# COMBINED SUMMARY
# ============================================================================

cat("\n\n=== COMBINED SUMMARY: Seasonal confounding ===\n")
cat("\nStem data (Harvard Forest, n=", stem_n, "):\n")
cat("  Dormant exclusion (R\u00b2>0.5 or Chr95%): ", stem_pct_dorm_excl, "%\n")
cat("  Growing exclusion: ", stem_pct_grow_excl, "%\n")

cat("\nSoil data (YMF wetland, n=", n_soil, "):\n")
soil_excl_grow <- soil$fails_r2_05[soil$season == "Growing (May\u2013Oct)"]
soil_excl_dorm <- soil$fails_r2_05[soil$season == "Dormant (Nov\u2013Apr)"]
cat("  Dormant exclusion (R\u00b2>0.5): ",
    round(100 * mean(soil_excl_dorm), 1), "%\n")
cat("  Growing exclusion: ",
    round(100 * mean(soil_excl_grow), 1), "%\n")

message("\n=== All seasonal confounding analyses complete ===")
