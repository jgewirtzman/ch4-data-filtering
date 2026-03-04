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
df$inst_label <- ifelse(df$year == 2025, "LI-7810 (2025)",
                         "GLA131 (2023\u201324)")

# Season
df$season <- ifelse(df$month_num %in% 5:10, "Growing (May\u2013Oct)",
                     "Dormant (Nov\u2013Apr)")


# ============================================================================
# FIGURE A: SEASONAL FRACTION FLAGGED BY MONTH PER FILTER CRITERION
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

# Compute per-month fraction flagged for each criterion
monthly_flagged <- lapply(names(key_filters), function(fn) {
  fails <- key_filters[[fn]](df)
  out <- df %>%
    mutate(fails_filter = fails) %>%
    group_by(month_num, month_lab) %>%
    summarise(
      n_total  = n(),
      n_flagged = sum(fails_filter, na.rm = TRUE),
      pct_flagged = 100 * n_flagged / n_total,
      mean_airt = mean(airt, na.rm = TRUE),
      .groups = "drop"
    )
  out$filter <- fn
  out
})
monthly_flagged_df <- do.call(rbind, monthly_flagged)

# Order filters for legend
monthly_flagged_df$filter <- factor(monthly_flagged_df$filter,
                                     levels = names(key_filters))

# Add background shading for dormant months
dormant_rects <- data.frame(
  xmin = c(0.5, 10.5),
  xmax = c(4.5, 12.5),
  ymin = -Inf, ymax = Inf
)

p_seasonal_frac <- ggplot(monthly_flagged_df,
                           aes(x = month_num, y = pct_flagged,
                               color = filter, group = filter)) +
  geom_rect(data = dormant_rects,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE, fill = "grey90", alpha = 0.5) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +
  scale_color_brewer(palette = "Dark2", name = "Filter criterion") +
  labs(
    x = NULL,
    y = "% of measurements excluded",
    title = "Quality filter exclusion rates by month",
    subtitle = paste0("Stem flux dataset (n=", nrow(df),
                      "). Grey bands = dormant season (Nov\u2013Apr).\n",
                      "Exclusion rates peak in cold months when true fluxes approach zero.")
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

# Print summary: growing vs dormant exclusion rates
cat("\n=== Growing vs dormant exclusion rates ===\n")
season_summary <- monthly_flagged_df %>%
  mutate(season = ifelse(month_num %in% 5:10, "Growing", "Dormant")) %>%
  group_by(filter, season) %>%
  summarise(
    mean_pct_flagged = round(mean(pct_flagged), 1),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = season, values_from = mean_pct_flagged,
              names_prefix = "pct_")
print(as.data.frame(season_summary), row.names = FALSE)


# ============================================================================
# FIGURE B: FLUX MAGNITUDE VS R² COLORED BY TEMPERATURE
# ============================================================================

message("\n=== Figure B: |Flux| vs R\u00b2 by temperature ===")

# Temperature bins for coloring
df$temp_bin <- cut(df$airt,
                    breaks = c(-Inf, 0, 10, 20, Inf),
                    labels = c("< 0\u00b0C", "0\u201310\u00b0C",
                               "10\u201320\u00b0C", "> 20\u00b0C"))

p_r2_flux <- ggplot(df, aes(x = abs(CH4_flux_nmolpm2ps), y = CH4_r2)) +
  geom_point(aes(color = airt), alpha = 0.5, size = 1.2) +
  geom_hline(yintercept = 0.5, linetype = "dashed",
             color = "red", linewidth = 0.4) +
  geom_hline(yintercept = 0.7, linetype = "dashed",
             color = "darkred", linewidth = 0.4) +
  annotate("text", x = max(abs(df$CH4_flux_nmolpm2ps)) * 0.95,
           y = 0.52, label = "R\u00b2 = 0.5", hjust = 1, size = 3, color = "red") +
  annotate("text", x = max(abs(df$CH4_flux_nmolpm2ps)) * 0.95,
           y = 0.72, label = "R\u00b2 = 0.7", hjust = 1, size = 3, color = "darkred") +
  scale_x_continuous(
    trans = scales::pseudo_log_trans(sigma = 0.01),
    breaks = c(0, 0.01, 0.1, 1, 10, 50),
    labels = c("0", "0.01", "0.1", "1", "10", "50")
  ) +
  scale_color_viridis_c(option = "plasma", name = "Air temp (\u00b0C)") +
  facet_wrap(~ inst_label, ncol = 2) +
  labs(
    x = expression("|"*CH[4]~flux*"|"~(nmol~m^{-2}~s^{-1})),
    y = expression(CH[4]~R^2),
    title = expression(CH[4]~R^2~is~correlated~with~flux~magnitude~and~temperature),
    subtitle = paste0("Low R\u00b2 measurements are systematically colder and lower-flux.\n",
                      "R\u00b2 filters preferentially exclude dormant/cold-season measurements.")
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "right",
    strip.text = element_text(size = 11, face = "bold"),
    plot.title = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 9, color = "grey40")
  )

ggsave(file.path(fig_dir, "r2_vs_flux_by_temperature.png"),
       p_r2_flux, width = 12, height = 6, dpi = 300)
ggsave(file.path(fig_dir, "r2_vs_flux_by_temperature.pdf"),
       p_r2_flux, width = 12, height = 6)
message("Saved: r2_vs_flux_by_temperature")

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

# Use calendar year 2024 (GLA131 only, all 12 months represented)
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
    subtitle = paste0("GLA131 stem data, calendar year 2024 (n=", n_2024,
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

message("\n=== Seasonal confounding analysis complete ===")
