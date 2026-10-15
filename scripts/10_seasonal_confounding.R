# ============================================================================
# 10_seasonal_confounding.R — Seasonal structure in quality filter exclusions
#
# Demonstrates that quality filters (R², MDF, SNR) do not exclude measurements
# randomly but preferentially exclude cold-season, low-flux conditions ---
# introducing structured, ecologically correlated bias.
#
# Stem and wetland-soil analyses (summary tables printed to the console):
#   A. Fraction flagged by season and by temperature bin per filter criterion
#   B. R² vs flux sign and temperature
#   C. Filter status (R² and MDF) by season
#   D. Temperature of excluded vs retained measurements
#   E. Cumulative annual budget under different filter approaches
# The plot data are saved for 19_main_figures.R (Figs 6 and 7), 20_generality.R
# (soil) and 22_transform_inference.R (df_gla).
#
# Loads: outputs/03_mdf_results.RData, outputs/01_loaded_data.RData
# Saves: outputs/10_seasonal_confounding.RData
# ============================================================================

source("scripts/00_setup.R")

message("\n=== Loading MDF results ===")
load(file.path(output_dir, "03_mdf_results.RData"))

# --- Prepare date columns ---
df$date_parsed <- as.Date(df$date)
df$month_num   <- as.numeric(format(df$date_parsed, "%m"))
df$month_lab   <- factor(month.abb[df$month_num], levels = month.abb)
df$year_month  <- format(df$date_parsed, "%Y-%m")

# Instrument label
df$inst_label <- ifelse(df$year == 2025, "Precise analyzer",
                         "Noisier analyzer")

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
# A: FRACTION FLAGGED BY SEASON PER FILTER CRITERION
# ============================================================================

message("\n=== A: Seasonal fraction flagged ===")

# Select key filter criteria to show
key_filters <- list(
  "CH4 R\u00b2 > 0.5"     = function(d) d$CH4_r2 <= 0.5,
  "CH4 R\u00b2 > 0.7"     = function(d) d$CH4_r2 <= 0.7,
  "Manufacturer MDF"       = function(d) ifelse(!is.na(d$CH4_below_MDF_manuf),
                                                  d$CH4_below_MDF_manuf, FALSE),
  "Campaign \u03c3 95%"           = function(d) ifelse(!is.na(d$CH4_below_MDF_wass95),
                                                  d$CH4_below_MDF_wass95, FALSE),
  "Per-closure \u03c3 \u00d73t 95%"       = function(d) ifelse(!is.na(d$CH4_below_MDF_chr95),
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
# A2: FRACTION FLAGGED BY 10°C TEMPERATURE BIN
# ============================================================================

message("\n=== A2: Fraction flagged by temperature bin ===")

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

cat("\n=== Exclusion rates by temperature bin ===\n")
temp_summary <- temp_flagged_df %>%
  group_by(filter, temp_bin) %>%
  summarise(pct = round(mean(pct_flagged), 1), .groups = "drop") %>%
  pivot_wider(names_from = temp_bin, values_from = pct)
print(as.data.frame(temp_summary), row.names = FALSE)


# ============================================================================
# B: R² VS FLUX SIGN AND TEMPERATURE
# ============================================================================

message("\n=== B: R\u00b2 vs flux and temperature ===")

# Use Noisier analyzer only (Precise analyzer covers only May-Oct, biasing the seasonal pattern)
df_gla <- df %>% filter(inst_label != "Precise analyzer")
message("  Noisier analyzer measurements: ", nrow(df_gla))

# Mean temperature per R\u00b2 bin x sign
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
# C: FILTER STATUS (R² AND MDF) BY SEASON
# ============================================================================

message("\n=== C: Filter status by season ===")

# Classify each measurement by filter status
# Use two key filters: R² > 0.5 and the reference (campaign-sigma 95%) MDF
df$fails_r2 <- df$CH4_r2 <= 0.5
df$fails_mdf <- ifelse(!is.na(df$CH4_below_MDF_wass95),
                        df$CH4_below_MDF_wass95, FALSE)

df$filter_status <- case_when(
  df$fails_r2 & df$fails_mdf  ~ "Fails both",
  df$fails_r2 & !df$fails_mdf ~ "Fails R\u00b2 only",
  !df$fails_r2 & df$fails_mdf ~ "Fails MDF only",
  TRUE                         ~ "Passes all"
)

df$filter_status <- factor(df$filter_status,
                            levels = c("Passes all", "Fails MDF only",
                                       "Fails R\u00b2 only", "Fails both"))

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
# D: TEMPERATURE OF EXCLUDED VS RETAINED MEASUREMENTS
# ============================================================================

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
# E: CUMULATIVE ANNUAL BUDGET UNDER DIFFERENT FILTER APPROACHES
# ============================================================================

message("\n=== E: Cumulative annual budget ===")

# Use calendar year 2024 (Noisier analyzer only, all 12 months represented)
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
  "Campaign \u03c3 95%\n(remove)"  = function(d) {
    keep <- !ifelse(!is.na(d$CH4_below_MDF_wass95), d$CH4_below_MDF_wass95, FALSE)
    d$CH4_flux_nmolpm2ps[keep]
  },
  "Campaign \u03c3 95%\n(keep all, flag)" = function(d) d$CH4_flux_nmolpm2ps
)

# Compute monthly mean flux under each treatment
monthly_budget <- lapply(names(treatments), function(trt) {
  vals <- treatments[[trt]](df_2024)

  if (trt %in% c("No filter (keep all)", "Campaign \u03c3 95%\n(keep all, flag)")) {
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
      keep <- !ifelse(!is.na(df_2024$CH4_below_MDF_wass95),
                       df_2024$CH4_below_MDF_wass95, FALSE)
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
# SOIL A: FRACTION FLAGGED BY SEASON (R²-based filters only)
# ============================================================================

message("\n=== Soil A: Seasonal fraction flagged ===")

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

# Seasonal exclusion summary
cat("\n=== Soil: Exclusion rates by season ===\n")
soil_season_sum <- soil_seasonal_df %>%
  group_by(filter, season4) %>%
  summarise(pct = round(mean(pct_flagged), 1), .groups = "drop") %>%
  pivot_wider(names_from = season4, values_from = pct)
print(as.data.frame(soil_season_sum), row.names = FALSE)


# ============================================================================
# SOIL A2: FRACTION FLAGGED BY 10°C TEMPERATURE BIN
# ============================================================================

message("\n=== Soil A2: Fraction flagged by temperature bin ===")

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

cat("\n=== Soil: Exclusion rates by temperature bin ===\n")
soil_temp_sum <- soil_temp_flagged_df %>%
  group_by(filter, temp_bin) %>%
  summarise(pct = round(mean(pct_flagged), 1), .groups = "drop") %>%
  pivot_wider(names_from = temp_bin, values_from = pct)
print(as.data.frame(soil_temp_sum), row.names = FALSE)


# ============================================================================
# SOIL B: R² VS FLUX SIGN AND TEMPERATURE
# ============================================================================

message("\n=== Soil B: R\u00b2 vs flux and temperature ===")

# Mean temperature per R\u00b2 bin x sign
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

cat("\n=== Soil: R\u00b2 vs |flux| and temperature ===\n")
cat("Correlation |flux| vs R\u00b2:",
    round(cor(abs(soil$flux), soil$r2, use = "complete"), 3), "\n")
cat("Correlation flux (signed) vs R\u00b2:",
    round(cor(soil$flux, soil$r2, use = "complete"), 3), "\n")
cat("Correlation temp vs R\u00b2:",
    round(cor(soil$temp_C, soil$r2, use = "complete"), 3), "\n")


# ============================================================================
# SOIL C: R² FILTER STATUS BY SEASON
# ============================================================================

message("\n=== Soil C: R\u00b2 filter status by season ===")

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

cat("\n=== Soil: Filter status by season ===\n")
soil_status_season <- soil %>%
  group_by(season, filter_status) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(season) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  ungroup()
print(as.data.frame(soil_status_season), row.names = FALSE)


# ============================================================================
# SOIL D: TEMPERATURE OF EXCLUDED VS RETAINED MEASUREMENTS
# ============================================================================

cat("\n=== Soil: Temperature of excluded vs retained ===\n")
soil_excl <- soil$fails_r2_05
cat("  Retained (R2>0.5): n=", sum(!soil_excl),
    ", mean temp=", round(mean(soil$temp_C[!soil_excl], na.rm=TRUE), 1), "\u00b0C",
    ", median flux=", round(median(soil$flux[!soil_excl]), 4), "\n")
cat("  Excluded (R2<=0.5): n=", sum(soil_excl),
    ", mean temp=", round(mean(soil$temp_C[soil_excl], na.rm=TRUE), 1), "\u00b0C",
    ", median flux=", round(median(soil$flux[soil_excl]), 4), "\n")


# ============================================================================
# SOIL E: CUMULATIVE ANNUAL BUDGET UNDER R² FILTERS
# ============================================================================

message("\n=== Soil E: Cumulative annual budget ===")

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

  # Date midpoint of each month
  monthly$date_mid <- as.Date(paste0(monthly$year_month, "-15"))
  monthly
})
soil_budget_df <- do.call(rbind, soil_budget)
soil_budget_df$treatment <- factor(soil_budget_df$treatment,
                                    levels = names(soil_treatments))

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

# Plot-ready data for the assembled main figures (19_main_figures.R)
save(seasonal_flagged_df, temp_flagged_df, soil_seasonal_df, soil_temp_flagged_df,
     monthly_budget_df, soil_budget_df, df_gla, soil,
     file = file.path(output_dir, "10_seasonal_confounding.RData"))
