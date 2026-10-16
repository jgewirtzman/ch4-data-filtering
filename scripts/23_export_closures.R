# ============================================================================
# 23_export_closures.R — Closure-level results table for the stems
#
# One row per stem closure analysed in the paper (1547), with the flux, the
# precision, the MDFs, the detection class, the regression-based detection and
# the closure screens. Written to results/stem_closures_processed.csv (tracked,
# archived with the code); column definitions in results/README.md.
#
# Loads:  outputs/03_mdf_results.RData
# Writes: results/stem_closures_processed.csv
# ============================================================================

source("scripts/00_setup.R")
load(file.path(output_dir, "03_mdf_results.RData"))

dir.create("results", showWarnings = FALSE)

# Regression interval of the slope (Cowan et al. 2025) with the campaign sigma, for
# n evenly spaced readings dt seconds apart: sum((t - tbar)^2) = dt^2 n (n^2 - 1) / 12
sxx <- with(df, dt_s^2 * n_pts * (n_pts^2 - 1) / 12)
reg_half <- with(df, ifelse(n_pts > 2, qt(0.975, n_pts - 2), NA_real_) * sigma_camp_CH4 / sqrt(sxx) * flux_term)

out <- df %>%
  transmute(
    closure_id         = stem_id,
    date               = as.character(date),
    site               = location,
    plot               = PLOT,
    tree               = Tree,
    species            = SPECIES,
    analyzer           = if_else(year < 2025, "noisier (LGR GLA131)", "precise (LI-COR LI-7810)"),
    air_temp_C         = airt,
    n_readings         = n_pts,
    logging_interval_s = dt_s,
    closure_s          = t_sec,
    flux_term          = flux_term,
    CH4_flux           = CH4_flux_nmolpm2ps,
    CH4_model          = CH4_model,
    CH4_r2             = CH4_r2,
    CO2_flux           = CO2_flux_umolpm2ps,
    sigma_closure_ppb  = sc_CH4,
    sigma_campaign_ppb = sigma_camp_CH4,
    MDF_reference      = CH4_MDF_ref,
    MDF_datasheet      = CH4_MDF_manufacturer,
    detection_class    = case_when(is.na(CH4_MDF_ref) ~ NA_character_,
                                   CH4_flux_nmolpm2ps >  CH4_MDF_ref ~ "emission",
                                   CH4_flux_nmolpm2ps < -CH4_MDF_ref ~ "uptake",
                                   TRUE ~ "below MDF"),
    reg_half_95        = reg_half,
    lm_p               = lm_p_CH4,
    resid_ac1          = lm_ac1_CH4,
    across(c(qc_c0, qc_co2_tracer, qc_convex, qc_min_obs, qc_ambient, qc_noisy, qc_any), ~ .x)
  ) %>%
  arrange(date, closure_id)

write.csv(out, file.path("results", "stem_closures_processed.csv"), row.names = FALSE, na = "")
message("Wrote results/stem_closures_processed.csv: ", nrow(out), " closures, ", ncol(out), " columns")
print(table(out$analyzer, out$detection_class, useNA = "ifany"))
