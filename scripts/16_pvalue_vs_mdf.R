# ============================================================================
# 16_pvalue_vs_mdf.R — Regression p-value vs MDF as the detection criterion
#
# The slope t-test of a 1 Hz closure (|flux| / SE, or goFlux's LM.p.val) treats
# consecutive residuals as independent. They are not (instrument drift, chamber
# mixing), so nearly every closure "detects" a flux by p < 0.05. This script
# quantifies the over-detection relative to the reference MDF (campaign sigma,
# 95 %) on the stem and canopy datasets, and inspects the negative fluxes:
# which are significant by p, which clear the MDF, and whether the concurrent
# CO2 flux is positive (a coherent closure).
#
# Loads: outputs/03_mdf_results.RData
# Saves: outputs/16_pvalue_vs_mdf.csv, outputs/16_negatives.csv,
#        outputs/16_pvalue_vs_mdf.RData (plotted as Fig. 4 by 19_main_figures.R)
# ============================================================================

source("scripts/00_setup.R")

message("\n=== Loading MDF results from 03 ===")
load(file.path(output_dir, "03_mdf_results.RData"))

# ── Stem dataset: OLS slope p-value refitted on the raw trace (script 02) ───
# (the SE stored in the source flux table is not reliable for the noisier analyzer;
# lm_p_CH4 and the residual lag-1 autocorrelation come straight from the concentration trace)
stem <- df %>%
  filter(!is.na(CH4_MDF_ref), !is.na(lm_p_CH4)) %>%
  mutate(
    dataset   = "HF stem",
    analyzer  = if_else(year == 2025, "Precise analyzer", "Noisier analyzer"),
    t_stat    = NA_real_,
    p_val     = lm_p_CH4,
    ac1       = lm_ac1_CH4,
    detect_p  = p_val < 0.05,
    detect_mdf = !CH4_below_MDF_ref,
    flux      = CH4_flux_nmolpm2ps,
    co2_flux  = CO2_flux_umolpm2ps
  ) %>%
  select(dataset, analyzer, UniqueID, flux, co2_flux, CH4_MDF_ref, t_stat, p_val, ac1,
         detect_p, detect_mdf, CH4_det_class)

# ── Canopy datasets: goFlux linear-model p-value ───────────────────────────
canopy_tbl <- function(out, dataset) {
  out %>%
    filter(!is.na(CH4_MDF_ref), !is.na(CH4_LM.p.val)) %>%
    transmute(
      dataset, analyzer = paste0("A-", sub("LGR", "", instrument)),
      UniqueID, flux = CH4_best.flux, co2_flux = CO2_best.flux, CH4_MDF_ref,
      t_stat = NA_real_, p_val = CH4_LM.p.val, ac1 = NA_real_,
      detect_p = CH4_LM.p.val < 0.05, detect_mdf = !CH4_below_MDF_ref,
      CH4_det_class)
}
canopy <- bind_rows(canopy_tbl(hf_out, "HF canopy"), canopy_tbl(ymf_out, "YMF oak"))
canopy$analyzer[canopy$dataset == "YMF oak"] <- "A-3 (YMF)"

all_det <- bind_rows(stem, canopy)

# ── Agreement table ────────────────────────────────────────────────────────
agree <- all_det %>%
  group_by(dataset, analyzer) %>%
  summarise(
    n            = n(),
    pct_detect_p   = round(100 * mean(detect_p), 1),
    pct_detect_mdf = round(100 * mean(detect_mdf), 1),
    n_p_only     = sum(detect_p & !detect_mdf),    # significant by p, below MDF
    n_mdf_only   = sum(!detect_p & detect_mdf),    # above MDF, not significant
    n_both       = sum(detect_p & detect_mdf),
    n_neither    = sum(!detect_p & !detect_mdf),
    .groups = "drop"
  )
overall <- all_det %>%
  group_by(dataset) %>%
  summarise(analyzer = "all", n = n(),
            pct_detect_p = round(100 * mean(detect_p), 1),
            pct_detect_mdf = round(100 * mean(detect_mdf), 1),
            n_p_only = sum(detect_p & !detect_mdf), n_mdf_only = sum(!detect_p & detect_mdf),
            n_both = sum(detect_p & detect_mdf), n_neither = sum(!detect_p & !detect_mdf),
            .groups = "drop")
agree <- bind_rows(agree, overall)
cat("\n=== Residual lag-1 autocorrelation of the slope fit (stem traces) ===\n")
print(stem %>% group_by(analyzer) %>%
        summarise(n = n(), ac1_median = round(median(ac1, na.rm = TRUE), 2),
                  ac1_q25 = round(quantile(ac1, 0.25, na.rm = TRUE), 2),
                  ac1_q75 = round(quantile(ac1, 0.75, na.rm = TRUE), 2), .groups = "drop") %>%
        as.data.frame(), row.names = FALSE)
cat("\n=== Detection by regression p < 0.05 vs reference MDF (campaign sigma 95%) ===\n")
print(as.data.frame(agree), row.names = FALSE)
write.csv(agree, file.path(output_dir, "16_pvalue_vs_mdf.csv"), row.names = FALSE)

# ── Negative fluxes: noise, or coherent uptake? ────────────────────────────
neg <- all_det %>%
  filter(flux < 0) %>%
  mutate(co2_rising = !is.na(co2_flux) & co2_flux > 0) %>%
  group_by(dataset, analyzer) %>%
  summarise(
    n_neg          = n(),
    n_p_sig        = sum(detect_p),
    n_above_mdf    = sum(detect_mdf),
    n_above_mdf_co2_rising = sum(detect_mdf & co2_rising),
    pct_above_mdf  = round(100 * mean(detect_mdf), 1),
    .groups = "drop"
  )
cat("\n=== Negative CH4 fluxes: significance vs MDF vs CO2 coherence ===\n")
print(as.data.frame(neg), row.names = FALSE)
write.csv(neg, file.path(output_dir, "16_negatives.csv"), row.names = FALSE)

save(all_det, agree, neg, file = file.path(output_dir, "16_pvalue_vs_mdf.RData"))
message("\n=== p-value vs MDF comparison complete ===")
