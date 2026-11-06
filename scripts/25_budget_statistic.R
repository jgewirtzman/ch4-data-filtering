# ============================================================================
# 25_budget_statistic.R — Which statistic scales to a budget, and what removing
#                         large (hot-spot) fluxes does to it
#
# A budget is a sum (flux x area x time), so it scales from the arithmetic mean of all
# closures. Here the same closures are summarised by the median, by the geometric mean of
# the positive fluxes (back-transformed mean of log10), and by the mean after an outlier
# rule (Tukey fences, 1.5 x IQR beyond the quartiles) has removed the largest values, and
# compared with the mean of all closures and with removal of closures below the reference
# MDF, with R2 <= 0.5 or with a slope p-value >= 0.05. Removal for closure-level evidence of failure
# (flat CO2 tracer; stems only, Sect. S5) is shown for contrast.
#
# Loads:  outputs/03_mdf_results.RData, 10_seasonal_confounding.RData, 20_generality.RData,
#         data/from_agsoil_study/, data/from_coastal_wetland_study/coastal_closures.csv (if present)
# Writes: outputs/25_budget_statistic.RData (stem_cum, soil_cum, conc, stat_cmp, hot_trees, stem_trees),
#         outputs/25_budget_statistic_table.csv
# ============================================================================

source("scripts/00_setup.R")
ld <- function(f) { e <- new.env(); load(file.path(output_dir, f), envir = e); e }

gmean_pos <- function(v) { v <- v[v > 0]; if (length(v)) 10^mean(log10(v)) else NA_real_ }
in_fences <- function(v) { q <- quantile(v, c(0.25, 0.75)); i <- diff(q); v >= q[1] - 1.5 * i & v <= q[2] + 1.5 * i }
stats_of <- function(v, keep = in_fences(v)) c(mean = mean(v), median = median(v), gmean = gmean_pos(v), trimmed = mean(v[keep]))
stat_lab <- c(mean = "Mean, all closures", median = "Median", gmean = "Geometric mean, > 0",
              trimmed = "Mean, outliers removed")
to_mmol <- function(x, date) x * lubridate::days_in_month(date) * 86400 / 1e6   # nmol m-2 s-1 over a month -> mmol m-2

# ── Stems, 2024 (noisier analyzer; the closures of Fig. 7b), month by month ──
a03 <- ld("03_mdf_results.RData")$df
a10 <- ld("10_seasonal_confounding.RData")
# outliers are flagged against the fences of the whole dataset (all closures of that analyzer, or of
# that part of the soil gradient), as an outlier rule would be applied in practice, then each month
# is summarised
s24 <- a10$df_gla %>% mutate(in_f = in_fences(CH4_flux_nmolpm2ps)) %>% filter(year == 2024)
stem_cum <- s24 %>% group_by(month_num) %>%
  summarise(as_tibble(as.list(stats_of(CH4_flux_nmolpm2ps, in_f))), .groups = "drop") %>%
  tidyr::pivot_longer(-month_num, names_to = "stat", values_to = "flux") %>%
  mutate(date = as.Date(sprintf("2024-%02d-01", month_num)), flux = coalesce(flux, 0)) %>%
  group_by(stat) %>% arrange(month_num, .by_group = TRUE) %>%
  mutate(cum = cumsum(to_mmol(flux, date))) %>% ungroup()

# ── Soils along the upland-wetland gradient, campaign by campaign (Fig. 7c) ──
# soils: fences over the whole gradient, as an outlier rule would be applied to the dataset
soil_cum <- a10$soil %>% mutate(in_f = in_fences(flux)) %>%
  group_by(year_month) %>%
  summarise(as_tibble(as.list(stats_of(flux, in_f))), .groups = "drop") %>%
  select(-gmean) %>%    # the positive fluxes of a net sink say nothing about its budget
  left_join(a10$soil_budget_df %>% distinct(year_month, date_mid), by = "year_month") %>%
  tidyr::pivot_longer(-c(year_month, date_mid), names_to = "stat", values_to = "flux") %>%
  group_by(stat) %>% arrange(date_mid, .by_group = TRUE) %>%
  mutate(cum = cumsum(to_mmol(flux, date_mid))) %>% ungroup()

# ── Every dataset: closure-level fluxes, MDF and (stems) the CO2-tracer flag ──
sm <- ld("20_generality.RData")$soil_sigma
ag <- read.csv(file.path(agsoil_dir, "agsoil_closures_deidentified.csv"))
ag <- ag[!duplicated(ag[, c("CH4_flux", "N2O_flux", "CO2_flux")]), ]
sm_p <- merge(sm, a10$soil[, c("unique_id", "CH4_LM_pval")], by.x = "UniqueID", by.y = "unique_id", all.x = TRUE)
st <- a03 %>% filter(t_src == "measured")
sets <- list(
  "Stems, noisier analyzer" = with(st[st$inst_label == "Noisier analyzer", ], list(v = CH4_flux_nmolpm2ps, mdf = CH4_MDF_ref, r2 = CH4_r2, p = lm_p_CH4, bad = qc_co2_tracer %in% TRUE)),
  "Stems, precise analyzer" = with(st[st$inst_label == "Precise analyzer", ], list(v = CH4_flux_nmolpm2ps, mdf = CH4_MDF_ref, r2 = CH4_r2, p = lm_p_CH4, bad = qc_co2_tracer %in% TRUE)),
  "Soils, upland"  = with(sm_p[sm_p$landscape_position %in% c("U", "I"), ], list(v = flux, mdf = mdf_ref, r2 = r2, p = CH4_LM_pval)),
  "Soils, wetland" = with(sm_p[sm_p$landscape_position %in% c("WD", "WS"), ], list(v = flux, mdf = mdf_ref, r2 = r2, p = CH4_LM_pval)),
  # CO2 tracer: CO2 rose clearly in every agricultural closure (>= 3 umol m-2 s-1, linear-fit R2 > 0.9 in the
  # source table), so none is invalid by the stems' criterion
  "Agricultural soil, CH4" = list(v = ag$CH4_flux, mdf = ag$CH4_MDF, r2 = ag$CH4_r2, p = ag$CH4_LM_p, bad = !(ag$CO2_flux > 0)),
  "Agricultural soil, N2O" = list(v = ag$N2O_flux, mdf = ag$N2O_MDF, r2 = ag$N2O_r2, p = ag$N2O_LM_p, bad = !(ag$CO2_flux > 0)))
coastal_file <- file.path(data_dir, "from_coastal_wetland_study", "coastal_closures.csv")
if (file.exists(coastal_file)) {
  cw <- read.csv(coastal_file)
  for (g in c("CH4", "N2O")) sets[[paste("Coastal wetland,", g)]] <- with(cw[cw$gas == g, ], list(v = best.flux, mdf = MDF_ref, r2 = LM.r2, p = lm_p))
}

# each statistic or removal as a percentage of the mean of all closures (signed: a negative
# value means the summary has the opposite sign to the mean)
stat_cmp <- bind_rows(lapply(names(sets), function(nm) {
  s <- sets[[nm]]; ok <- is.finite(s$v); v <- s$v[ok]; m <- mean(v)
  r2 <- s$r2[ok]; pv <- s$p[ok]
  x <- c(stats_of(v)[c("median", "gmean", "trimmed")], below = mean(v[abs(v) >= s$mdf[ok]]),
         zero = mean(ifelse(abs(v) < s$mdf[ok], 0, v)), zero_median = median(ifelse(abs(v) < s$mdf[ok], 0, v)),
         r2 = mean(v[is.finite(r2) & r2 > 0.5]), pval = mean(v[is.finite(pv) & pv < 0.05]))
  if (!is.null(s$bad)) x["invalid"] <- mean(v[!s$bad[ok]])
  if (m < 0) x["gmean"] <- NA          # a net sink: the positive fluxes alone are not its budget
  tibble(dataset = nm, n = length(v), mean = m, proc = names(x), value = unname(x), pct = 100 * value / m,
         pct_removed = c(median = NA, gmean = 100 * mean(v <= 0), trimmed = 100 * mean(!in_fences(v)),
                         below = 100 * mean(abs(v) < s$mdf[ok]), zero = 100 * mean(abs(v) < s$mdf[ok]),
                zero_median = 100 * mean(abs(v) < s$mdf[ok]), r2 = 100 * mean(!(is.finite(r2) & r2 > 0.5)),
                pval = 100 * mean(!(is.finite(pv) & pv < 0.05)), invalid = if (is.null(s$bad)) NA else 100 * mean(s$bad[ok]))[names(x)])
}))

# share of the net total carried by the largest closures (emission-dominated datasets)
conc <- bind_rows(lapply(intersect(c("Stems, noisier analyzer", "Stems, precise analyzer", "Agricultural soil, N2O",
                                     "Coastal wetland, CH4"), names(sets)), function(nm) {
  v <- sort(sets[[nm]]$v[is.finite(sets[[nm]]$v)], decreasing = TRUE)
  tibble(dataset = nm, frac_closures = seq_along(v) / length(v), frac_total = cumsum(v) / sum(v))
}))

# are the stem outliers recurrent hot trees? (noisier analyzer, its own fence)
hot_trees <- local({
  d <- st %>% filter(inst_label == "Noisier analyzer") %>% mutate(tree = paste(PLOT, Tree), v = CH4_flux_nmolpm2ps,
                     hi = v > quantile(v, 0.75) + 1.5 * IQR(v))
  tr <- d %>% group_by(tree) %>% summarise(nhi = sum(hi), .groups = "drop")
  rep <- tr$tree[tr$nhi >= 2]
  list(n_trees = nrow(tr), n_hi = sum(d$hi), n_hot_trees = length(rep),
       share_hi_on_repeat = mean(d$tree[d$hi] %in% rep),
       share_total_repeat = sum(d$v[d$tree %in% rep]) / sum(d$v),
       share_closures_repeat = mean(d$tree %in% rep),
       hi_co2_flagged = sum(d$hi & d$qc_co2_tracer %in% TRUE),
       hi_species = table(d$species_label[d$hi]))
})

# every stem closure with its tree, for drawing the distribution by tree: below its own MDF, or above
# the upper Tukey fence of all stem closures (the high outliers above)
stem_trees <- st %>% filter(inst_label == "Noisier analyzer") %>%
  transmute(tree = paste(PLOT, Tree), species = species_label, flux = CH4_flux_nmolpm2ps,
                               below = abs(flux) < CH4_MDF_ref,
                               high = !in_fences(flux)) %>%
  group_by(tree) %>% mutate(tree_median = median(flux)) %>% ungroup()

# every soil closure with its collar, for the same drawing: below its MDF, or outside the fences of
# the whole gradient
soil_collars <- sm %>% transmute(collar = sub("^[0-9]+_(.*)_[0-9]+$", "\\1", UniqueID), flux, below = abs(flux) < mdf_ref,
                                 part = if_else(landscape_position %in% c("U", "I"), "Upland", "Wetland")) %>%
  mutate(out = !in_fences(flux), lo_f = quantile(flux, 0.25) - 1.5 * IQR(flux),
         hi_f = quantile(flux, 0.75) + 1.5 * IQR(flux)) %>%
  group_by(collar) %>% mutate(collar_median = median(flux)) %>% ungroup()

message("Stems 2024 budget (mmol m-2) by statistic:")
print(stem_cum %>% filter(month_num == 12) %>% select(stat, cum))
message("Soil budget (mmol m-2) by statistic:")
print(soil_cum %>% group_by(stat) %>% slice_max(date_mid, n = 1) %>% select(stat, cum))
print(stat_cmp %>% mutate(across(where(is.numeric), ~ signif(.x, 3))), n = 60)
str(hot_trees)
write.csv(stat_cmp, file.path(output_dir, "25_budget_statistic_table.csv"), row.names = FALSE)
save(stem_cum, soil_cum, stat_cmp, conc, hot_trees, stem_trees, soil_collars, stat_lab, file = file.path(output_dir, "25_budget_statistic.RData"))
message("Saved 25_budget_statistic")
