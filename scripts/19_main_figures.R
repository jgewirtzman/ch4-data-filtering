# ============================================================================
# 19_main_figures.R — Assemble the main-text figures 1, 3, 4, 6, 7, 8, 9 and 10 and
#                     SI Figs S3 and S6 from the plot-ready data saved by the
#                     analysis scripts, in the paper style
#
# Figure 2 (sigma explained) is built by 17_sigma_figures.R, Figure 5 (CO2 as a
# closure tracer) by 15_co2_validity_figure.R. 11_organize_figures.R copies all
# figures into figures/main_text/ and figures/SI/ under their published numbers.
#
# Loads: outputs/03_mdf_results.RData, 07_treatment.RData, 08_transformation.RData,
#        09_detection_rates.RData, 10_seasonal_confounding.RData, 16_pvalue_vs_mdf.RData,
#        18_known_truth.RData, 20_generality.RData, 21_logging_interval.RData,
#        22_transform_inference.RData, 25_budget_statistic.RData, outputs/precision_comparison.csv
# Writes (figures/, .png/.pdf): processing_chain (Fig. 1), precision_detection (Fig. 3),
#         regression_vs_mdf (Fig. 4), r2_seasonal (Fig. 6), treatment_budget (Fig. 7),
#         transformation (Fig. 8), budget_statistic (Fig. 9), generality (Fig. 10), si_logging_interval (Fig. S3),
#         si_transform_demo (Fig. S6)
# ============================================================================

source("scripts/00_setup.R")
W <- 7.5                          # full page width (in)
th <- theme_paper(base_size = 8)
tg <- function(p, l) p + labs(tag = l) + theme(plot.tag = element_text(face = "bold", size = 9))
add_tags <- function(p) p + plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(face = "bold", size = 9))
save_main <- function(p, name, h) {
  ggsave(file.path(fig_dir, paste0(name, ".png")), p, width = W, height = h, dpi = 300)
  ggsave(file.path(fig_dir, paste0(name, ".pdf")), p, width = W, height = h)
  message("Saved ", name)
}
ld <- function(f) { e <- new.env(); load(file.path(output_dir, f), envir = e); e }

a03 <- ld("03_mdf_results.RData")
stem <- a03$df %>% mutate(analyzer = if_else(year == 2025, "Precise analyzer", "Noisier analyzer"))

# ── Figure 3: realized precision and what it does to detection ───────────────
gen <- ld("20_generality.RData")$gen
# coastal wetland closures (24_coastal_wetland.R), when present
cw_file <- file.path(data_dir, "from_coastal_wetland_study", "coastal_closures.csv")
cw <- if (file.exists(cw_file)) read.csv(cw_file) %>% filter(gas %in% c("CH4", "N2O")) else NULL
pc <- read.csv(file.path(output_dir, "precision_comparison.csv"))
camp <- function(src) pc$CH4_precision_ppb[pc$source == src & grepl("^Campaign", pc$method)]
prec <- tibble::tribble(
  ~row,                                        ~gas,   ~spec,          ~sigma,
  "Noisier CH₄ analyzer, unit 1 \u00b7 stems",       "CH4",  PREC_CH4_LGR,   median(stem$sigma_camp_CH4[stem$analyzer == "Noisier analyzer"], na.rm = TRUE),
  "Noisier CH₄ analyzer, unit 1 \u00b7 canopy",      "CH4",  PREC_CH4_LGR,   camp("A-1"),
  "Noisier CH₄ analyzer, unit 2 \u00b7 canopy",      "CH4",  PREC_CH4_LGR,   camp("A-3"),
  "Noisier CH₄ analyzer \u00b7 soils",               "CH4",  PREC_CH4_LGR,   gen$sigma[gen$dataset == "Soils, upland"],
  "Precise CH₄ analyzer \u00b7 stems",               "CH4",  PREC_CH4_7810,  median(stem$sigma_camp_CH4[stem$analyzer == "Precise analyzer"], na.rm = TRUE),
  "Precise CH₄ analyzer \u00b7 agricultural soil",   "CH4",  PREC_CH4_7810,  gen$sigma[gen$dataset == "Agricultural soil, CH4"],
  "N₂O analyzer \u00b7 agricultural soil",           "N2O",  PREC_N2O_7820,  gen$sigma[gen$dataset == "Agricultural soil, N2O"])
if (any(gen$dataset == "Coastal wetland, CH4")) prec <- bind_rows(prec, tibble::tribble(
  ~row,                                        ~gas,   ~spec,          ~sigma,
  "Precise CH₄ analyzer \u00b7 coastal wetland",    "CH4",  PREC_CH4_7810,  gen$sigma[gen$dataset == "Coastal wetland, CH4"],
  "N₂O analyzer \u00b7 coastal wetland",            "N2O",  PREC_N2O_7820,  gen$sigma[gen$dataset == "Coastal wetland, N2O"]))
prec <- prec %>% mutate(ratio = sigma / spec, row = factor(row, levels = rev(row)))
p3a <- ggplot(prec, aes(y = row)) +
  geom_segment(aes(x = spec, xend = sigma, yend = row), colour = AXIS, linewidth = 1.2) +
  geom_point(aes(x = spec), shape = 21, size = 2.2, colour = COL$datasheet, fill = "white", stroke = 0.7) +
  geom_point(aes(x = sigma, colour = ratio >= 1), size = 2.6) +
  scale_colour_manual(values = c(`TRUE` = COL$r2, `FALSE` = COL$mdf), guide = "none") +
  geom_text(aes(x = 12, label = ifelse(ratio >= 1, sprintf("%.1f\u00d7 worse", ratio), sprintf("%.0f\u00d7 better", 1 / ratio))), hjust = 0, size = 2.3, colour = INK2) +
  scale_x_log10(breaks = c(0.03, 0.1, 0.3, 1, 3, 10), labels = c("0.03", "0.1", "0.3", "1", "3", "10")) +
  coord_cartesian(xlim = c(0.04, 11), clip = "off") +
  labs(x = "\u03c3 (ppb); open: datasheet, filled: realized", y = NULL) +
  th + theme(panel.grid.major.y = element_blank(), plot.margin = margin(4, 48, 4, 4))

# Fig. S3: realized sigma against the logging interval
lg <- ld("21_logging_interval.RData")$logsum %>% mutate(rel = sigma / sigma1)
p_s3 <- ggplot(lg %>% filter(dt > 1 | how == "1 s record"), aes(dt, rel)) +
  geom_line(data = data.frame(dt = c(1, 10)) %>% mutate(rel = 1 / sqrt(dt)), aes(dt, rel), colour = MUTED,
            linetype = "dashed", inherit.aes = FALSE) +
  geom_line(data = lg %>% filter(how != "snapshot"), aes(colour = analyzer)) +
  geom_line(data = bind_rows(lg %>% filter(how == "1 s record"), lg %>% filter(how == "snapshot")),
            aes(colour = analyzer), linetype = "22") +
  geom_point(aes(colour = analyzer, shape = how == "snapshot"), size = 1.6) +
  annotate("text", x = 10, y = 1 / sqrt(10), label = "datasheet\nrescaling", hjust = 1.1, vjust = 1.2, size = 2.2, colour = MUTED, lineheight = 0.9) +
  scale_x_log10(breaks = c(1, 2, 5, 10)) + scale_y_log10(breaks = c(0.3, 0.5, 1, 1.5)) +
  scale_colour_manual(values = PAL_ANALYZER, name = NULL) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 1), labels = c("averaged", "one reading per interval"), name = NULL) +
  labs(x = "Logging interval (s)", y = expression(sigma~"relative to 1 s record")) +
  th + theme(legend.position = "bottom", legend.box = "vertical", legend.spacing.y = unit(1, "pt"))

rd <- ld("09_detection_rates.RData")$rate_df %>% filter(dataset == "HF stem") %>%
  mutate(method = as.character(method), method = sub("Manufacturer\n\\(goFlux\\)", "Manufacturer", method),
         family = case_when(grepl("Manufacturer", method) ~ "Datasheet",
                            grepl("Campaign", method) ~ "Campaign \u03c3 \u00d7 z",
                            TRUE ~ "Per-closure \u03c3 \u00d7 3t"),
         level = case_when(grepl("90%", method) ~ "90", grepl("95%", method) ~ "95",
                           grepl("99%", method) ~ "99", TRUE ~ "spec"),
         family = factor(family, levels = c("Datasheet", "Campaign \u03c3 \u00d7 z", "Per-closure \u03c3 \u00d7 3t")),
         instrument = factor(ifelse(grepl("Noisier", instrument), "Noisier", "Precise")))
rd <- rd %>% filter(family != "Per-closure σ × 3t") %>% droplevels()
rd <- rd %>% mutate(instrument = paste("Stems,", tolower(as.character(instrument)), "analyzer"))
if (!is.null(cw)) rd <- bind_rows(rd, bind_rows(lapply(c("CH4", "N2O"), function(g) { d <- cw[cw$gas == g, ]
  lab <- paste("Coastal wetland,", ifelse(g == "CH4", "CH\u2084", "N\u2082O"))
  bind_rows(tibble(instrument = lab, family = "Datasheet", level = "spec", pct_below = 100 * mean(abs(d$best.flux) < d$MDF_ds)),
            tibble(instrument = lab, family = "Campaign \u03c3 \u00d7 z", level = c("90", "95", "99"),
                   pct_below = sapply(c(qnorm(0.95), qnorm(0.975), qnorm(0.995)), function(z)
                     100 * mean(abs(d$best.flux) < d$MDF_ref * z / qnorm(0.975))))) })))
rd <- rd %>% mutate(family = factor(family, levels = c("Datasheet", "Campaign \u03c3 \u00d7 z")),
                    instrument = factor(instrument, levels = unique(instrument)))
rd_ref <- rd %>% filter(level %in% c("spec", "95"))
rd_rng <- rd %>% filter(level != "spec") %>% group_by(instrument, family) %>%
  summarise(lo = min(pct_below), hi = max(pct_below), .groups = "drop")
fam_col <- c("Datasheet" = COL$datasheet, "Campaign \u03c3 \u00d7 z" = COL$mdf)
fam_lab <- c("Datasheet" = "Datasheet \u03c3", "Campaign \u03c3 \u00d7 z" = "Measured \u03c3, z = 1.96")
p3b <- ggplot(rd_ref, aes(family, pct_below, fill = family)) +
  geom_col(width = 0.7) +
  geom_errorbar(data = rd_rng, aes(x = family, ymin = lo, ymax = hi), inherit.aes = FALSE, width = 0.25,
                linewidth = 0.35, colour = INK) +
  geom_text(data = rd_ref %>% left_join(rd_rng, by = c("instrument", "family")) %>% mutate(y = pmax(pct_below, hi, na.rm = TRUE)),
            aes(y = y, label = sprintf("%.0f%%", pct_below)), vjust = -0.7, size = 2.4, colour = INK) +
  facet_wrap(~ instrument, nrow = 1, labeller = label_wrap_gen(18)) +
  scale_fill_manual(values = fam_col, labels = fam_lab, name = "MDF built from") +
  scale_y_continuous(limits = c(0, 70), breaks = c(0, 20, 40, 60), expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "% below MDF") +
  th + theme(axis.text.x = element_blank()) + legend_top()

fig3 <- add_tags(free(p3a) / p3b + plot_layout(heights = c(1, 0.9)))
# (Fig. 3 of earlier drafts; its panels are now part of the merged precision figure)
save_main(p_s3 + theme(legend.position = "right"), "si_logging_interval", 3.2)

# ── Figure 4: regression-based detection vs the MDF ─────────────────────────
a16 <- ld("16_pvalue_vs_mdf.RData")
pv <- a16$all_det %>%
  mutate(snr_mdf = abs(flux) / CH4_MDF_ref, p_val = pmax(p_val, 1e-100),
         class = case_when(detect_p & detect_mdf ~ "Both", detect_p & !detect_mdf ~ "Regression only",
                           !detect_p & detect_mdf ~ "MDF only", TRUE ~ "Neither"),
         class = factor(class, levels = c("Both", "Regression only", "MDF only", "Neither")),
         panel = factor(paste("Stems,", tolower(analyzer)))) %>% filter(dataset == "HF stem") %>%
  filter(is.finite(snr_mdf), snr_mdf > 0)
cnt <- pv %>% count(panel, class) %>% group_by(panel) %>% mutate(pct = 100 * n / sum(n)) %>%
  filter(class %in% c("Regression only", "MDF only"))
# ratios beyond the axis (0.01-300) are drawn at its edge; counts use every closure
p4 <- ggplot(pv %>% mutate(snr_plot = pmin(pmax(snr_mdf, 0.01), 300)), aes(snr_plot, -log10(p_val), colour = class)) +
  annotate("rect", xmin = 0.01, xmax = 1, ymin = -log10(0.05), ymax = Inf, fill = TINT_OK) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") + geom_vline(xintercept = 1, linetype = "dashed") +
  geom_point(shape = 16, size = 0.8, alpha = 0.7) +
  facet_wrap(~ panel, nrow = 1) +
  scale_x_log10(limits = c(0.01, 300), breaks = c(0.01, 0.1, 1, 10, 300), labels = c("≤0.01", "0.1", "1", "10", "≥300")) +
  scale_y_sqrt(breaks = c(0, -log10(0.05), 5, 25, 50, 100),
               labels = c("1", "0.05", "1e-5", "1e-25", "1e-50", "≤1e-100")) +
  scale_colour_manual(values = c("Both" = INK, "Regression only" = COL$reg, "MDF only" = COL$mdf, "Neither" = GREY),
                      name = "Detected by") +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  labs(x = "|flux| / realized MDF", y = "Slope p-value") + th +
  theme(legend.position = "bottom", legend.justification = "center")
save_main(p4, "regression_vs_mdf", 3.4)

# ── Figure 6: filters on R² remove periods of small flux ────────────────────────
a10 <- ld("10_seasonal_confounding.RData")
flux_breaks <- c(-10, -1, -0.1, -0.01, 0, 0.01, 0.1, 1, 10, 100)
flux_labs <- c("-10", "-1", "-0.1", "", "0", "", "0.1", "1", "10", "100")
signed_tr <- scales::pseudo_log_trans(sigma = 0.01, base = 10)   # linear within +/-0.01, log beyond
sx <- scale_x_continuous(trans = signed_tr, breaks = flux_breaks, labels = flux_labs)
PAL_R2FATE <- c("Kept" = INK, "Removed by R² ≤ 0.5" = "#c2c2be")   # kept: solid dark; removed: faint grey crosses
SHP_R2FATE <- c("Kept" = 16, "Removed by R² ≤ 0.5" = 4)
r2_panel <- function(d, ttl, note_x, ylab = expression(CH[4]~R^2), xlab = LAB_FLUX) {
  d <- d %>% mutate(fate = factor(if_else(r2 <= 0.5, "Removed by R² ≤ 0.5", "Kept"), levels = names(PAL_R2FATE))) %>%
    arrange(fate)
  ggplot(d, aes(flux, r2, colour = fate)) +
    geom_vline(xintercept = 0, colour = AXIS) +
    geom_hline(yintercept = 0.5, colour = INK2, linewidth = 0.4) +
    geom_point(aes(shape = fate), size = 0.9, stroke = 0.35) +
    sx + scale_y_continuous(breaks = c(0, 0.5, 1)) +
    scale_colour_manual(values = PAL_R2FATE, guide = "none") + scale_shape_manual(values = SHP_R2FATE, guide = "none") +
    labs(title = ttl, x = xlab, y = ylab) + th
}
p6a <- r2_panel(a10$df_gla %>% transmute(flux = CH4_flux_nmolpm2ps, r2 = CH4_r2) %>% filter(is.finite(r2)),
                "Stems, noisier analyzer", 3)
p6b <- r2_panel(a10$soil %>% transmute(flux, r2) %>% filter(is.finite(r2)),
                "Soil gradient", 3)
p6w <- if (!is.null(cw)) r2_panel(cw %>% filter(gas == "N2O") %>% transmute(flux = best.flux, r2 = LM.r2),
                                  "Coastal wetland, N\u2082O", 3, expression(N[2]*O~R^2),
                                  expression(N[2]*O~flux~(nmol~m^-2~s^-1))) +
  scale_x_continuous(trans = scales::pseudo_log_trans(sigma = 0.001, base = 10),
                     breaks = c(-0.1, -0.01, -0.001, 0, 0.001, 0.01, 0.1), labels = c("-0.1", "-0.01", "", "0", "", "0.01", "0.1")) else NULL

fate <- a10$df_gla %>% filter(is.finite(CH4_r2), is.finite(airt)) %>%
  mutate(fate = factor(if_else(CH4_r2 <= 0.5, "Removed by R² ≤ 0.5", "Kept"), levels = names(PAL_R2FATE))) %>%
  arrange(fate)
bin_means <- fate %>% mutate(bin = cut(airt, seq(-10, 35, 5))) %>% group_by(bin) %>%
  summarise(t = mean(airt), `All closures` = mean(CH4_flux_nmolpm2ps),
            `Kept by R² > 0.5` = mean(CH4_flux_nmolpm2ps[CH4_r2 > 0.5]), n = n(), .groups = "drop") %>%
  filter(n >= 15) %>% tidyr::pivot_longer(c(`All closures`, `Kept by R² > 0.5`), names_to = "what", values_to = "m")
p6c <- ggplot(fate, aes(airt, CH4_flux_nmolpm2ps)) +
  geom_hline(yintercept = 0, colour = AXIS) +
  geom_point(data = fate %>% filter(fate != "Kept"), aes(colour = fate, shape = fate), size = 0.55, stroke = 0.3) +
  geom_point(data = fate %>% filter(fate == "Kept"), aes(colour = fate, shape = fate), size = 0.6, alpha = 0.75) +
  geom_line(data = bin_means, aes(t, m, group = what), colour = "white", linewidth = 2.4) +
  geom_line(data = bin_means, aes(t, m, linetype = what), colour = COL$mdf, linewidth = 1) +
  geom_point(data = bin_means, aes(t, m, fill = what), shape = 21, colour = COL$mdf, size = 2.4, stroke = 0.9) +
  scale_y_continuous(trans = signed_tr, breaks = flux_breaks, labels = flux_labs) +
  scale_colour_manual(values = PAL_R2FATE, breaks = names(PAL_R2FATE), name = NULL, guide = guide_legend(override.aes = list(size = 2, alpha = 1, stroke = 0.7), order = 1)) +
  scale_shape_manual(values = SHP_R2FATE, breaks = names(SHP_R2FATE), name = NULL, guide = guide_legend(order = 1)) +
  scale_linetype_manual(values = c("solid", "22"), name = "Mean flux, 5 °C bins", guide = guide_legend(order = 2)) +
  scale_fill_manual(values = c(COL$mdf, "white"), name = "Mean flux, 5 °C bins", guide = guide_legend(order = 2)) +
  labs(
       x = "Air temperature (°C)", y = LAB_FLUX) +
  th + theme(legend.position = "bottom", legend.justification = "center", legend.box = "vertical",
             legend.box.just = "center", legend.spacing.y = unit(2, "pt"))

# (d) share removed by temperature bin: each dataset (stems, soils) x each filter (R2, MDF)
s_soil <- ld("10_seasonal_confounding.RData")$soil
s_soil$UniqueID <- s_soil$unique_id
soil_mdf <- ld("20_generality.RData")$soil_sigma %>% select(UniqueID, flux_m = flux, mdf_ref) %>%
  inner_join(s_soil %>% select(UniqueID, temp_bin), by = "UniqueID") %>% filter(!is.na(temp_bin)) %>%
  group_by(temp_bin) %>% summarise(pct_flagged = 100 * mean(abs(flux_m) < mdf_ref), .groups = "drop") %>%
  mutate(filter = "Campaign \u03c3 95%")
lab_f <- c("CH4 R² > 0.5" = "R² ≤ 0.5", "Campaign \u03c3 95%" = "Below MDF")
tb <- bind_rows(a10$temp_flagged_df %>% mutate(filter = as.character(filter), ds = "Stems"),
                a10$soil_temp_flagged_df %>% mutate(filter = as.character(filter), ds = "Soils"),
                soil_mdf %>% mutate(ds = "Soils")) %>%
  filter(filter %in% names(lab_f)) %>%
  mutate(filter = factor(lab_f[filter], levels = lab_f), ds = factor(ds, levels = c("Stems", "Soils")),
         bin = as.character(temp_bin))
tb$bin_mid <- sapply(regmatches(tb$bin, gregexpr("-?[0-9.]+", tb$bin)), function(v) mean(as.numeric(v)))
p6d <- ggplot(tb, aes(bin_mid, pct_flagged, colour = filter, linetype = ds, shape = ds, group = interaction(filter, ds))) +
  geom_line(linewidth = 0.7) + geom_point(size = 1.8, fill = "white", stroke = 0.7) +
  scale_colour_manual(values = c("R² ≤ 0.5" = COL$r2, "Below MDF" = COL$mdf), name = "Meeting criterion:") +
  scale_linetype_manual(values = c(Stems = "solid", Soils = "22"), name = NULL) +
  scale_shape_manual(values = c(Stems = 16, Soils = 21), name = NULL) +
  scale_y_continuous(limits = c(0, 100)) +
  labs(x = "Air temperature, 10 °C bins (°C)", y = "Closures meeting criterion (%)") +
  th + legend_in(0.98, 0.98) + theme(legend.box = "vertical", legend.box.just = "left", legend.spacing.y = unit(0, "pt"),
             legend.margin = margin(0, 0, 0, 0)) +
  guides(colour = guide_legend(order = 1, override.aes = list(linetype = "solid", shape = NA)),
         linetype = guide_legend(order = 2, override.aes = list(colour = INK)), shape = guide_legend(order = 2))
# letters are fixed so that they match the manuscript when the coastal panel (c) is absent:
# the public-input run leaves its place empty
top6 <- tg(p6a, "a") | tg(p6b, "b") | (if (is.null(p6w)) plot_spacer() else tg(p6w, "c"))
fig6 <- top6 / (tg(p6c, "d") | tg(p6d, "e")) + plot_layout(heights = c(1, 1.15))
save_main(fig6, "r2_seasonal", 6.4)

# ── Figure 7: below-detection treatment and annual budgets ──────────────────
a07 <- ld("07_treatment.RData")
pcl <- a07$pct_change %>% filter(!grepl("3t", filter)) %>%   # datasheet and campaign-sigma MDFs only
  select(filter, treatment, pct_change_mean, pct_change_median) %>%
  tidyr::pivot_longer(c(pct_change_mean, pct_change_median), names_to = "stat", values_to = "pct") %>%
  mutate(stat = recode(stat, pct_change_mean = "Mean", pct_change_median = "Median"),
         treatment = recode(treatment, "Set below MDF to zero" = "Set to zero", "Remove below MDF" = "Removed"),
         filter = sub("Manufacturer MDF", "Datasheet", sub("Campaign σ ", "Campaign σ\n", filter)))
lv <- c("Datasheet", paste0("Campaign σ\n", c("90%", "95%", "99%")))
pcl$filter <- factor(pcl$filter, levels = intersect(lv, unique(pcl$filter)))
short_lab <- setNames(c("Spec", "90%", "95%", "99%"), lv)
p7a <- ggplot(pcl, aes(filter, pct, fill = treatment)) +
  geom_hline(yintercept = 0, colour = INK, linewidth = 0.3) +
  geom_col(position = position_dodge(0.75), width = 0.7) +
  facet_wrap(~ stat, nrow = 1) +
  scale_fill_manual(values = c("Set to zero" = COL$zero, "Removed" = COL$mdf), name = "Below-MDF closures") +
  scale_x_discrete(labels = short_lab) +
  labs(
       x = "MDF used to flag", y = "Change from keeping all (%)") +
  th + theme(legend.position = "top", legend.justification = "left")
cum <- a10$monthly_budget_df %>% filter(treatment != "Campaign σ 95%\n(keep all, flag)") %>%
  mutate(treatment = recode(as.character(treatment), "No filter (keep all)" = "Keep all (flag)",
                            "R² > 0.5 (remove)" = "Remove R² ≤ 0.5",
                            "Campaign σ 95%\n(remove)" = "Remove below MDF")) %>%
  # treatments are compared over the months with closures under every treatment (no month is
  # filled: after R2 removal January 2024 has none, so January is left out of all three)
  group_by(month_num) %>% filter(all(is.finite(mean_flux))) %>% ungroup() %>%
  # monthly mean (nmol m-2 s-1) x seconds in the month -> mmol m-2, cumulated
  group_by(treatment) %>% arrange(month_num, .by_group = TRUE) %>%
  mutate(cum_flux = cumsum(mean_flux * lubridate::days_in_month(as.Date(sprintf("2024-%02d-01", month_num))) * 86400 / 1e6)) %>%
  ungroup()
fin <- cum %>% filter(month_num == 12)
base <- fin$cum_flux[fin$treatment == "Keep all (flag)"]
fin <- fin %>% mutate(lab = if_else(treatment == "Keep all (flag)", treatment, sprintf("%s (+%d%%)", treatment, as.integer(round(100 * (cum_flux / base - 1) + 1e-9)))))
pal_cum <- c("Keep all (flag)" = COL$keep, "Remove R² ≤ 0.5" = COL$r2, "Remove below MDF" = COL$mdf)
p7b <- ggplot(cum, aes(month_num, cum_flux, colour = treatment)) +
  annotate("rect", xmin = c(0.5, 10.5), xmax = c(4.5, 12.5), ymin = -Inf, ymax = Inf, fill = TINT_NEU) +
  geom_line(linewidth = 0.7) + geom_point(size = 1) +
  scale_colour_manual(values = pal_cum, labels = setNames(fin$lab, fin$treatment), breaks = c("Keep all (flag)", "Remove below MDF", "Remove R² ≤ 0.5"), name = NULL) +
  scale_x_continuous(breaks = 1:12, labels = substr(month.abb, 1, 1)) +
  scale_y_continuous(expand = expansion(mult = c(0.03, 0.12))) +
  guides(colour = guide_legend(ncol = 1)) +
  labs(title = "Stems (2024): net emission",
       x = NULL, y = expression(Cumulative~CH[4]~flux~(mmol~m^-2))) + th + legend_in(0.02, 0.98)
sb <- a10$soil_budget_df %>% mutate(treatment = recode(as.character(treatment), "No filter (keep all)" = "Keep all (flag)",
                                                       "R² > 0.3 (remove)" = "Remove R² ≤ 0.3",
                                                       "R² > 0.5 (remove)" = "Remove R² ≤ 0.5",
                                                       "R² > 0.7 (remove)" = "Remove R² ≤ 0.7"))
sb <- sb %>% filter(treatment %in% c("Keep all (flag)", "Remove R² ≤ 0.5"))
# Remove below the reference MDF (soil MDF from 20), cumulated the same way as in 10:
# the sum of campaign (month) means, months without a campaign counting as zero
s_mdf <- a10$soil %>% inner_join(ld("20_generality.RData")$soil_sigma %>% select(unique_id = UniqueID, mdf_ref),
                                 by = "unique_id") %>% filter(abs(flux) >= mdf_ref) %>%
  group_by(year_month) %>% summarise(mean_flux = mean(flux), .groups = "drop")
s_mdf <- sb %>% filter(treatment == "Keep all (flag)") %>% select(year_month, date_mid) %>%
  left_join(s_mdf, by = "year_month") %>% arrange(date_mid) %>%
  mutate(cum_flux = cumsum(mean_flux), treatment = "Remove below MDF")
sb <- bind_rows(sb, s_mdf) %>% group_by(year_month) %>% filter(all(is.finite(mean_flux))) %>% ungroup() %>%
  group_by(treatment) %>% arrange(date_mid, .by_group = TRUE) %>%
  mutate(cum_flux = cumsum(mean_flux * lubridate::days_in_month(date_mid) * 86400 / 1e6)) %>% ungroup()
# start every curve at zero, half a month before the first campaign
sb <- bind_rows(sb %>% distinct(treatment) %>% mutate(date_mid = min(sb$date_mid) - 15, cum_flux = 0), sb)
sfin <- sb %>% group_by(treatment) %>% slice_max(date_mid, n = 1) %>% ungroup()
sbase <- sfin$cum_flux[sfin$treatment == "Keep all (flag)"]
sfin <- sfin %>% mutate(lab = if_else(treatment == "Keep all (flag)", treatment, sprintf("%s (%+.0f%%)", treatment, 100 * (cum_flux / sbase - 1))))
pal_sb <- c("Keep all (flag)" = COL$keep, "Remove below MDF" = COL$mdf, "Remove R² ≤ 0.5" = COL$r2)
p7c <- ggplot(sb, aes(date_mid, cum_flux, colour = treatment)) +
  geom_hline(yintercept = 0, colour = AXIS) +
  geom_line(linewidth = 0.7) + geom_point(size = 1) +
  scale_colour_manual(values = pal_sb, labels = setNames(sfin$lab, sfin$treatment), name = NULL,
                      breaks = c("Keep all (flag)", "Remove below MDF", "Remove R² ≤ 0.5")) +
  guides(colour = guide_legend(ncol = 1)) +
  scale_x_date(date_breaks = "3 months", date_labels = "%b %y") +
  labs(title = "Soil gradient (2020\u20132021): net uptake", x = NULL,
       y = expression(Cumulative~CH[4]~flux~(mmol~m^-2))) + th + legend_in(0.02, 0.02)
sim <- ld("18_known_truth.RData")$score %>% filter(analyzer == "Noisier analyzer") %>%
  tidyr::pivot_longer(c(mean, budget), names_to = "stat", values_to = "err") %>%
  mutate(stat = recode(stat, mean = "Mean", budget = "Budget"),
         treatment = recode(as.character(treatment), "Set below MDF to zero" = "Set to zero", "Remove below MDF" = "Remove below MDF"),
         treatment = factor(treatment, levels = rev(c("Keep all (flag)", "Set to zero", "Remove below MDF", "Remove R² ≤ 0.5"))))
pal_sim <- c("Keep all (flag)" = COL$keep, "Set to zero" = COL$zero, "Remove below MDF" = COL$mdf, "Remove R² ≤ 0.5" = COL$r2)
p7d <- ggplot(sim, aes(err, treatment, colour = treatment, shape = stat)) +
  geom_vline(xintercept = 0, colour = AXIS) +
  stat_summary(fun = median, fun.min = function(v) quantile(v, .05), fun.max = function(v) quantile(v, .95),
               position = position_dodge(width = 0.55), size = 0.3, linewidth = 0.5) +
  scale_colour_manual(values = pal_sim, guide = "none") +
  scale_shape_manual(values = c(Mean = 16, Budget = 1), name = NULL) +
  labs(title = "Simulated stems", x = "Error vs true value (%)", y = NULL) +
  th + legend_top() + theme(panel.grid.major.y = element_blank())
# (the budgets of p7b and p7c are drawn in Fig. 9 a-b, with the other summaries)
fig7 <- add_tags((p7a | p7d) + plot_layout(widths = c(2.6, 1)))
save_main(p7a, "si_treatment_threshold", 3.0)

# ── Figure 8: transformation and what it does to inference ──────────────────
a08 <- ld("08_transformation.RData")
s_tr <- c("Forest stems" = median(ld("03_mdf_results.RData")$df$CH4_MDF_ref, na.rm = TRUE),
          "Upland\u2013wetland soils" = median(ld("20_generality.RData")$soil_sigma$mdf_ref, na.rm = TRUE))
lv_tr <- c("Raw", "log10", "arcsinh(F / s)")
tl <- a08$transform_long %>% mutate(transform = sub(" \\(.*", "", transform)) %>% filter(transform != "arcsinh")
tl <- bind_rows(tl, tl %>% filter(transform == "Raw") %>% mutate(transform = "arcsinh(F / s)", value = mapply(goFlux::flux.asinh, value, s_tr[as.character(dataset)]))) %>%
  filter(!is.na(value)) %>% mutate(transform = factor(transform, levels = lv_tr), dataset = factor(dataset, levels = names(s_tr)))
pcnt <- a08$panel_counts %>% mutate(transform = sub(" \\(.*", "", transform)) %>% filter(transform != "arcsinh")
pcnt <- bind_rows(pcnt, pcnt %>% filter(transform == "Raw") %>% mutate(transform = "arcsinh(F / s)")) %>%
  mutate(transform = factor(transform, levels = lv_tr), dataset = factor(dataset, levels = names(s_tr)),
         label = if_else(n_dropped > 0, sprintf("%d of %d kept", n_total - n_dropped, n_total), sprintf("all %d kept", n_total)))
p8a <- ggplot(tl, aes(value)) +
  geom_histogram(bins = 40, fill = "#6f6f6c", colour = NA) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
  ggh4x::facet_grid2(dataset ~ transform, scales = "free", independent = "all") +
  geom_text(data = pcnt, aes(x = Inf, y = Inf, label = label), hjust = 1.05, vjust = 1.3, size = 2.2,
            colour = INK2, inherit.aes = FALSE) +
  labs(x = expression(CH[4]~flux~(nmol~m^-2~s^-1)*", as transformed"), y = "Closures") + th + theme(strip.text.y = element_text(angle = 0, hjust = 0))
dm <- tibble(x = seq(-6, 6, by = 0.02)) %>%
  mutate(`arcsinh(F / s)` = asinh(x), `sign log(2|F / s|)` = ifelse(abs(x) > 0.6, sign(x) * log(2 * abs(x)), NA), `y = F / s` = x) %>%
  tidyr::pivot_longer(-x, names_to = "fn", values_to = "y") %>% filter(!is.na(y)) %>%
  mutate(fn = factor(fn, levels = c("arcsinh(F / s)", "sign log(2|F / s|)", "y = F / s")))
p8b <- ggplot(dm, aes(x, y, colour = fn, linetype = fn)) +
  annotate("rect", xmin = -1, xmax = 1, ymin = -Inf, ymax = Inf, fill = TINT_NEU) +
  geom_hline(yintercept = 0, colour = AXIS) + geom_vline(xintercept = 0, colour = AXIS) + geom_line(linewidth = 0.6) +
  scale_colour_manual(values = c(INK, MUTED, MUTED), name = NULL) +
  scale_linetype_manual(values = c("solid", "dashed", "dotted"), name = NULL) +
  coord_equal(xlim = c(-6, 6), ylim = c(-6, 6), expand = FALSE) +
  scale_x_continuous(breaks = seq(-6, 6, 3)) + scale_y_continuous(breaks = seq(-6, 6, 3)) +
  labs(x = "F / s (flux / median MDF)", y = "Transformed value") + th + legend_in(0.98, 0.02) +
  theme(legend.justification = c(1, 0), legend.direction = "vertical", legend.key.width = unit(4, "mm"),
        legend.text = element_text(size = 6), legend.background = element_rect(fill = "white", colour = NA))
a22 <- ld("22_transform_inference.RData")
meth_lv <- c("Raw", "log10, positives", "log10(x + 1)", "arcsinh(x)", "arcsinh(x / MDF)", "Ranks")
meth_lab <- c("Raw" = "Raw", "log10, positives" = "log10, positives only", "log10(x + 1)" = "log10(flux + 1)",
              "arcsinh(x)" = "arcsinh(flux)", "arcsinh(x / MDF)" = "arcsinh(F / s)", "Ranks" = "Ranks")
pw <- a22$grp_sim %>% tidyr::pivot_longer(-method, names_to = "scen", values_to = "pct") %>%
  mutate(scen = factor(recode(scen, false_pos = "No true difference", power_2fold = "Fluxes 2\u00d7 larger",
                              power_uptake_share = "More near-zero and uptake"),
                       levels = c("Fluxes 2\u00d7 larger", "More near-zero and uptake", "No true difference")),
         method = factor(method, levels = rev(meth_lv)))
# one panel per scenario: a test should score low on the first (false positives) and high on the others
pw <- pw %>% mutate(scen = factor(scen, levels = c("No true difference", "Fluxes 2\u00d7 larger", "More near-zero and uptake"),
                                  labels = c("No true difference:\nfalse alarms\n(ideal: 5 %)", "True difference:\nemission in one\ngroup 2\u00d7 larger",
                                             "True difference:\nmore near-zero and\nuptake fluxes\n(mean 33 % lower)")),
                    fillv = if_else(scen == levels(scen)[1], NA_real_, pct))
p8c <- ggplot(pw, aes(scen, method)) +
  geom_tile(aes(fill = fillv), colour = "white", linewidth = 1.2) +
  geom_text(aes(label = sprintf("%.0f %%", pct), colour = !is.na(fillv) & pct > 65), size = 2.6) +
  scale_fill_gradient(low = "#eef2f5", high = "#2f5d7c", limits = c(0, 100), na.value = "#f4f4f2", guide = "none") +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = INK), guide = "none") +
  scale_x_discrete(position = "top") + scale_y_discrete(labels = meth_lab) +
  labs(x = NULL, y = NULL, title = "Simulated comparisons of two groups: % called different (p < 0.05)") +
  th + theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.x.top = element_text(size = 6.5, lineheight = 0.9),
             plot.title = element_text(size = 8, hjust = 0.5))
fig8 <- add_tags(wrap_plots(p8a, wrap_plots(p8b, p8c, nrow = 1, widths = c(1, 1.6)), ncol = 1, heights = c(1, 1)))
save_main(fig8, "transformation", 6.4)

# SI: whether and how much, two kinds of difference (known-truth simulation)
whd <- a22$wh %>% filter(method != "Recommended") %>%
  mutate(method = factor(method, levels = rev(names(meth_lab))),
         scenario = factor(scenario, levels = unique(a22$wh$scenario), labels = c("B: every flux 2\u00d7 larger", "B: more hot spots (mean 2\u00d7, typical flux unchanged)")))
ps_a <- ggplot(whd, aes(detected, method)) +
  geom_vline(xintercept = 5, colour = AXIS, linetype = "dashed") +
  geom_segment(aes(x = 0, xend = detected, yend = method), colour = AXIS, linewidth = 0.8) +
  geom_point(size = 2.4, colour = INK) +
  facet_wrap(~ scenario, nrow = 1) + scale_y_discrete(labels = meth_lab) + scale_x_continuous(limits = c(0, 100)) +
  labs(x = "Is a difference detected? (% of simulations, p < 0.05)", y = NULL) + th
ps_b <- ggplot(whd %>% filter(!is.na(ratio_med)), aes(ratio_med, method)) +
  geom_vline(xintercept = 2, colour = COL$r2, linewidth = 0.5) +
  geom_errorbar(aes(xmin = ratio_lo, xmax = ratio_hi), orientation = "y", width = 0, colour = INK2) +
  geom_point(size = 2.4, colour = INK) +
  facet_wrap(~ scenario, nrow = 1) +
  scale_y_discrete(labels = c(meth_lab, "Raw" = "Raw (ratio of means)")) +
  scale_x_log10(breaks = c(0.5, 1, 2, 4)) +
  labs(x = "How much larger is B reported to be? (true: 2\u00d7)", y = NULL) + th
sfig_tr <- add_tags(wrap_plots(ps_a, ps_b, ncol = 1))
ggsave(file.path(fig_dir, "si_transform_demo.png"), sfig_tr, width = W, height = 5.2, dpi = 300)
ggsave(file.path(fig_dir, "si_transform_demo.pdf"), sfig_tr, width = W, height = 5.2)

# ── Figure 1: the processing chain (conceptual) ─────────────────────────────
steps <- c("Precision", "Detection", "Validity", "Filters", "Treatment", "Transformation", "Budgets")
rows  <- c("Common practice", "Problem", "Recommended", "Shown in")
txt <- rbind(
  c("Take \u03c3 from the datasheet",
    "Slope p-value alone, or an MDF with unstated z and t",
    "Reject closures on R\u00b2, p or CO\u2082 alone",
    "Drop closures with low R\u00b2 or p \u2265 0.05",
    "Remove fluxes below the MDF, or set them to zero",
    "Log-transform",
    "Scale from the median or a back-transformed mean; remove outliers"),
  c("Real \u03c3 can be several times worse or better, which changes how many fluxes count as detected",
    "p-values are over-confident when readings are correlated; neither sees leaks or drift",
    "Valid small fluxes are rejected; leaks still pass",
    "Drops small fluxes, and with them whole seasons and places",
    "Removing inflates means and budgets; zeroing skews medians",
    "Drops every zero and negative flux",
    "Can change a budget's size and sign; removes real hot spots"),
  c("Estimate \u03c3 from your own closures",
    "MDF = 1.96 \u03c3 / t \u00d7 flux term with your \u03c3; also report the regression interval",
    "Remove only closures with evidence of failure, not for their flux size; review suspect traces",
    "Do not filter on R\u00b2 or p; report what any filter removes",
    "Keep every measured value, flagged as emission, uptake or below MDF",
    "Test on arcsinh(F / s), one s for all groups, or ranks; sizes from raw values",
    "Scale from the mean of every valid closure; no outlier removal by size"),
  c("Fig 2", "Fig 3", "Fig 4", "Fig 5", "Fig 6", "Fig 7", "Fig 8"))
cells <- expand.grid(row = seq_along(rows), col = seq_along(steps)) %>%
  mutate(label = mapply(function(r, c) txt[r, c], row, col),
         label = vapply(label, function(z) paste(strwrap(z, 19), collapse = "\n"), ""),
         x = col, y = -row)
fill_row <- c(TINT_NEU, "#f7ebe4", "#e8f1ec", "white")
hdr_col <- rep("#4f6576", 7)
arrow_df <- data.frame(x = seq_along(steps))
p1 <- ggplot() +
  geom_tile(data = cells, aes(x, y, fill = factor(row)), width = 0.96, height = 0.94, colour = NA) +
  geom_text(data = cells, aes(x, y, label = label, fontface = ifelse(row == 4, "bold", "plain")),
            size = 2.15, lineheight = 0.92, colour = INK, family = PAPER_FONT) +
  geom_polygon(data = do.call(rbind, lapply(seq_along(steps), function(i)
    data.frame(i = i, x = i + c(-0.48, 0.40, 0.50, 0.40, -0.48, -0.40), y = -0.3 + c(-0.2, -0.2, 0, 0.2, 0.2, 0)))),
    aes(x, y, group = i, fill = factor(i + 10)), colour = NA) +
  geom_text(data = data.frame(x = seq_along(steps), y = -0.3, l = steps),
            aes(x, y, label = l), colour = "white", fontface = "bold", size = 2.15, family = PAPER_FONT) +
  geom_text(data = data.frame(y = -seq_along(rows), l = rows), aes(0.42, y, label = l), hjust = 1,
            size = 2.3, colour = INK2, fontface = "bold", family = PAPER_FONT) +
  annotate("text", x = 4, y = -4.55, label = "Recommendations in brief: Box 1",
           size = 2.3, colour = INK2, fontface = "italic", family = PAPER_FONT) +
  scale_fill_manual(values = c(setNames(fill_row, 1:4), setNames(hdr_col, 11:17)), guide = "none") +
  coord_cartesian(xlim = c(-0.55, 7.5), ylim = c(-4.7, -0.1), clip = "off") +
  theme_void(base_family = PAPER_FONT) + theme(plot.margin = margin(4, 4, 4, 4))
save_main(p1, "processing_chain", 4.3)

# ── Figure 10: the same chain across gases, surfaces and analyzers ───────────
gen <- ld("20_generality.RData")$gen
sw <- ld("18_known_truth.RData")$sweep_sum
short <- c("Stems, noisier analyzer" = "Stems, noisier analyzer", "Stems, precise analyzer" = "Stems, precise analyzer",
           "Soils, upland" = "Upland soils", "Soils, wetland" = "Wetland soils",
           "Agricultural soil, CH4" = "Agricultural soil, CH₄", "Agricultural soil, N2O" = "Agricultural soil, N₂O",
           "Coastal wetland, CH4" = "Coastal wetland, CH₄", "Coastal wetland, N2O" = "Coastal wetland, N₂O")
pal_rm <- c("Remove below MDF" = COL$mdf, "Remove p ≥ 0.05" = COL$reg, "Remove R² ≤ 0.5" = COL$r2)
pts <- gen %>% transmute(lab = short[dataset], x = pct_below_ref,
                         `Mean|Remove below MDF` = mean_change_mdf, `Mean|Remove R² ≤ 0.5` = mean_change_r2,
                         `Median|Remove below MDF` = median_change_mdf, `Median|Remove R² ≤ 0.5` = median_change_r2) %>%
  tidyr::pivot_longer(-c(lab, x), names_to = "k", values_to = "chg") %>%
  tidyr::separate(k, c("stat", "filt"), sep = "\\|") %>% filter(is.finite(chg))
# the slope p-value rule (from 25_budget_statistic.R), for the mean
pv <- ld("25_budget_statistic.RData")$stat_cmp %>% filter(proc == "pval") %>%
  transmute(lab = short[dataset], stat = "Mean", filt = "Remove p ≥ 0.05", chg = abs(pct) - 100) %>%
  left_join(gen %>% transmute(lab = short[dataset], x = pct_below_ref), by = "lab")
pts <- bind_rows(pts, pv) %>% mutate(filt = factor(filt, levels = names(pal_rm)))
labs_df <- pts %>% filter(stat == "Mean") %>% group_by(stat, lab, x) %>% summarise(y = max(chg), .groups = "drop")
swl <- sw %>% select(pct_below, mean_change_mdf, mean_change_r2, median_change_mdf, median_change_r2) %>%
  tidyr::pivot_longer(-pct_below, names_to = "k", values_to = "chg") %>%
  mutate(stat = ifelse(grepl("^mean", k), "Mean", "Median"),
         filt = ifelse(grepl("mdf$", k), "Remove below MDF", "Remove R² ≤ 0.5"), x = pct_below)
p9 <- ggplot() +
  geom_hline(yintercept = 0, colour = AXIS) +
  geom_line(data = swl, aes(x, chg, colour = filt), linewidth = 0.6, alpha = 0.45) +
  geom_segment(data = pts %>% group_by(stat, lab, x) %>% summarise(lo = min(chg), hi = max(chg), .groups = "drop"),
               aes(x = x, xend = x, y = lo, yend = hi), colour = AXIS, linewidth = 0.4) +
  geom_point(data = pts, aes(x, chg, colour = filt), size = 2) +
  ggrepel::geom_text_repel(data = labs_df %>% filter(x >= 10), aes(x, y, label = lab), colour = INK2, size = 2.3, family = PAPER_FONT,
                           nudge_y = 0.25, direction = "both", min.segment.length = 0.2, seed = 3, box.padding = 0.6,
                           point.padding = 0.5, force = 4, max.overlaps = Inf) +
  ggrepel::geom_text_repel(data = labs_df %>% filter(x < 10), aes(x, y, label = lab), colour = INK2, size = 2.3, family = PAPER_FONT,
                           nudge_x = 14, nudge_y = 0.15, direction = "y", hjust = 0, min.segment.length = 0, segment.colour = AXIS,
                           seed = 3, box.padding = 0.35, max.overlaps = Inf) +
  facet_wrap(~ stat, nrow = 1) +
  scale_colour_manual(values = pal_rm, name = NULL,
                      guide = guide_legend(override.aes = list(linetype = "solid", alpha = 1, linewidth = 0.6))) +
  scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 5, base = 10), breaks = c(-30, 0, 10, 30, 100, 300)) +
  scale_x_continuous(limits = c(0, 80)) +
  labs(x = "Closures below the MDF (%)", y = "Increase over keeping every closure (%)") +
  th + legend_top()
save_main(p9, "si_generality", 3.6)


# ── Figure 9: the statistic a budget scales from, and removing large fluxes ──
a25 <- ld("25_budget_statistic.RData")
COL_STAT <- c(mean = INK, below = COL$mdf, r2 = COL$r2, trimmed = COL$outlier, median = COL$median, gmean = COL$gmean)
lab_stat <- c(mean = "Mean, all closures", below = "Mean, below MDF removed", r2 = "Mean, R² ≤ 0.5 removed",
              trimmed = "Mean, outliers removed", median = "Median", gmean = "Geometric mean, > 0")
# removal budgets as in Fig. 7 (cum, sb), the other summaries from 25; one data frame per dataset
trt_key <- c("Keep all (flag)" = "mean", "Remove below MDF" = "below", "Remove R² ≤ 0.5" = "r2")
stem_all <- bind_rows(cum %>% transmute(month_num, cum = cum_flux, stat = trt_key[as.character(treatment)]),
                      a25$stem_cum %>% filter(stat != "mean") %>% select(month_num, cum, stat))
# Fig. 8 compares statistics of all closures, available in every month: its mean comes from 25
stem_all8 <- a25$stem_cum %>% select(month_num, cum, stat)
soil_all <- bind_rows(sb %>% transmute(date_mid, cum = cum_flux, stat = trt_key[as.character(treatment)]),
                      a25$soil_cum %>% filter(stat != "mean") %>% select(date_mid, cum, stat),
                      tidyr::expand_grid(date_mid = min(sb$date_mid), cum = 0, stat = c("median", "trimmed")))
cum_panel <- function(d, xvar, title, legend = TRUE) {
  fin <- d %>% group_by(stat) %>% slice_max(.data[[xvar]], n = 1) %>% ungroup()
  ref <- fin$cum[fin$stat == "mean"]
  # each line's final budget relative to the mean, as a label at its end
  fin <- fin %>% mutate(lab = if_else(stat == "mean", "100%", sprintf("%.0f%%", 100 * cum / ref)))
  ggplot(d, aes(.data[[xvar]], cum, colour = stat)) +
    geom_hline(yintercept = 0, colour = AXIS) +
    geom_line(linewidth = 0.7) + geom_point(size = 1) +
    ggrepel::geom_text_repel(data = fin, aes(label = lab), size = 2.4, hjust = 0, nudge_x = 0.08 * as.numeric(diff(range(d[[xvar]]))),
                             xlim = c(as.numeric(max(d[[xvar]])) + 0.05 * as.numeric(diff(range(d[[xvar]]))), NA), direction = "y", segment.colour = AXIS, min.segment.length = 0, seed = 1, show.legend = FALSE,
                             family = PAPER_FONT) +
    scale_colour_manual(values = COL_STAT, labels = lab_stat[intersect(names(COL_STAT), unique(d$stat))],
                        breaks = intersect(names(COL_STAT), unique(d$stat)), name = NULL,
                        guide = if (legend) "legend" else "none") +
    coord_cartesian(clip = "off") +
    labs(title = title, x = NULL, y = expression(Cumulative~CH[4]~flux~(mmol~m^-2))) + th +
    theme(legend.position = "bottom", plot.margin = margin(4, 30, 4, 14))
}
p10a <- cum_panel(stem_all8, "month_num", "Stems, January\u2013December 2024") +
  scale_x_continuous(breaks = 1:12, labels = substr(month.abb, 1, 1), expand = expansion(mult = c(0.02, 0.16)))
soil_all8 <- bind_rows(a25$soil_cum %>% select(date_mid, cum, stat),
                       tidyr::expand_grid(date_mid = min(a25$soil_cum$date_mid) - 15, cum = 0, stat = c("mean", "median", "trimmed")))
p10b <- cum_panel(soil_all8, "date_mid", "Soils, 2020–21", legend = FALSE) +
  scale_x_date(date_breaks = "4 months", date_labels = "%b %y", expand = expansion(mult = c(0.02, 0.18)))
short_d <- c("Stems, noisier analyzer" = "Stems, noisier", "Stems, precise analyzer" = "Stems, precise",
             "Soils, upland" = "Upland soils", "Soils, wetland" = "Wetland soils",
             "Agricultural soil, CH4" = "Agric. soil, CH₄", "Agricultural soil, N2O" = "Agric. soil, N₂O",
             "Coastal wetland, CH4" = "Coastal, CH₄", "Coastal wetland, N2O" = "Coastal, N₂O")
pal_conc <- c("Stems, noisier analyzer" = INK, "Stems, precise analyzer" = COL$B,
              "Agricultural soil, N2O" = COL$mdf, "Coastal wetland, CH4" = COL$reg)
p10conc <- ggplot(a25$conc, aes(100 * frac_closures, 100 * frac_total, colour = dataset)) +
  geom_hline(yintercept = c(0, 100), colour = AXIS) +
  geom_line(data = data.frame(x = c(0.5, 100)), aes(x, x), colour = MUTED, linetype = "dashed", inherit.aes = FALSE) +
  annotate("text", x = 40, y = 30, label = "equal shares", colour = MUTED, size = 2.2, angle = 12) +
  geom_line(linewidth = 0.7) +
  scale_colour_manual(values = pal_conc, labels = short_d[names(pal_conc)], name = NULL) +
  scale_x_log10(breaks = c(1, 3, 10, 30, 100), labels = c("1", "3", "10", "30", "100")) +
  coord_cartesian(xlim = c(0.5, 100), ylim = c(0, 165)) + scale_y_continuous(breaks = seq(0, 150, 50)) +
  guides(colour = guide_legend(ncol = 1)) +
  labs(title = NULL,
       x = "Largest closures (%, log scale)", y = "Net emission (% of total)") + th + legend_in(0.02, 0.98) +
  theme(legend.justification = c(0, 1))
proc_lab <- c(trimmed = "Remove outliers\n(beyond 1.5 IQR)", median = "Median", gmean = "Geometric mean,\nfluxes > 0",
              invalid = "Remove invalid only\n(CO₂ not rising)")
sc <- a25$stat_cmp %>% filter(is.finite(pct), proc %in% c("trimmed", "median", "gmean")) %>%
  mutate(lab = factor(short_d[dataset], levels = rev(short_d)),
         proc = factor(proc, levels = names(proc_lab)),
         chg = pct - 100,
         gap = 14,                                    # arrowhead stops short of the point
         arrow_end = pct - sign(chg) * gap, moved = abs(chg) > gap + 4)
# one lollipop per dataset and summary: an arrow from the sample mean (100 %) to a black point at the
# summary, red where it is larger than the mean and blue where smaller, darker the further it moves
p10d <- ggplot(sc, aes(y = lab)) +
  geom_vline(xintercept = 0, colour = AXIS, linewidth = 0.3) +
  geom_vline(xintercept = 100, colour = INK, linewidth = 0.4) +
  geom_segment(data = filter(sc, moved), aes(x = 100, xend = arrow_end, yend = lab, colour = proc), linewidth = 0.8,
               arrow = arrow(length = unit(1.5, "mm"), type = "closed")) +
  geom_point(aes(x = pct), colour = INK, size = 1.4) +
  geom_text(aes(x = pct, label = if_else(round(chg) == 0, "0 %", sprintf("%+.0f %%", chg)), hjust = if_else(chg < -0.5, 1.35, -0.35)),
            size = 2.2, colour = INK2) +
  facet_wrap(~ proc, nrow = 1, labeller = as_labeller(proc_lab)) +
  coord_cartesian(clip = "off") +
  # one colour per choice, as in the budget lines of (a, b)
  scale_colour_manual(values = c(trimmed = COL_STAT[["trimmed"]], median = COL_STAT[["median"]],
                                 gmean = COL_STAT[["gmean"]], invalid = INK2), guide = "none") +
  scale_y_discrete(limits = rev(unname(short_d))) +
  scale_x_continuous(breaks = c(0, 100, 200), expand = expansion(mult = c(0.3, 0.25))) +
  labs(x = "% of the sample mean of all closures", y = NULL) +
  th + theme(panel.grid.major.y = element_line(colour = TINT_NEU), strip.text = element_text(size = 6.5, vjust = 0),
             legend.position = "right", legend.key.height = unit(5, "mm"), legend.key.width = unit(2.5, "mm"),
             axis.text.x = element_text(size = 6))
# (c, d) every closure by tree (stems, noisier analyzer) or collar (soils), ordered by median flux:
# grey below its MDF, purple outside the outlier fences; dashed, the fences
cls_lev <- c("Below its MDF", "Other", "Outside outlier fences")
cls_col <- c("Below its MDF" = GREY_D, "Other" = INK, "Outside outlier fences" = COL_STAT[["trimmed"]])
dist_panel <- function(d, unit, fences, title, xlab, legend = TRUE) {
  ggplot(d, aes(reorder(.data[[unit]], med), flux, colour = cls)) +
    geom_hline(yintercept = 0, colour = AXIS) +
    geom_point(position = position_jitter(width = 0.25, height = 0, seed = 1), size = 0.5, alpha = 0.8) +
    geom_segment(data = fences, aes(x = x0, xend = x1, y = y, yend = y), colour = COL_STAT[["trimmed"]],
                 linetype = "dashed", linewidth = 0.4, inherit.aes = FALSE) +
    scale_colour_manual(values = cls_col, name = NULL, guide = if (legend) "legend" else "none") +
    scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 0.05, base = 10), breaks = c(-10, -1, 0, 1, 10, 50)) +
    guides(colour = guide_legend(override.aes = list(size = 1.8, alpha = 1), nrow = 1)) +
    labs(title = title, x = xlab, y = expression(CH[4]~flux~(nmol~m^-2~s^-1))) +
    th + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid.major.x = element_blank())
}
# soils: one panel per part of the gradient, each with its own fences

st_d <- a25$stem_trees %>% mutate(med = tree_median, cls = factor(case_when(high ~ cls_lev[3], below ~ cls_lev[1], TRUE ~ cls_lev[2]), levels = cls_lev))
st_f <- with(a25$stem_trees, tibble(x0 = 0.5, x1 = n_distinct(tree) + 0.5,
                                    y = c(quantile(flux, 0.25) - 1.5 * IQR(flux), quantile(flux, 0.75) + 1.5 * IQR(flux))))
p10c <- dist_panel(st_d, "tree", st_f, "Stems", "Trees, sorted by median flux")
# soils: collars grouped by part of the gradient (upland first), each part with its own fences
so_d <- a25$soil_collars %>% mutate(med = collar_median + if_else(part == "Wetland", 1e3, 0),
                                    cls = factor(case_when(out ~ cls_lev[3], below ~ cls_lev[1], TRUE ~ cls_lev[2]), levels = cls_lev))
n_up <- n_distinct(so_d$collar[so_d$part == "Upland"]); n_all <- n_distinct(so_d$collar)
so_f <- a25$soil_collars %>% distinct(part, lo_f, hi_f) %>% tidyr::pivot_longer(c(lo_f, hi_f), values_to = "y") %>%
  mutate(x0 = if_else(part == "Upland", 0.5, n_up + 0.5), x1 = if_else(part == "Upland", n_up + 0.5, n_all + 0.5))
so_d$med <- so_d$med - if_else(so_d$part == "Wetland", 1e3, 0)    # sort every collar by its own median
so_f2 <- a25$soil_collars %>% distinct(lo_f, hi_f) %>% tidyr::pivot_longer(everything(), values_to = "y") %>%
  mutate(x0 = 0.5, x1 = n_distinct(so_d$collar) + 0.5)
p10d_dist <- dist_panel(so_d, "collar", so_f2, "Soils", "Collars, sorted by median flux", legend = FALSE) +
  annotate("text", x = c(1, n_distinct(so_d$collar)), y = 60, label = c("upland", "wetland"), hjust = c(0, 1),
           size = 2.3, colour = INK2)
row_budget <- (tg(p10a, "c") | tg(p10b, "d")) + plot_layout(guides = "collect") & theme(legend.position = "bottom", legend.justification = "center")
row_dist <- ((tg(p10c, "a") | tg(p10d_dist, "b")) + plot_layout(widths = c(1, 1))) + plot_layout(guides = "collect") & theme(legend.position = "bottom", legend.justification = "center")
fig_bs <- wrap_elements(full = row_dist) / wrap_elements(full = row_budget) / wrap_elements(full = tg(p10d, "e")) +
  plot_layout(heights = c(1, 1, 1.15))
save_main(fig_bs, "budget_statistic", 8.2)

# Fig. S7: share of net emission carried by the largest closures (moved from Fig. 9)
ggsave(file.path(fig_dir, "si_concentration.png"), p10conc, width = 3.6, height = 3.2, dpi = 300)
ggsave(file.path(fig_dir, "si_concentration.pdf"), p10conc, width = 3.6, height = 3.2)


# ── Proposed merged figures (previews) ───────────────────────────────────────
# Precision and detection: how sigma is estimated (17_sigma_figures.R panels) and what it does
sp <- readRDS(file.path(output_dir, "17_sigma_panels.rds"))
tx <- theme(axis.text.x = element_text())
pm_a <- sp$a + labs(x = "Time in closure (s)") + tx
pm_b <- sp$b + labs(x = "Time in closure (s)") + tx
pm_c <- sp$c
# (d) as Fig. 3a, with the analyzer of panels a-c marked: its realized sigma is the median of every
# closure's sigma over each field day
hl <- "Precise CH₄ analyzer \u00b7 stems"
p3a_hl <- ggplot(prec, aes(y = row)) +
  geom_segment(aes(x = spec, xend = sigma, yend = row), colour = AXIS, linewidth = 1.2) +
  geom_point(aes(x = spec, shape = "Datasheet"), size = 2.2, colour = COL$datasheet, fill = "white", stroke = 0.7) +
  geom_point(aes(x = sigma, fill = "Realized"), shape = 21, size = 2.6, colour = "white", stroke = 0.2) +
  geom_point(data = filter(prec, row == hl), aes(x = sigma, colour = "Analyzer of panels a–c"), shape = 21, size = 4.6,
             fill = NA, stroke = 0.6) +
  geom_text(aes(x = 12, label = ifelse(ratio >= 1, sprintf("%.1f\u00d7 worse", ratio), sprintf("%.0f\u00d7 better", 1 / ratio))),
            hjust = 0, size = 2.3, colour = INK2) +
  scale_shape_manual(values = c(Datasheet = 21), name = NULL) +
  scale_fill_manual(values = c("Realized" = INK), name = NULL) +
  scale_colour_manual(values = c("Analyzer of panels a–c" = INK), name = NULL) +
  guides(shape = guide_legend(order = 1), fill = guide_legend(order = 2, override.aes = list(size = 2.6)),
         colour = guide_legend(order = 3)) +
  scale_x_log10(breaks = c(0.03, 0.1, 0.3, 1, 3, 10), labels = c("0.03", "0.1", "0.3", "1", "3", "10")) +
  coord_cartesian(xlim = c(0.04, 11), clip = "off") +
  labs(x = "Precision, \u03c3 (ppb)", y = NULL) +
  th + theme(panel.grid.major.y = element_blank(), plot.margin = margin(4, 48, 4, 4), legend.position = "bottom",
             legend.direction = "vertical", legend.box = "horizontal")
# rows aligned separately (wrap_elements), so the long analyzer names of (d) do not push (a) right;
# panel letters set by hand
p3b_pm <- p3b + facet_wrap(~ instrument, nrow = 1, labeller = as_labeller(function(x)
  sub("Stems, (\\w+) analyzer", "Stems,\n\\1", sub("Coastal wetland, ", "Coastal,\n", x)))) +
  scale_fill_manual(values = fam_col, labels = c("Datasheet" = "Datasheet σ", "Campaign σ × z" = "Realized σ, z = 1.96"),
                    name = NULL) +
  theme(legend.position = "bottom", legend.direction = "vertical", strip.text = element_text(size = 6.5))
fig_pm <- wrap_elements(full = tg(pm_a, "a") | tg(pm_b, "b") | tg(pm_c, "c")) /
  wrap_elements(full = (tg(p3a_hl, "d") | tg(p3b_pm, "e")) + plot_layout(widths = c(1.35, 1))) +
  plot_layout(heights = c(1, 1.15))
save_main(fig_pm, "precision_detection", 5.8)

# Treatment: simulation against the truth (arrows from the truth), group comparisons, and the bias
# against the share below the MDF across datasets
kt <- ld("18_known_truth.RData")
trt_lab <- c("Keep all (flag)" = "Keep all", "Set below MDF to zero" = "Set to zero",
             "Remove below MDF" = "Remove below MDF", "Remove R² ≤ 0.5" = "Remove R² ≤ 0.5")
sim6 <- kt$score %>% filter(analyzer == "Noisier analyzer") %>%
  tidyr::pivot_longer(c(mean, median, budget), names_to = "stat", values_to = "err") %>%
  group_by(treatment, stat) %>%
  summarise(med = median(err), lo = quantile(err, 0.05), hi = quantile(err, 0.95), .groups = "drop") %>%
  mutate(trt = factor(trt_lab[as.character(treatment)], levels = rev(trt_lab)),
         stat = factor(stat, levels = c("mean", "median", "budget"), labels = c("Mean", "Median", "Annual budget")),
         ) %>% mutate(gap = 18, moved = abs(med) > 2 * gap, arrow_end = med - sign(med) * gap,
               val = if_else(abs(med) < 1, "0", sprintf("%+.0f", med)))
# one colour per processing choice, the same in the arrows and in the budget lines
COL_PROC <- c("Set below MDF to zero" = COL$zero, "Remove below MDF" = COL$mdf, "Remove p ≥ 0.05" = COL$reg,
              "Remove R² ≤ 0.5" = COL$r2, "Keep all" = INK)
proc_cols <- scale_colour_manual(values = COL_PROC, guide = "none")
cols5 <- c("Keep all", "Set below MDF to zero", "Remove below MDF", "Remove p ≥ 0.05", "Remove R² ≤ 0.5")
arrow_panel <- function(d, xlab, xref = 0, xlim = NULL, breaks = waiver()) {
  ggplot(d, aes(y = row)) +
    geom_vline(xintercept = xref, colour = INK, linewidth = 0.4) +
    { if ("lo" %in% names(d)) geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0, colour = TINT_NEU, linewidth = 2.6) } +
    geom_segment(data = filter(d, moved), aes(x = xref, xend = arrow_end, yend = row, colour = proc), linewidth = 0.8,
                 arrow = arrow(length = unit(1.5, "mm"), type = "closed")) +
    geom_point(aes(x = med), colour = INK, size = 1.4) +
    geom_text(aes(x = med, label = val, hjust = if_else(med < xref - 0.5, 1.35, -0.35)), size = 2.2, colour = INK2) +
    facet_wrap(~ proc, nrow = 1, drop = FALSE, labeller = as_labeller(c("Keep all" = "Keep all\n ",
      "Set below MDF to zero" = "Set below\nMDF to zero", "Remove below MDF" = "Remove\nbelow MDF",
      "Remove p ≥ 0.05" = "Remove\np ≥ 0.05", "Remove R² ≤ 0.5" = "Remove\nR² ≤ 0.5"))) + proc_cols +
    coord_cartesian(clip = "off") +
    scale_x_continuous(limits = xlim, breaks = breaks) +
    labs(x = xlab, y = NULL) + th + theme(panel.grid.major.y = element_line(colour = TINT_NEU))
}
# (a) simulation: rows are the statistics, columns the choices, as in (b)
trt_proc <- c("Keep all (flag)" = "Keep all", "Set below MDF to zero" = "Set below MDF to zero",
              "Remove below MDF" = "Remove below MDF", "Remove p ≥ 0.05" = "Remove p ≥ 0.05", "Remove R² ≤ 0.5" = "Remove R² ≤ 0.5")
sim6 <- sim6 %>% mutate(proc = factor(trt_proc[as.character(treatment)], levels = cols5),
                        row = factor(paste("Simulated", tolower(as.character(stat))),
                                     levels = rev(c("Simulated mean", "Simulated median", "Simulated annual budget"))),
                        gap = 11, moved = abs(med) > 20, arrow_end = med - sign(med) * gap,
                        val = if_else(round(med) == 0, "0 %", sprintf("%+.0f %%", med)))
p6a <- arrow_panel(sim6, "Simulated stems: difference from the true value (%)", xlim = c(-250, 650), breaks = c(0, 200, 400))
# (b) every dataset, each rule as a percentage of keeping every closure; rows ordered by the share below the MDF
gen_b <- a25$stat_cmp %>% filter(proc %in% c("zero", "below", "pval", "r2"), is.finite(pct)) %>%
  left_join(gen %>% select(dataset, pct_below_ref), by = "dataset") %>%
  mutate(row = sprintf("%s\n%.0f%% below MDF", short_d[dataset], pct_below_ref),
         proc = factor(c(zero = "Set below MDF to zero", below = "Remove below MDF", pval = "Remove p ≥ 0.05", r2 = "Remove R² ≤ 0.5")[proc],
                       levels = cols5),
         med = pct, moved = abs(med - 100) > 20, arrow_end = med - sign(med - 100) * 11,
         val = if_else(round(med - 100) == 0, "0 %", sprintf("%+.0f %%", med - 100)))
gen_b <- bind_rows(gen_b, gen_b %>% distinct(row, pct_below_ref) %>%
                     mutate(proc = factor("Keep all", levels = cols5), med = 100, moved = FALSE, val = ""))
p6b <- arrow_panel(gen_b, "Mean as % of the mean of all closures",
                   xref = 100, xlim = c(40, 480), breaks = c(100, 200, 300, 400)) +
  scale_y_discrete(limits = gen_b %>% distinct(row, pct_below_ref) %>% arrange(pct_below_ref) %>% pull(row))
# (c, d) field budgets, removing small fluxes
p6c <- cum_panel(stem_all %>% filter(stat %in% c("mean", "below", "r2")), "month_num", "Stems, February\u2013December 2024") +
  scale_x_continuous(breaks = 1:12, labels = substr(month.abb, 1, 1), expand = expansion(mult = c(0.02, 0.16)))
p6d <- cum_panel(soil_all %>% filter(stat %in% c("mean", "below", "r2")), "date_mid", "Soils, 2020–21", legend = FALSE) +
  scale_x_date(date_breaks = "4 months", date_labels = "%b %y", expand = expansion(mult = c(0.02, 0.18)))
fig_tr <- wrap_elements(full = (tg(p6c, "a") | tg(p6d, "b")) + plot_layout(guides = "collect") & theme(legend.position = "bottom", legend.justification = "center")) /
  wrap_elements(full = tg(p6b, "c")) / wrap_elements(full = tg(p6a, "d")) +
  plot_layout(heights = c(1, 1.6, 0.6))
save_main(fig_tr, "treatment_budget", 8.4)
