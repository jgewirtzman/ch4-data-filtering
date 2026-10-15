# ============================================================================
# 19_main_figures.R — Assemble the main-text figures 1, 3, 4, 6, 7, 8 and 9 and
#                     SI Figs S3 and S6 from the plot-ready data saved by the
#                     analysis scripts, in the paper style
#
# Figure 2 (sigma explained) is built by 17_sigma_figures.R, Figure 5 (CO2 as a
# closure tracer) by 15_co2_validity_figure.R. 11_organize_figures.R copies all
# figures into figures/main_text/ and figures/SI/ under their published numbers.
#
# Loads: outputs/02_allan_deviation.RData, 02b_goflux_fluxes.RData, 03_mdf_results.RData,
#        07_treatment.RData, 08_transformation.RData, 09_detection_rates.RData,
#        10_seasonal_confounding.RData, 16_pvalue_vs_mdf.RData, 18_known_truth.RData,
#        20_generality.RData, 21_logging_interval.RData, 22_transform_inference.RData,
#        outputs/precision_comparison.csv
# Writes (figures/, .png/.pdf): main_fig1_processing_chain, main_fig2_precision_detection
#         (Fig. 3), main_fig3_regression_vs_mdf (Fig. 4), main_fig5_r2_seasonal (Fig. 6),
#         main_fig6_treatment_budget (Fig. 7), main_fig9_transformation (Fig. 8),
#         main_fig8_generality (Fig. 9), si_logging_interval (Fig. S3),
#         si_transform_demo (Fig. S6)
# ============================================================================

source("scripts/00_setup.R")
W <- 7.5                          # full page width (in)
th <- theme_paper(base_size = 8)
add_tags <- function(p) p + plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(face = "bold", size = 9))
save_main <- function(p, name, h) {
  ggsave(file.path(fig_dir, paste0(name, ".png")), p, width = W, height = h, dpi = 300)
  ggsave(file.path(fig_dir, paste0(name, ".pdf")), p, width = W, height = h)
  message("Saved ", name)
}
ld <- function(f) { e <- new.env(); load(file.path(output_dir, f), envir = e); e }

a02 <- ld("02_allan_deviation.RData"); a03 <- ld("03_mdf_results.RData")
stem <- a03$df %>% mutate(analyzer = if_else(year == 2025, "Precise analyzer", "Noisier analyzer"))

# ── Figure 3: realized precision and what it does to detection ───────────────
gen <- ld("20_generality.RData")$gen
pc <- read.csv(file.path(output_dir, "precision_comparison.csv"))
camp <- function(src) pc$CH4_precision_ppb[pc$source == src & grepl("^Campaign", pc$method)]
prec <- tibble::tribble(
  ~row,                                        ~gas,   ~spec,          ~sigma,
  "Noisier CH₄ analyzer, unit 1 \u00b7 stems",       "CH4",  PREC_CH4_LGR,   median(stem$sigma_camp_CH4[stem$analyzer == "Noisier analyzer"], na.rm = TRUE),
  "Noisier CH₄ analyzer, unit 1 \u00b7 canopy",      "CH4",  PREC_CH4_LGR,   camp("A-1"),
  "Noisier CH₄ analyzer, unit 2 \u00b7 canopy",      "CH4",  PREC_CH4_LGR,   camp("A-3"),
  "Noisier CH₄ analyzer \u00b7 soils",               "CH4",  PREC_CH4_LGR,   gen$sigma[gen$dataset == "Soils, upland"],
  "Precise CH₄ analyzer \u00b7 stems",               "CH4",  PREC_CH4_7810,  median(stem$sigma_camp_CH4[stem$analyzer == "Precise analyzer"], na.rm = TRUE),
  "Precise CH₄ analyzer \u00b7 farm soil",   "CH4",  0.6,            gen$sigma[gen$dataset == "Agricultural soil, CH4"],
  "N₂O analyzer \u00b7 farm soil",           "N2O",  0.4,            gen$sigma[gen$dataset == "Agricultural soil, N2O"]) %>%
  mutate(ratio = sigma / spec, row = factor(row, levels = rev(row)))
p2a <- ggplot(prec, aes(y = row)) +
  geom_segment(aes(x = spec, xend = sigma, yend = row), colour = AXIS, linewidth = 1.2) +
  geom_point(aes(x = spec), shape = 21, size = 2.2, colour = COL$datasheet, fill = "white", stroke = 0.7) +
  geom_point(aes(x = sigma, colour = ratio >= 1), size = 2.6) +
  scale_colour_manual(values = c(`TRUE` = COL$r2, `FALSE` = COL$mdf), guide = "none") +
  geom_text(aes(x = 12, label = ifelse(ratio >= 1, sprintf("%.1f\u00d7 worse", ratio), sprintf("%.0f\u00d7 better", 1 / ratio))), hjust = 0, size = 2.3, colour = INK2) +
  scale_x_log10(breaks = c(0.1, 0.3, 1, 3, 10), labels = c("0.1", "0.3", "1", "3", "10")) +
  coord_cartesian(xlim = c(0.07, 11), clip = "off") +
  labs(x = "\u03c3 (ppb); open: datasheet, filled: realized", y = NULL) +
  th + theme(panel.grid.major.y = element_blank(), plot.margin = margin(4, 48, 4, 4))

fl_ids <- ld("02b_goflux_fluxes.RData")$df$fl_UniqueID
sig_i <- bind_rows(tibble(analyzer = "Noisier analyzer", s = a02$allan_df_lgr$sc_CH4[a02$allan_df_lgr$UniqueID %in% fl_ids]),
                   tibble(analyzer = "Precise analyzer", s = a02$allan_df_7810$sc_CH4)) %>% filter(is.finite(s), s > 0)
sig_med <- sig_i %>% group_by(analyzer) %>% summarise(m = median(s), .groups = "drop")
p2b <- ggplot(sig_i, aes(s, fill = analyzer)) +
  geom_histogram(bins = 60, colour = NA, alpha = 0.85) +
  geom_vline(data = data.frame(x = c(PREC_CH4_LGR, PREC_CH4_7810), analyzer = c("Noisier analyzer", "Precise analyzer")),
             aes(xintercept = x, colour = analyzer), linetype = "dashed", linewidth = 0.4) +
  scale_x_log10(breaks = c(0.1, 0.3, 1, 3, 10), labels = c("0.1", "0.3", "1", "3", "10"),
                expand = expansion(mult = 0.08)) +
  scale_fill_manual(values = PAL_ANALYZER, name = NULL) + scale_colour_manual(values = PAL_ANALYZER, guide = "none") +
  labs(
       x = expression(sigma[i]~(ppb)), y = "Closures") +
  th + legend_top()

lg <- ld("21_logging_interval.RData")$logsum %>% mutate(rel = sigma / sigma1)
p2c <- ggplot(lg %>% filter(dt > 1 | how == "1 s record"), aes(dt, rel)) +
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
rd_ref <- rd %>% filter(level %in% c("spec", "95"))
rd_rng <- rd %>% filter(level != "spec") %>% group_by(instrument, family) %>%
  summarise(lo = min(pct_below), hi = max(pct_below), .groups = "drop")
fam_col <- c("Datasheet" = COL$datasheet, "Campaign \u03c3 \u00d7 z" = COL$mdf, "Per-closure \u03c3 \u00d7 3t" = RAMP_X3T[2])
fam_lab <- c("Datasheet" = "Datasheet \u03c3", "Campaign \u03c3 \u00d7 z" = "Measured \u03c3, z = 1.96")
p2d <- ggplot(rd_ref, aes(family, pct_below, fill = family)) +
  geom_col(width = 0.7) +
  geom_errorbar(data = rd_rng, aes(x = family, ymin = lo, ymax = hi), inherit.aes = FALSE, width = 0.25,
                linewidth = 0.35, colour = INK) +
  geom_text(data = rd_ref %>% left_join(rd_rng, by = c("instrument", "family")) %>% mutate(y = pmax(pct_below, hi, na.rm = TRUE)),
            aes(y = y, label = sprintf("%.0f%%", pct_below)), vjust = -0.7, size = 2.4, colour = INK) +
  facet_wrap(~ instrument, nrow = 1) +
  scale_fill_manual(values = fam_col, labels = fam_lab, name = "MDF built from") +
  scale_y_continuous(limits = c(0, 70), breaks = c(0, 20, 40, 60), expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "% below MDF") +
  th + theme(axis.text.x = element_blank()) + legend_top()

fig2 <- add_tags(free(p2a) / p2d + plot_layout(heights = c(1, 0.9)))
save_main(fig2, "main_fig2_precision_detection", 5.6)
save_main(p2c + theme(legend.position = "right"), "si_logging_interval", 3.2)

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
p3 <- ggplot(pv, aes(snr_mdf, -log10(p_val), colour = class)) +
  annotate("rect", xmin = 0, xmax = 1, ymin = -log10(0.05), ymax = Inf, fill = TINT_OK) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") + geom_vline(xintercept = 1, linetype = "dashed") +
  geom_point(aes(shape = class), size = 0.8, alpha = 0.7) +
  facet_wrap(~ panel, nrow = 1) +
  scale_x_log10(breaks = c(0.01, 0.1, 1, 10, 100), labels = c("0.01", "0.1", "1", "10", "100")) +
  scale_y_sqrt(breaks = c(0, -log10(0.05), 5, 25, 50, 100),
               labels = c("1", "0.05", "1e-5", "1e-25", "1e-50", "≤1e-100")) +
  scale_colour_manual(values = c("Both" = INK, "Regression only" = COL$reg, "MDF only" = COL$mdf, "Neither" = GREY),
                      name = "Detected by") +
  scale_shape_manual(values = c("Both" = 16, "Regression only" = 16, "MDF only" = 16, "Neither" = 16), name = "Detected by") +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1)), shape = guide_legend(override.aes = list(size = 2))) +
  labs(x = "|flux| / MDF (campaign σ, z = 1.96)", y = "Slope p-value") + th +
  coord_cartesian(xlim = c(0.01, 300))
save_main(p3, "main_fig3_regression_vs_mdf", 3.0)

# ── Figure 6: filters on R² remove periods of small flux ────────────────────────
a10 <- ld("10_seasonal_confounding.RData")
flux_breaks <- c(-10, -1, -0.1, -0.01, 0, 0.01, 0.1, 1, 10, 100)
flux_labs <- c("-10", "-1", "-0.1", "", "0", "", "0.1", "1", "10", "100")
signed_tr <- scales::pseudo_log_trans(sigma = 0.01, base = 10)   # linear within +/-0.01, log beyond
sx <- scale_x_continuous(trans = signed_tr, breaks = flux_breaks, labels = flux_labs)
FIG6_STYLE <- Sys.getenv("FIG6_STYLE", "G")    # E: kept solid dark, removed open terracotta; G: removed faint grey crosses
PAL_R2FATE <- if (FIG6_STYLE == "G") c("Kept" = INK, "Removed by R² ≤ 0.5" = "#c2c2be") else c("Kept" = INK, "Removed by R² ≤ 0.5" = COL$r2)
SHP_R2FATE <- if (FIG6_STYLE == "G") c("Kept" = 16, "Removed by R² ≤ 0.5" = 4) else c("Kept" = 16, "Removed by R² ≤ 0.5" = 1)
r2_panel <- function(d, ttl, note_x) {
  d <- d %>% mutate(fate = factor(if_else(r2 <= 0.5, "Removed by R² ≤ 0.5", "Kept"), levels = names(PAL_R2FATE))) %>%
    arrange(fate)
  ggplot(d, aes(flux, r2, colour = fate)) +
    geom_vline(xintercept = 0, colour = AXIS) +
    geom_hline(yintercept = 0.5, colour = INK2, linewidth = 0.4) +
    geom_point(aes(shape = fate), size = 0.9, stroke = 0.35) +
    sx + scale_y_continuous(breaks = c(0, 0.5, 1)) +
    scale_colour_manual(values = PAL_R2FATE, guide = "none") + scale_shape_manual(values = SHP_R2FATE, guide = "none") +
    labs(title = ttl, x = LAB_FLUX, y = expression(CH[4]~R^2)) + th
}
p5a <- r2_panel(a10$df_gla %>% transmute(flux = CH4_flux_nmolpm2ps, r2 = CH4_r2) %>% filter(is.finite(r2)),
                "Stems", 3)
p5b <- r2_panel(a10$soil %>% transmute(flux, r2) %>% filter(is.finite(r2)),
                "Upland\u2013wetland soils", 3)

fate <- a10$df_gla %>% filter(is.finite(CH4_r2), is.finite(airt)) %>%
  mutate(fate = factor(if_else(CH4_r2 <= 0.5, "Removed by R² ≤ 0.5", "Kept"), levels = names(PAL_R2FATE))) %>%
  arrange(fate)
bin_means <- fate %>% mutate(bin = cut(airt, seq(-10, 35, 5))) %>% group_by(bin) %>%
  summarise(t = mean(airt), `All closures` = mean(CH4_flux_nmolpm2ps),
            `Kept by R² > 0.5` = mean(CH4_flux_nmolpm2ps[CH4_r2 > 0.5]), n = n(), .groups = "drop") %>%
  filter(n >= 15) %>% tidyr::pivot_longer(c(`All closures`, `Kept by R² > 0.5`), names_to = "what", values_to = "m")
p5c <- ggplot(fate, aes(airt, CH4_flux_nmolpm2ps)) +
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
  th + legend_top() + theme(legend.box = "vertical", legend.box.just = "left", legend.spacing.y = unit(2, "pt"))

# (d) share removed by temperature bin: each dataset (stems, soils) x each filter (R2, MDF)
s_soil <- ld("10_seasonal_confounding.RData")$soil
s_soil$UniqueID <- s_soil[[intersect(c("unique_id", "UniqueID"), names(s_soil))[1]]]
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
p5d <- ggplot(tb, aes(bin_mid, pct_flagged, colour = filter, linetype = ds, shape = ds, group = interaction(filter, ds))) +
  geom_line(linewidth = 0.7) + geom_point(size = 1.8, fill = "white", stroke = 0.7) +
  scale_colour_manual(values = c("R² ≤ 0.5" = COL$r2, "Below MDF" = COL$mdf), name = "Removed:") +
  scale_linetype_manual(values = c(Stems = "solid", Soils = "22"), name = NULL) +
  scale_shape_manual(values = c(Stems = 16, Soils = 21), name = NULL) +
  scale_y_continuous(limits = c(0, 100)) +
  labs(x = "Air temperature, 10 °C bins (°C)", y = "Closures removed (%)") +
  th + legend_in(0.98, 0.98) + theme(legend.box = "vertical", legend.box.just = "left", legend.spacing.y = unit(0, "pt"),
             legend.margin = margin(0, 0, 0, 0)) +
  guides(colour = guide_legend(order = 1, override.aes = list(linetype = "solid", shape = NA)),
         linetype = guide_legend(order = 2, override.aes = list(colour = INK)), shape = guide_legend(order = 2))
fig5 <- add_tags((p5a | p5b) / (p5c | p5d) + plot_layout(heights = c(1, 1.15)))
save_main(fig5, "main_fig5_r2_seasonal", 6.4)

# ── Figure 7: below-detection treatment and annual budgets ──────────────────
a07 <- ld("07_treatment.RData")
pcl <- a07$pct_change %>% select(filter, treatment, pct_change_mean, pct_change_median) %>%
  tidyr::pivot_longer(c(pct_change_mean, pct_change_median), names_to = "stat", values_to = "pct") %>%
  mutate(stat = recode(stat, pct_change_mean = "Mean", pct_change_median = "Median"),
         treatment = recode(treatment, "Set below MDF to zero" = "Set to zero", "Remove below MDF" = "Remove"),
         filter = sub("Manufacturer MDF", "Datasheet", sub("Campaign σ ", "Campaign σ\n", sub("Per-closure σ ×3t ", "×3t\n", filter))))
lv <- c("Datasheet", paste0("Campaign σ\n", c("90%", "95%", "99%")), paste0("×3t\n", c("90%", "95%", "99%")))
pcl <- pcl %>% filter(!grepl("3t", filter))
pcl$filter <- factor(pcl$filter, levels = intersect(lv, unique(pcl$filter)))
short_lab <- setNames(c("Spec", "90%", "95%", "99%", "3t 90%", "3t 95%", "3t 99%"), lv)
p6a <- ggplot(pcl, aes(filter, pct, fill = treatment)) +
  geom_hline(yintercept = 0, colour = INK, linewidth = 0.3) +
  geom_col(position = position_dodge(0.75), width = 0.7) +
  facet_wrap(~ stat, nrow = 1) +
  scale_fill_manual(values = c("Set to zero" = COL$zero, "Remove" = COL$mdf), name = "Below-MDF closures") +
  scale_x_discrete(labels = short_lab) +
  labs(
       x = "MDF used to flag", y = "% change") +
  th + theme(legend.position = "top", legend.justification = "left")
cum <- a10$monthly_budget_df %>% filter(treatment != "Campaign σ 95%\n(keep all, flag)") %>%
  mutate(treatment = recode(as.character(treatment), "No filter (keep all)" = "Keep all (flag)",
                            "R² > 0.5 (remove)" = "Remove R² ≤ 0.5",
                            "Campaign σ 95%\n(remove)" = "Remove below MDF"))
fin <- cum %>% filter(month_num == 12)
base <- fin$cum_flux[fin$treatment == "Keep all (flag)"]
fin <- fin %>% mutate(lab = if_else(treatment == "Keep all (flag)", treatment, sprintf("%s (+%d%%)", treatment, as.integer(round(100 * (cum_flux / base - 1) + 1e-9)))))
pal_cum <- c("Keep all (flag)" = COL$keep, "Remove R² ≤ 0.5" = COL$r2, "Remove below MDF" = COL$mdf)
p6b <- ggplot(cum, aes(month_num, cum_flux, colour = treatment)) +
  annotate("rect", xmin = c(0.5, 10.5), xmax = c(4.5, 12.5), ymin = -Inf, ymax = Inf, fill = TINT_NEU) +
  geom_line(linewidth = 0.7) + geom_point(size = 1) +
  scale_colour_manual(values = pal_cum, labels = setNames(fin$lab, fin$treatment), breaks = c("Keep all (flag)", "Remove below MDF", "Remove R² ≤ 0.5"), name = NULL) +
  scale_x_continuous(breaks = 1:12, labels = substr(month.abb, 1, 1)) +
  scale_y_continuous(expand = expansion(mult = c(0.03, 0.12))) +
  guides(colour = guide_legend(ncol = 1)) +
  labs(title = "Stems, 2024",
       x = NULL, y = "Cumulative flux") + th + legend_in(0.02, 0.98)
sb <- a10$soil_budget_df %>% mutate(treatment = recode(as.character(treatment), "No filter (keep all)" = "Keep all (flag)",
                                                       "R² > 0.3 (remove)" = "Remove R² ≤ 0.3",
                                                       "R² > 0.5 (remove)" = "Remove R² ≤ 0.5",
                                                       "R² > 0.7 (remove)" = "Remove R² ≤ 0.7"))
sb <- sb %>% filter(treatment %in% c("Keep all (flag)", "Remove R² ≤ 0.5"))
sfin <- sb %>% group_by(treatment) %>% slice_max(date_mid, n = 1) %>% ungroup()
sbase <- sfin$cum_flux[sfin$treatment == "Keep all (flag)"]
sfin <- sfin %>% mutate(lab = if_else(treatment == "Keep all (flag)", treatment, sprintf("%s (%+.0f%%)", treatment, 100 * (cum_flux / sbase - 1))))
pal_sb <- c("Keep all (flag)" = COL$keep, "Remove R² ≤ 0.3" = R2_LIGHT, "Remove R² ≤ 0.5" = COL$r2, "Remove R² ≤ 0.7" = R2_DARK)
p6c <- ggplot(sb, aes(date_mid, cum_flux, colour = treatment)) +
  geom_line(linewidth = 0.7) + geom_point(size = 1) +
  scale_colour_manual(values = pal_sb, labels = setNames(sfin$lab, sfin$treatment), name = NULL) +
  guides(colour = guide_legend(ncol = 1)) +
  scale_x_date(date_breaks = "3 months", date_labels = "%b %y") +
  labs(title = "Upland\u2013wetland soils", x = NULL, y = "Cumulative flux") + th + legend_in(0.02, 0.02)
sim <- ld("18_known_truth.RData")$score %>% filter(analyzer == "Noisier analyzer") %>%
  tidyr::pivot_longer(c(mean, budget), names_to = "stat", values_to = "err") %>%
  mutate(stat = recode(stat, mean = "Mean", budget = "Budget"),
         treatment = recode(as.character(treatment), "Set below MDF to zero" = "Set to zero", "Remove below MDF" = "Remove below MDF"),
         treatment = factor(treatment, levels = rev(c("Keep all (flag)", "Set to zero", "Remove below MDF", "Remove R² ≤ 0.5"))))
pal_sim <- c("Keep all (flag)" = COL$keep, "Set to zero" = COL$zero, "Remove below MDF" = COL$mdf, "Remove R² ≤ 0.5" = COL$r2)
p6d <- ggplot(sim, aes(err, treatment, colour = treatment, shape = stat)) +
  geom_vline(xintercept = 0, colour = AXIS) +
  stat_summary(fun = median, fun.min = function(v) quantile(v, .05), fun.max = function(v) quantile(v, .95),
               position = position_dodge(width = 0.55), size = 0.3, linewidth = 0.5) +
  scale_colour_manual(values = pal_sim, guide = "none") +
  scale_shape_manual(values = c(Mean = 16, Budget = 1), name = NULL) +
  labs(title = "Simulated stems", x = "Error vs true value (%)", y = NULL) +
  th + legend_top() + theme(panel.grid.major.y = element_blank())
fig6 <- add_tags(((p6a | p6d) + plot_layout(widths = c(2.6, 1))) / (p6b | p6c) + plot_layout(heights = c(1, 1.1)))
save_main(fig6, "main_fig6_treatment_budget", 6.2)

# ── Figure 8: transformation and what it does to inference ──────────────────
a08 <- ld("08_transformation.RData")
s_tr <- c("Forest stems" = median(ld("03_mdf_results.RData")$df$CH4_MDF_ref, na.rm = TRUE),
          "Upland\u2013wetland soils" = median(ld("20_generality.RData")$soil_sigma$mdf_ref, na.rm = TRUE))
lv_tr <- c("Raw", "log10", "arcsinh(flux / MDF)")
tl <- a08$transform_long %>% mutate(transform = sub(" \\(.*", "", transform)) %>% filter(transform != "arcsinh")
tl <- bind_rows(tl, tl %>% filter(transform == "Raw") %>% mutate(transform = "arcsinh(flux / MDF)", value = asinh(value / s_tr[as.character(dataset)]))) %>%
  filter(!is.na(value)) %>% mutate(transform = factor(transform, levels = lv_tr), dataset = factor(dataset, levels = names(s_tr)))
pcnt <- a08$panel_counts %>% mutate(transform = sub(" \\(.*", "", transform)) %>% filter(transform != "arcsinh")
pcnt <- bind_rows(pcnt, pcnt %>% filter(transform == "Raw") %>% mutate(transform = "arcsinh(flux / MDF)")) %>%
  mutate(transform = factor(transform, levels = lv_tr), dataset = factor(dataset, levels = names(s_tr)),
         label = if_else(n_dropped > 0, sprintf("%d of %d kept", n_total - n_dropped, n_total), sprintf("all %d kept", n_total)))
p7a <- ggplot(tl, aes(value)) +
  geom_histogram(bins = 40, fill = "#6f6f6c", colour = NA) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
  ggh4x::facet_grid2(dataset ~ transform, scales = "free", independent = "all") +
  geom_text(data = pcnt, aes(x = Inf, y = Inf, label = label), hjust = 1.05, vjust = 1.3, size = 2.2,
            colour = INK2, inherit.aes = FALSE) +
  labs(x = "Flux, as transformed", y = "Closures") + th + theme(strip.text.y = element_text(angle = 0, hjust = 0))
dm <- tibble(x = seq(-6, 6, by = 0.02)) %>%
  mutate(`arcsinh(flux / MDF)` = asinh(x), `sign log(2|flux / MDF|)` = ifelse(abs(x) > 0.6, sign(x) * log(2 * abs(x)), NA), `flux / MDF` = x) %>%
  tidyr::pivot_longer(-x, names_to = "fn", values_to = "y") %>% filter(!is.na(y)) %>%
  mutate(fn = factor(fn, levels = c("arcsinh(flux / MDF)", "sign log(2|flux / MDF|)", "flux / MDF")))
p7b <- ggplot(dm, aes(x, y, colour = fn, linetype = fn)) +
  annotate("rect", xmin = -1, xmax = 1, ymin = -Inf, ymax = Inf, fill = TINT_NEU) +
  geom_hline(yintercept = 0, colour = AXIS) + geom_vline(xintercept = 0, colour = AXIS) + geom_line(linewidth = 0.6) +
  scale_colour_manual(values = c(INK, MUTED, MUTED), name = NULL) +
  scale_linetype_manual(values = c("solid", "dashed", "dotted"), name = NULL) +
  coord_equal(xlim = c(-6, 6), ylim = c(-6, 6), expand = FALSE) +
  scale_x_continuous(breaks = seq(-6, 6, 3)) + scale_y_continuous(breaks = seq(-6, 6, 3)) +
  labs(x = "Flux / MDF", y = "Transformed value") + th + legend_top() + theme(legend.direction = "vertical")
a22 <- ld("22_transform_inference.RData")
meth_lv <- c("Raw", "log10, positives", "log10(x + 1)", "arcsinh(x)", "arcsinh(x / MDF)", "Ranks")
meth_lab <- c("Raw" = "Raw", "log10, positives" = "log10, positives only", "log10(x + 1)" = "log10(flux + 1)",
              "arcsinh(x)" = "arcsinh(flux)", "arcsinh(x / MDF)" = "arcsinh(flux / MDF)", "Ranks" = "Ranks")
pw <- a22$grp_sim %>% tidyr::pivot_longer(-method, names_to = "scen", values_to = "pct") %>%
  mutate(scen = factor(recode(scen, false_pos = "No true difference", power_2fold = "Fluxes 2\u00d7 larger",
                              power_uptake_share = "More near-zero and uptake"),
                       levels = c("Fluxes 2\u00d7 larger", "More near-zero and uptake", "No true difference")),
         method = factor(method, levels = rev(meth_lv)))
p7c <- ggplot(pw, aes(pct, method, colour = scen, shape = scen)) +
  geom_vline(xintercept = 5, colour = AXIS, linetype = "dashed") +
  geom_point(size = 2.2, stroke = 0.7) +
  scale_colour_manual(values = c(COL$mdf, COL$r2, GREY_D), name = NULL) +
  scale_shape_manual(values = c(16, 17, 1), name = NULL) +
  scale_y_discrete(labels = meth_lab) + scale_x_continuous(limits = c(0, 100)) +
  labs(x = "Simulated comparisons called different (%)", y = NULL) + th + legend_top() + theme(legend.direction = "vertical")
fig7 <- add_tags(wrap_plots(p7a, wrap_plots(p7b, p7c, nrow = 1, widths = c(1, 1.6)), ncol = 1, heights = c(1.1, 1)))
save_main(fig7, "main_fig9_transformation", 6.8)

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
  geom_errorbarh(aes(xmin = ratio_lo, xmax = ratio_hi), height = 0, colour = INK2) +
  geom_point(size = 2.4, colour = INK) +
  facet_wrap(~ scenario, nrow = 1) +
  scale_y_discrete(labels = c(meth_lab, "Raw" = "Raw (ratio of means)")) +
  scale_x_log10(breaks = c(0.5, 1, 2, 4)) +
  labs(x = "How much larger is B reported to be? (true: 2\u00d7)", y = NULL) + th
sfig_tr <- add_tags(wrap_plots(ps_a, ps_b, ncol = 1))
ggsave(file.path(fig_dir, "si_transform_demo.png"), sfig_tr, width = W, height = 5.2, dpi = 300)
ggsave(file.path(fig_dir, "si_transform_demo.pdf"), sfig_tr, width = W, height = 5.2)

# ── Figure 1: the processing chain (conceptual) ─────────────────────────────
steps <- c("Precision", "Detection", "Validity", "Filters", "Treatment", "Transformation")
rows  <- c("Common practice", "Consequence", "Recommended", "Evidence")
txt <- rbind(
  c("Datasheet σ at 1 s", "MDF from datasheet σ; slope p-value alone",
    "Fit statistics or CO₂ used as blanket filters", "R² threshold", "Remove, or set to zero", "log transformation"),
  c("Realized precision differs from the datasheet, in either direction", "Share below detection set by σ source and multiplier; p-values over-confident with correlated noise",
    "Valid small fluxes discarded; leaks pass fit-based checks", "Removes small fluxes of either sign, and with them the periods and places that produce them",
    "Removal inflates means and budgets; zeroing distorts medians", "log drops every value ≤ 0, the low end of the distribution"),
  c("σ from the closures: median of per-closure second-difference precision", "MDF = 1.96 σ / t × flux term, t in seconds; report the regression interval too",
    "Screen the incubation (seal, CO₂ tracer, ambient start), not the flux", "No fit-quality filter by default; report exclusion by season",
    "Keep every value; flag emission / uptake / below MDF", "arcsinh(flux / MDF) or ranks for tests; means and effects on raw values"),
  c("Figs 2\u20133", "Fig 4", "Fig 5", "Fig 6", "Fig 7", "Fig 8"))
cells <- expand.grid(row = seq_along(rows), col = seq_along(steps)) %>%
  mutate(label = mapply(function(r, c) txt[r, c], row, col),
         label = vapply(label, function(z) paste(strwrap(z, 22), collapse = "\n"), ""),
         x = col, y = -row)
fill_row <- c(TINT_NEU, "#f7ebe4", "#e8f1ec", "white")
hdr_col <- rep("#4f6576", 6)
arrow_df <- data.frame(x = seq_along(steps))
p1 <- ggplot() +
  geom_tile(data = cells, aes(x, y, fill = factor(row)), width = 0.96, height = 0.94, colour = NA) +
  geom_text(data = cells, aes(x, y, label = label, fontface = ifelse(row == 4, "bold", "plain")),
            size = 2.15, lineheight = 0.92, colour = INK, family = PAPER_FONT) +
  geom_polygon(data = do.call(rbind, lapply(seq_along(steps), function(i)
    data.frame(i = i, x = i + c(-0.48, 0.40, 0.50, 0.40, -0.48, -0.40), y = -0.3 + c(-0.2, -0.2, 0, 0.2, 0.2, 0)))),
    aes(x, y, group = i, fill = factor(i + 10)), colour = NA) +
  geom_text(data = data.frame(x = seq_along(steps), y = -0.3, l = paste0(seq_along(steps), "  ", steps)),
            aes(x, y, label = l), colour = "white", fontface = "bold", size = 2.15, family = PAPER_FONT) +
  geom_text(data = data.frame(y = -seq_along(rows), l = rows), aes(0.42, y, label = l), hjust = 1,
            size = 2.3, colour = INK2, fontface = "bold", family = PAPER_FONT) +
  annotate("text", x = 3.5, y = -4.55, label = "Across gases, surfaces and analyzers: Fig 9.  Recommendations in brief: Box 1",
           size = 2.3, colour = INK2, fontface = "italic", family = PAPER_FONT) +
  scale_fill_manual(values = c(setNames(fill_row, 1:4), setNames(hdr_col, 11:16)), guide = "none") +
  coord_cartesian(xlim = c(-0.55, 6.5), ylim = c(-4.7, -0.1), clip = "off") +
  theme_void(base_family = PAPER_FONT) + theme(plot.margin = margin(4, 4, 4, 4))
save_main(p1, "main_fig1_processing_chain", 4.0)

# ── Figure 9: the same chain across gases, surfaces and analyzers ───────────
gen <- ld("20_generality.RData")$gen
sw <- ld("18_known_truth.RData")$sweep_sum
short <- c("Stems, noisier analyzer" = "Stems, noisier analyzer", "Stems, precise analyzer" = "Stems, precise analyzer",
           "Soils, upland" = "Upland soils", "Soils, wetland" = "Wetland soils",
           "Agricultural soil, CH4" = "Farm soil, CH₄", "Agricultural soil, N2O" = "Farm soil, N₂O")
pal_rm <- c("Remove below MDF" = COL$mdf, "Remove R² ≤ 0.5" = COL$r2)
pts <- gen %>% transmute(lab = short[dataset], x = pct_below_ref,
                         `Mean|Remove below MDF` = mean_change_mdf, `Mean|Remove R² ≤ 0.5` = mean_change_r2,
                         `Median|Remove below MDF` = median_change_mdf, `Median|Remove R² ≤ 0.5` = median_change_r2) %>%
  tidyr::pivot_longer(-c(lab, x), names_to = "k", values_to = "chg") %>%
  tidyr::separate(k, c("stat", "filt"), sep = "\\|")
labs_df <- pts %>% filter(stat == "Mean") %>% group_by(stat, lab, x) %>% summarise(y = max(chg), .groups = "drop")
swl <- sw %>% select(pct_below, mean_change_mdf, mean_change_r2, median_change_mdf, median_change_r2) %>%
  tidyr::pivot_longer(-pct_below, names_to = "k", values_to = "chg") %>%
  mutate(stat = ifelse(grepl("^mean", k), "Mean", "Median"),
         filt = ifelse(grepl("mdf$", k), "Remove below MDF", "Remove R² ≤ 0.5"), x = pct_below)
p8 <- ggplot() +
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
  labs(x = "Closures below the MDF (%)", y = "Increase in the statistic (%)") +
  th + legend_top()
save_main(p8, "main_fig8_generality", 3.6)

