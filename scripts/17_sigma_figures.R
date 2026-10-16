# ============================================================================
# 17_sigma_figures.R — What the campaign sigma is, how it is derived and applied,
#                      and why the alternative estimators are not used
#
# Writes to figures/; 11_organize_figures.R copies them into main_text/ and SI/:
# Main figure  sigma_explained -> figures/main_text/fig2_sigma_explained.png
#   (a) three closures of one precise analyzer day (near-zero, moderate and
#       large flux): rise from the first reading vs time
#   (b) their first differences: the flux becomes a constant offset (slope x dt)
#       plus noise; a curved closure leaves a trend in them
#   (c) their second differences: centred on zero for any smooth trend; their
#       scaled MAD gives each closure's sigma_i = 1.4826 MAD / sqrt(6) (sigma_closure)
#   (d) sigma_i of every closure that day; the median is the campaign sigma
#
# Supplement
#   worked_example.json (every number of the worked example in Sect. S1)
#   worked_example_mdf.png (Fig. S1; the Supplement includes figures/worked_example_mdf.png
#     directly, 11_organize_figures.R also copies it to figS01)
#   sigma_alternatives -> figures/SI/figS02_sigma_alternatives.png
#
# Loads: outputs/02_allan_deviation.RData, outputs/02b_manID.RData,
#        outputs/02b_goflux_fluxes.RData,
#        analysis/precision_bakeoff/synthetic_results.rds,
#        data/from_hf_stem_study/figure_closures/fig2_sigma_day.csv and worked_example_closures.csv
#        (readings of the closures drawn, from 14_figure_closures.R; no raw-archive access)
# ============================================================================

source("scripts/00_setup.R")

load(file.path(output_dir, "02_allan_deviation.RData"))
load(file.path(output_dir, "02b_manID.RData"))          # manID: stem closure windows
bk <- "analysis/precision_bakeoff"
syn   <- readRDS(file.path(bk, "synthetic_results.rds"))

th <- theme_paper(base_size = 8)

# ---------------------------------------------------------------------------
# Example closure and its analyzer-day
# ---------------------------------------------------------------------------
win <- manID %>% filter(flag == 1) %>% arrange(UniqueID, POSIX.time)
per_closure <- win %>% group_by(UniqueID, analyzer, date_str) %>%
  summarise(n = n(), rise = last(CH4dry_ppb) - first(CH4dry_ppb),
            sigma_i = sigma_closure(CH4dry_ppb, POSIX.time), .groups = "drop")
# The worked example (Sect. S1): two near-linear noisier analyzer closures from one
# well-sampled day, a flux just below the MDF and a large flux (pick_worked_example(), 00_setup.R)
pick <- pick_worked_example(win)
ex_day <- pick$day
ex_lo <- pick$low
ex_hi <- pick$high
ex_id <- ex_hi$UniqueID
# Display the raw-record segments and per-closure sigma that 02 computed (the source of every
# campaign sigma and MDF in the paper), not goFlux's windows: goFlux's window functions truncate
# timestamps to whole seconds, which at ~0.995 s logging adds one duplicate second and one edge
# reading per window. The segments' readings were extracted from the raw record by
# 14_figure_closures.R (data/from_hf_stem_study/figure_closures/worked_example_closures.csv).
fit_env <- new.env(); load(file.path(output_dir, "02b_goflux_fluxes.RData"), envir = fit_env)
stem_fl <- setNames(fit_env$df$fl_UniqueID, fit_env$df$stem_id)
we_raw <- read.csv(file.path(fig_closures_dir, "worked_example_closures.csv"), stringsAsFactors = FALSE,
                   colClasses = c(time_s = "numeric", CH4d_ppm = "numeric"))
seg_02 <- function(stem) {
  s <- we_raw[we_raw$UniqueID == stem, ]
  if (nrow(s) == 0) stop("No readings for worked-example closure ", stem, " in ", fig_closures_dir,
                         "/worked_example_closures.csv; rerun 14_figure_closures.R with the raw records")
  posix <- .POSIXct(s$time_s, tz = "UTC")
  list(t = as.numeric(posix - posix[1], units = "secs"), x = s$CH4d_ppm * 1000, posix = posix)
}
sc_02 <- setNames(allan_df_lgr$sc_CH4, allan_df_lgr$UniqueID)
for (nm in c("ex_lo", "ex_hi")) {
  r <- get(nm); g <- seg_02(r$UniqueID)
  r$t <- list(g$t); r$x <- list(g$x); r$n <- length(g$x)
  r$slope <- unname(coef(lm(g$x ~ g$t))[2]); r$sigma_i <- unname(sc_02[stem_fl[[r$UniqueID]]])
  assign(nm, r)
}
day_rows <- per_closure %>% filter(analyzer == "LGR1", date_str == ex_day) %>%
  mutate(sigma_i = unname(sc_02[stem_fl[UniqueID]]))           # 02's per-closure sigma (MDF source)
sig_day <- mad_days$mad_sd_CH4[mad_days$instrument == "LGR1" & mad_days$date_str == ex_day][1]
wr_day <- mad_days$wr_sd_CH4[mad_days$instrument == "LGR1" & mad_days$date_str == ex_day][1]

# Walkthrough on the precise analyzer, whose noise (~0.1 ppb) is small enough that the flux
# is visible in the first differences: three closures of one day with near-zero, moderate and
# large fluxes; their records, first differences (offset by the flux, and by any curvature)
# and second differences (all centred on zero); then sigma_i of every closure that day and
# their median, the campaign sigma.
wk_day <- "2025-06-03"
# Every closure window of that day (02b), from
# data/from_hf_stem_study/figure_closures/fig2_sigma_day.csv (14_figure_closures.R)
wk <- read.csv(file.path(fig_closures_dir, "fig2_sigma_day.csv"), stringsAsFactors = FALSE,
               colClasses = c(time_s = "numeric", CH4dry_ppb = "numeric")) %>%
  mutate(POSIX.time = .POSIXct(time_s, tz = "UTC")) %>%
  arrange(UniqueID, POSIX.time) %>% group_by(UniqueID) %>%
  mutate(t = as.numeric(POSIX.time - POSIX.time[1], units = "secs")) %>% ungroup()
wk_s <- wk %>% group_by(UniqueID) %>%
  summarise(slope = unname(coef(lm(CH4dry_ppb ~ t))[2]),
            sigma_i = sigma_closure(CH4dry_ppb, POSIX.time), .groups = "drop")
wk_pick <- c("s1307" = "Near zero", "s1289" = "Moderate", "s1295" = "Large")
pal3 <- setNames(RAMP_BLUE, wk_pick)
ex3 <- wk %>% filter(UniqueID %in% names(wk_pick)) %>% mutate(cl = factor(wk_pick[UniqueID], levels = wk_pick)) %>%
  group_by(cl) %>% mutate(x0 = CH4dry_ppb - CH4dry_ppb[1]) %>% ungroup()
w1 <- ex3 %>% group_by(cl) %>% reframe(t = t[-1], d = diff(CH4dry_ppb))
w2 <- ex3 %>% group_by(cl) %>% reframe(t = t[-(1:2)], d = diff(CH4dry_ppb, differences = 2))
wk_sig <- median(wk_s$sigma_i)
wk_whole <- mad(unlist(lapply(split(wk$CH4dry_ppb, wk$UniqueID), diff))) / sqrt(2)   # first differences, all closures pooled
message("Walkthrough day ", wk_day, ": campaign sigma ", round(wk_sig, 3), " ppb; pooled first-difference ", round(wk_whole, 3), " ppb; slopes ",
        paste(round(wk_s$slope[match(names(wk_pick), wk_s$UniqueID)], 3), collapse = ", "),
        "; sigma_i ", paste(round(wk_s$sigma_i[match(names(wk_pick), wk_s$UniqueID)], 3), collapse = ", "))
ylim_d <- range(c(w1$d, w2$d))
pa <- ggplot(ex3, aes(t, x0, colour = cl)) + geom_line(linewidth = 0.5) +
  scale_colour_manual(values = pal3, name = "Flux") +
  labs(x = NULL, y = "Rise from start (ppb)") + th +
  theme(legend.position = c(0.02, 0.98), legend.justification = c(0, 1), legend.key.size = unit(8, "pt"),
        legend.background = element_blank(), axis.text.x = element_blank())
pb <- ggplot(w1, aes(t, d, colour = cl)) + geom_hline(yintercept = 0, colour = AXIS) + geom_line(linewidth = 0.35) +
  scale_colour_manual(values = pal3, guide = "none") + scale_y_continuous(limits = ylim_d) +
  labs(x = NULL, y = "First difference (ppb)") + th + theme(axis.text.x = element_blank())
pc <- ggplot(w2, aes(t, d, colour = cl)) + geom_hline(yintercept = 0, colour = AXIS) + geom_line(linewidth = 0.35) +
  scale_colour_manual(values = pal3, guide = "none") + scale_y_continuous(limits = ylim_d) +
  labs(x = "Time in closure (s)", y = "Second difference (ppb)") + th
wk_s <- wk_s %>% arrange(sigma_i) %>% mutate(rank = row_number(), cl = unname(wk_pick[UniqueID]))
pd <- ggplot(wk_s, aes(sigma_i, rank)) +
  geom_segment(aes(x = 0, xend = sigma_i, yend = rank), colour = "#e6e3dc", linewidth = 0.5) +
  geom_vline(xintercept = wk_sig, colour = INK, linewidth = 0.6) +
  geom_vline(xintercept = wk_whole, colour = MUTED, linewidth = 0.5, linetype = "22") +
  annotate("text", x = wk_whole, y = 1, label = "pooled first\ndifferences",
           hjust = 1.05, vjust = 0, size = 2.3, colour = INK2, lineheight = 0.9) +
  geom_point(colour = GREY, size = 1.5) +
  geom_point(data = wk_s %>% filter(!is.na(cl)), aes(colour = cl), size = 2.4) +
  annotate("text", x = wk_sig, y = max(wk_s$rank) + 1.3, label = "median", hjust = -0.1, size = 2.4, colour = INK2) +
  scale_colour_manual(values = pal3, guide = "none") +
  scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.08))) +
  scale_y_continuous(breaks = NULL, expand = expansion(add = c(0.6, 2.6))) +
  labs(x = expression(sigma[i]~"of each closure that day (ppb)"), y = NULL) + th

fig <- ((pa / pb / pc) | pd) + plot_layout(widths = c(1.7, 1)) +
  plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(face = "bold", size = 9))
ggsave(file.path(fig_dir, "sigma_explained.png"), fig, width = 7.5, height = 5.6, dpi = 300)
ggsave(file.path(fig_dir, "sigma_explained.pdf"), fig, width = 7.5, height = 5.6)
message("Saved fig_sigma_explained (example closure ", ex_id, ", ", ex_day, ")")

# ---------------------------------------------------------------------------
# Worked example (Supplement Sect. S1): every number from the two noisier
# analyzer closures picked above (ex_lo, ex_hi), written to figures/worked_example.json
# ---------------------------------------------------------------------------
cb <- as.data.frame(fit_env$ch4_best)
worked <- lapply(list(low = ex_lo, high = ex_hi), function(r) {
  id <- r$UniqueID; tt <- r$t[[1]]; xx <- r$x[[1]]
  f <- cb[cb$UniqueID == id, ][1, ]
  g <- manID[manID$UniqueID == id & manID$flag == 1, ][1, ]
  dd1 <- diff(xx); dd2 <- diff(xx, differences = 2)
  sxx <- sum((tt - mean(tt))^2); nn <- length(xx)
  se_plan <- f$det.prec / sqrt(sxx) * f$flux.term                 # Cowan et al. (2025), campaign sigma
  list(UniqueID = id, date = ex_day, n = nn, dt = logging_interval(as.POSIXct(tt, origin = "1970-01-01")),
       span = diff(range(tt)), t = f$det.t, rise = tail(xx, 1) - head(xx, 1), slope = r$slope,
       d1_median = median(dd1), d1_sd = sd(dd1), d1_mad_sigma = mad(dd1) / sqrt(2),
       d2_mad_raw = mad(dd2, constant = 1), sigma_i = r$sigma_i,
       sigma_campaign = f$det.prec, z = Z95, dC = Z95 * f$det.prec, rate = Z95 * f$det.prec / f$det.t,
       Vtot_L = g$Vtot, Area_cm2 = g$Area, Pcham_kPa = g$Pcham, Tcham_C = g$Tcham,
       flux_term = f$flux.term, MDF = f$det.MDF, flux = f$best.flux, model = f$model,
       ratio = f$best.flux / f$det.MDF, class = f$det.class,
       cowan_half = qt(0.975, nn - 2) * se_plan, lm_half = qt(0.975, nn - 2) * f$LM.SE, lm_p = f$LM.p.val)
})
worked$day <- list(n_closures = nrow(day_rows), sigma_median = sig_day, sigma_min = min(day_rows$sigma_i, na.rm = TRUE),
                   sigma_max = max(day_rows$sigma_i, na.rm = TRUE), whole_record = wr_day)
jsonlite::write_json(worked, file.path(fig_dir, "worked_example.json"), auto_unbox = TRUE, digits = 8, pretty = TRUE)

wedge <- do.call(rbind, lapply(c("low", "high"), function(k) {
  w <- worked[[k]]; r <- if (k == "low") ex_lo else ex_hi
  data.frame(closure = ifelse(k == "low", sprintf("Flux near the MDF (%.2f \u00d7 MDF)", w$ratio), sprintf("Large flux (%.0f \u00d7 MDF)", w$ratio)),
             t = r$t[[1]], y = r$x[[1]] - r$x[[1]][1], mdf_rate = w$MDF / w$flux_term, slope = r$slope) }))
wedge$closure <- factor(wedge$closure, levels = unique(wedge$closure))
we <- ggplot(wedge, aes(t, y)) +
  geom_ribbon(data = wedge, aes(x = t, ymin = -mdf_rate * t, ymax = mdf_rate * t), inherit.aes = FALSE,
              fill = TINT_NEU, alpha = 1) +
  geom_point(size = 0.5, colour = INK2) +
  geom_abline(aes(intercept = 0, slope = slope), colour = INK, linewidth = 0.6) +
  facet_wrap(~ closure, scales = "free_y") +
  labs(x = "Time in window (s)", y = expression(CH[4]~"change from first reading (ppb)"),
       title = "The MDF as a slope: grey wedge = ±1.96 sigma / t (concentration change a flux at the MDF would produce)") + th
ggsave(file.path(fig_dir, "worked_example_mdf.png"), we, width = 7.5, height = 3, dpi = 300)
message("Saved worked_example.json and worked_example_mdf.png")


# ---------------------------------------------------------------------------
# SI: all estimators on simulated records
# ---------------------------------------------------------------------------
keep_sc <- c("1 baseline, small fluxes, 1 s", "2 large heterogeneous fluxes, 1 s", "3 large fluxes, 10 s logging",
             "5 large + curvature, 10 s", "6 ebullition (30 % bubbly)", "8 low-noise analyzer (sigma 0.1)", "11 linear drift 0.01 ppb/s")
sc_lab <- c("Small fluxes", "Large, variable fluxes", "Large fluxes, 10 s logging", "Curved closures, 10 s logging",
            "Bubbles in 30% of closures", "Very precise analyzer", "Linear instrument drift")
sf <- syn %>% filter(scenario %in% keep_sc) %>%
  select(scenario, E1_rec_d1, E3_win_centred, E5_median_closure_d2) %>%
  tidyr::pivot_longer(-scenario, names_to = "est", values_to = "ratio") %>%
  mutate(est = factor(recode(est, E1_rec_d1 = "First differences, whole record",
                             E3_win_centred = "First differences, per closure",
                             E5_median_closure_d2 = "Second differences, closure median (used here)"),
                      levels = c("Second differences, closure median (used here)", "First differences, per closure",
                                 "First differences, whole record")),
         scenario = factor(sc_lab[match(scenario, keep_sc)], levels = rev(sc_lab)))
sfig <- ggplot(sf, aes(ratio, scenario, colour = est)) +
  geom_vline(xintercept = 1, colour = AXIS) +
  stat_summary(fun = median, fun.min = function(v) quantile(v, .1), fun.max = function(v) quantile(v, .9),
               position = position_dodge(width = 0.6), size = 0.3, linewidth = 0.5) +
  scale_colour_manual(values = c(INK, "#b3b3ae", "#6f6f6c"), name = NULL, guide = guide_legend(ncol = 1)) +
  scale_x_log10(breaks = c(0.8, 1, 1.25, 1.5, 2, 2.5)) +
  labs(x = expression("Estimated / true"~sigma~"(log scale)"), y = NULL) +
  th + legend_top()
ggsave(file.path(fig_dir, "sigma_alternatives.png"), sfig, width = 6, height = 4, dpi = 300)
ggsave(file.path(fig_dir, "sigma_alternatives.pdf"), sfig, width = 6, height = 4)
