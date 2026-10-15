# ============================================================================
# 18_known_truth_simulation.R — Filters and below-detection treatments scored
#                               against a KNOWN true flux
#
# Every bias in the paper's main results is measured against "keep every
# value", because the true flux of a real closure is unknown. Here the truth is
# known. Simulated field years of stem closures:
#   true flux: a seasonal (temperature-driven, Q10 = 2.5) lognormal emission
#     (median 0.12 nmol m-2 s-1 at 15 C, sdlog 1.3) with 10 % small true uptake,
#     calibrated so that the fractions below the MDF match the stem data
#     (Noisier analyzer: ~53 % overall, ~75 % in winter, ~30 % in summer; Precise analyzer: ~2 %)
#   each closure: 180 readings at 1 Hz, slope = flux / flux term, plus white
#     noise at the realized precision of Noisier analyzer (3.6 ppb) or B (0.097 ppb)
#   processing as in the paper: OLS slope, R2, p-value; sigma_i from second
#     differences, campaign sigma = monthly median; MDF = 1.96 sigma / t x flux term
# Scored: mean, median and annual budget (sum of monthly means) under keep /
# set-to-zero / remove-below-MDF / remove R2 <= 0.5, as % error from the truth;
# and detection calibration on zero-flux closures (white and red noise).
#
# Writes: outputs/18_known_truth.RData, figures/known_truth_simulation.png/.pdf
# ============================================================================

source("scripts/00_setup.R")
set.seed(20261013)

FLUX_TERM <- 2.2          # nmol m-2 s-1 per ppb s-1 (stem median, both analyzers)
N_OBS     <- 180          # readings per closure (1 Hz)
T_SEC     <- N_OBS        # window span + one logging interval
SIGMAS    <- c("Noisier analyzer" = 3.6, "Precise analyzer" = 0.097)
N_YEARS   <- 200          # replicate field years
N_MONTH   <- 80           # closures per month
temp_m    <- 11 + 14 * sin(2 * pi * ((1:12) - 4.5) / 12)   # ~ -3 C (Jan) to 25 C (Jul)

tt  <- seq_len(N_OBS) - 1
tc  <- tt - mean(tt); Sxx <- sum(tc^2)

true_flux <- function(n, temp) {
  med <- 0.12 * 2.5^((temp - 15) / 10)
  f <- rlnorm(n, log(med), 1.3)
  up <- runif(n) < 0.10
  f[up] <- -abs(rnorm(sum(up), 0.02, 0.01))
  f
}

# process a block of closures: rows = closures, cols = readings
process <- function(f, sigma, noise = c("white", "red"), phi = 0.6) {
  noise <- match.arg(noise); n <- length(f)
  eps <- matrix(rnorm(n * N_OBS, 0, sigma), n)
  if (noise == "red") {                       # AR(1) with the same marginal sd
    for (j in 2:N_OBS) eps[, j] <- phi * eps[, j - 1] + sqrt(1 - phi^2) * eps[, j]
  }
  x <- outer(f / FLUX_TERM, tt) + eps
  xm <- rowMeans(x)
  slope <- as.vector(x %*% tc) / Sxx
  fit_ss <- slope^2 * Sxx
  tot_ss <- rowSums((x - xm)^2)
  res_ss <- pmax(tot_ss - fit_ss, 1e-12)
  se <- sqrt(res_ss / (N_OBS - 2) / Sxx)
  d2 <- x[, 3:N_OBS] - 2 * x[, 2:(N_OBS - 1)] + x[, 1:(N_OBS - 2)]
  sig_i <- apply(d2, 1, function(v) mad(v)) / sqrt(6)
  data.frame(true = f, flux = slope * FLUX_TERM, r2 = fit_ss / tot_ss,
             p = 2 * pt(-abs(slope / se), N_OBS - 2), sig_i = sig_i)
}

# ── 1. Field years: treatments scored against the truth ─────────────────────
score <- list()
for (an in names(SIGMAS)) {
  message("Simulating ", an)
  for (y in seq_len(N_YEARS)) {
    d <- do.call(rbind, lapply(1:12, function(m) {
      f <- true_flux(N_MONTH, temp_m[m])
      cbind(month = m, process(f, SIGMAS[[an]]))
    }))
    sig_c <- ave(d$sig_i, d$month, FUN = median)               # campaign sigma per month
    mdf <- Z95 * sig_c / T_SEC * FLUX_TERM
    below <- abs(d$flux) < mdf
    lowr2 <- d$r2 <= 0.5
    est <- list(`Keep all (flag)` = d$flux,
                `Set below MDF to zero` = ifelse(below, 0, d$flux),
                `Remove below MDF` = ifelse(below, NA, d$flux),
                `Remove R² ≤ 0.5` = ifelse(lowr2, NA, d$flux))
    truth_mean <- mean(d$true); truth_med <- median(d$true)
    truth_bud <- sum(tapply(d$true, d$month, mean))
    for (k in names(est)) {
      v <- est[[k]]
      score[[length(score) + 1]] <- data.frame(
        analyzer = an, year = y, treatment = k,
        mean = 100 * (mean(v, na.rm = TRUE) / truth_mean - 1),
        median = 100 * (median(v, na.rm = TRUE) / truth_med - 1),
        budget = 100 * (sum(local({b <- tapply(v, d$month, mean, na.rm = TRUE); b[!is.finite(b)] <- 0; b})) / truth_bud - 1),   # empty month = 0, as in script 10
        removed_cold = 100 * mean(is.na(v[d$month %in% c(1:3, 11:12)])),
        removed_warm = 100 * mean(is.na(v[d$month %in% 6:8])),
        mdf = median(mdf))
    }
  }
}
score <- do.call(rbind, score)
score$treatment <- factor(score$treatment, levels = c("Keep all (flag)", "Set below MDF to zero", "Remove below MDF", "Remove R² ≤ 0.5"))
truth_summary <- score %>% group_by(analyzer, treatment) %>%
  summarise(across(c(mean, median, budget, removed_cold, removed_warm),
                   list(med = ~ median(.x, na.rm = TRUE), lo = ~ quantile(.x, 0.05, na.rm = TRUE), hi = ~ quantile(.x, 0.95, na.rm = TRUE))),
            mdf = median(mdf), .groups = "drop")
print(as.data.frame(truth_summary %>% select(analyzer, treatment, mean_med, median_med, budget_med, removed_cold_med, removed_warm_med)), digits = 3)

# ── 2. Detection calibration: zero-flux closures and a power curve ──────────
calib <- list()
for (noise in c("white", "red")) {
  d0 <- process(rep(0, 20000), 1, noise)
  mdf0 <- Z95 * median(d0$sig_i) / T_SEC * FLUX_TERM
  calib[[noise]] <- data.frame(noise = noise, criterion = c("|flux| > MDF", "p < 0.05"),
                               false_pos = 100 * c(mean(abs(d0$flux) > mdf0), mean(d0$p < 0.05)))
}
calib <- do.call(rbind, calib)
print(calib)
pow <- do.call(rbind, lapply(c("white", "red"), function(noise) {
  ref <- process(rep(0, 4000), 1, noise); mdf0 <- Z95 * median(ref$sig_i) / T_SEC * FLUX_TERM
  do.call(rbind, lapply(c(0, 0.05, 0.1, 0.2, 0.3, 0.5, 0.75, 1, 1.5, 2), function(r) {
    d1 <- process(rep(r * mdf0, 4000), 1, noise)
    data.frame(noise = noise, ratio = r,
               crit = c("|flux| > MDF", "p < 0.05"),
               detected = 100 * c(mean(abs(d1$flux) > mdf0), mean(d1$p < 0.05)))
  }))
}))

# ── 3. Two groups with different flux levels: does processing distort a comparison? ──
grp <- do.call(rbind, lapply(seq_len(100), function(y) {
  do.call(rbind, lapply(c(low = 1, high = 3), function(mult) {
    d <- do.call(rbind, lapply(1:12, function(m) cbind(month = m, process(true_flux(N_MONTH / 2, temp_m[m]) * mult, SIGMAS[["Noisier analyzer"]]))))
    mdf <- Z95 * ave(d$sig_i, d$month, FUN = median) / T_SEC * FLUX_TERM
    data.frame(year = y, group = ifelse(mult == 1, "low", "high"), truth = mean(d$true),
               keep = mean(d$flux), remove_mdf = mean(d$flux[abs(d$flux) >= mdf]), remove_r2 = mean(d$flux[d$r2 > 0.5]))
  }))
}))
ratio <- grp %>% tidyr::pivot_longer(c(truth, keep, remove_mdf, remove_r2), names_to = "what", values_to = "m") %>%
  tidyr::pivot_wider(names_from = group, values_from = m) %>% mutate(ratio = high / low) %>%
  group_by(what) %>% summarise(ratio_med = median(ratio), lo = quantile(ratio, .05), hi = quantile(ratio, .95), .groups = "drop")
print(ratio)

# ── 4. Noise sweep: removal bias against the share of closures below the MDF ──
sweep <- do.call(rbind, lapply(c(0.05, 0.1, 0.25, 0.5, 1, 2, 3.6, 6, 10), function(sg) {
  do.call(rbind, lapply(seq_len(20), function(y) {
    d <- do.call(rbind, lapply(1:12, function(m) cbind(month = m, process(true_flux(N_MONTH, temp_m[m]), sg))))
    mdf <- Z95 * ave(d$sig_i, d$month, FUN = median) / T_SEC * FLUX_TERM
    below <- abs(d$flux) < mdf
    data.frame(sigma = sg, pct_below = 100 * mean(below),
               mean_change_mdf = 100 * (mean(d$flux[!below]) / mean(d$flux) - 1),
               mean_change_r2 = 100 * (mean(d$flux[d$r2 > 0.5]) / mean(d$flux) - 1),
               median_change_mdf = 100 * (abs(median(d$flux[!below])) / abs(median(d$flux)) - 1),
               median_change_r2 = 100 * (abs(median(d$flux[d$r2 > 0.5])) / abs(median(d$flux)) - 1))
  }))
}))
sweep_sum <- sweep %>% group_by(sigma) %>% summarise(across(everything(), median), .groups = "drop")
print(sweep_sum)

save(score, truth_summary, calib, pow, grp, ratio, sweep, sweep_sum, file = file.path(output_dir, "18_known_truth.RData"))

# ── SI figure ───────────────────────────────────────────────────────────────
sl <- score %>% tidyr::pivot_longer(c(mean, median, budget), names_to = "stat", values_to = "err") %>%
  mutate(stat = factor(recode(stat, mean = "Mean", median = "Median", budget = "Annual budget"),
                       levels = c("Mean", "Median", "Annual budget")))
pal_t <- c("Keep all (flag)" = COL$keep, "Set below MDF to zero" = COL$zero,
           "Remove below MDF" = COL$mdf, "Remove R² ≤ 0.5" = COL$r2)
s1 <- ggplot(sl, aes(treatment, err, colour = treatment)) +
  geom_hline(yintercept = 0, colour = AXIS) +
  stat_summary(fun = median, fun.min = function(v) quantile(v, .05), fun.max = function(v) quantile(v, .95),
               size = 0.35, linewidth = 0.6) +
  facet_grid(analyzer ~ stat) +
  scale_colour_manual(values = pal_t, name = NULL) +
  scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 5, base = 10),
                     breaks = c(-100, 0, 100, 300)) +
  labs(x = NULL, y = "Error relative to the true value (%)") +
  theme_paper(base_size = 8) + theme(axis.text.x = element_blank()) + legend_top() +
  guides(colour = guide_legend(nrow = 2))
s2 <- ggplot(pow, aes(ratio, detected, colour = crit, linetype = noise)) +
  geom_hline(yintercept = c(5), colour = AXIS) +
  geom_line(linewidth = 0.6) + geom_point(size = 1) +
  scale_colour_manual(values = c("|flux| > MDF" = COL$mdf, "p < 0.05" = COL$reg), name = NULL) +
  scale_linetype_manual(values = c(white = "solid", red = "22"), name = NULL, breaks = c("white", "red"),
                        labels = c(white = "white noise", red = "correlated noise")) +
  labs(x = "True flux / MDF", y = "% of closures detected") +
  theme_paper(base_size = 8) + legend_top() + theme(legend.box = "vertical", legend.box.just = "left", legend.spacing.y = unit(2, "pt")) +
  guides(colour = guide_legend(nrow = 2, order = 1), linetype = guide_legend(nrow = 2, order = 2))
sfig <- (s1 | s2) + plot_layout(widths = c(2.2, 1)) + plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 9))
ggsave(file.path(fig_dir, "known_truth_simulation.png"), sfig, width = 7.5, height = 4.2, dpi = 300)
ggsave(file.path(fig_dir, "known_truth_simulation.pdf"), sfig, width = 7.5, height = 4.2)
message("Saved known_truth_simulation")
