# Synthetic field days with known noise and drift; scores every estimator.
# Run from the repo root: Rscript analysis/precision_bakeoff/synthetic.R
source("analysis/precision_bakeoff/estimators.R")
set.seed(20261002)

# One field day: n closures of length Tw (s), separated by ambient gaps, logged every dt s.
sim_day <- function(n = 30, Tw = 180, dt = 1, sigma = 1, amb = 2000,
                    slope = function(n) rnorm(n, 0, 0.02),   # ppb/s; true flux slopes
                    curv = 0,                                  # HM saturation: k (1/s); 0 = linear
                    bubbles = 0, bubble_size = 200,            # expected steps per closure; ppb
                    conc_noise = 0,                            # sigma scales as (C/amb)^conc_noise
                    drift = 0, rw = 0,                         # linear drift ppb/s; random-walk sd per sqrt(s)
                    ar = 0,                                    # AR(1) coefficient of the noise
                    alternate = FALSE, place_jump = 30) {
  s_true <- slope(n); rows <- list(); t0 <- 0
  for (i in seq_len(n)) {
    gap <- runif(1, 90, 240)
    ta <- seq(t0, t0 + gap - dt, by = dt)
    rows[[length(rows) + 1]] <- data.frame(t = ta, sig = 0, id = NA, win = FALSE, gap = i, sid = i, part = "amb")
    tw <- seq(t0 + gap, t0 + gap + Tw - dt, by = dt); tau <- tw - tw[1]
    sig <- if (curv > 0) s_true[i] / curv * (1 - exp(-curv * tau)) else s_true[i] * tau
    if (bubbles > 0) { nb <- rpois(1, bubbles); if (nb) for (b in runif(nb, 0, Tw)) sig <- sig + (tau >= b) * rexp(1, 1 / bubble_size) }
    sig <- sig + place_jump * (tau < 15) * (1 - tau / 15)             # chamber-placement transient
    rows[[length(rows) + 1]] <- data.frame(t = tw, sig = sig, id = i, win = tau >= 15, gap = NA, sid = i, part = "win")
    t0 <- max(tw) + dt
  }
  r <- do.call(rbind, rows); N <- nrow(r)
  base <- amb + drift * r$t + if (rw > 0) cumsum(rnorm(N, 0, rw * sqrt(dt))) else 0
  C <- base + r$sig
  e <- rnorm(N); if (ar > 0) e <- as.numeric(stats::filter(e * sqrt(1 - ar^2), ar, method = "recursive"))
  r$x <- C + sigma * (pmax(C, 1) / amb)^conc_noise * e
  if (alternate) r$x[seq(2, N, by = 2)] <- r$x[seq(1, N - 1, by = 2)]      # every other row carried forward
  # shoulders: last 60 s of ambient before a closure, first 60 s after it
  r$shoulder <- NA_character_
  for (i in seq_len(n)) {
    w <- which(r$id == i)
    pre <- which(r$part == "amb" & r$sid == i & r$t >= r$t[min(w)] - 60)
    post <- which(r$part == "amb" & r$sid == i + 1 & r$t <= r$t[max(w)] + 60)
    r$shoulder[pre] <- "pre"; r$shoulder[post] <- "post"; r$sid[post] <- i
  }
  r$id[!r$win] <- NA
  list(r = r, s_true = s_true, sigma = sigma, drift = drift)
}

big <- function(n) sample(c(1, -1), n, TRUE, c(.85, .15)) * rlnorm(n, log(0.5), 1.5)   # wide, up to ~50 ppb/s
scen <- list(
  `1 baseline, small fluxes, 1 s`      = list(),
  `2 large heterogeneous fluxes, 1 s`  = list(slope = big),
  `3 large fluxes, 10 s logging`       = list(slope = big, dt = 10),
  `4 large + strong curvature, 1 s`    = list(slope = big, curv = 1/120),
  `5 large + curvature, 10 s`          = list(slope = big, curv = 1/120, dt = 10),
  `6 ebullition (30 % bubbly)`         = list(slope = big, bubbles = 0.6),
  `7 noise grows with conc.`           = list(slope = big, conc_noise = 1),
  `8 low-noise analyzer (sigma 0.1)`   = list(slope = function(n) abs(rnorm(n, 0, 0.3)), sigma = 0.1),
  `9 red noise (AR 0.5)`               = list(ar = 0.5),
  `10 alternating-gas logger, 1 s`     = list(alternate = TRUE),
  `11 linear drift 0.01 ppb/s`         = list(drift = 0.01),
  `12 random-walk drift`               = list(rw = 0.05)
)
est <- list(E1_rec_d1 = sig_rec_d1, E2_rec_d2 = sig_rec_d2, E3_win_centred = sig_win_cent,
            E4_win_d2 = sig_win_d2, E5_median_closure_d2 = sig_med_d2, E6_win_quadres = sig_win_quad)
reps <- 12
res <- do.call(rbind, lapply(names(scen), function(s) {
  do.call(rbind, lapply(seq_len(reps), function(k) {
    d <- do.call(sim_day, scen[[s]]); r <- d$r
    data.frame(scenario = s, rep = k, as.list(sapply(est, function(f) f(r) / d$sigma)),
               white = diag_white(r), leak = diag_leak(r), conc = diag_conc(r), zero = diag_zero(r),
               drift_true = d$drift, drift_shoulder = drift_shoulder(r), drift_baseline = drift_baseline(r))
  }))
}))
saveRDS(res, "analysis/precision_bakeoff/synthetic_results.rds")
agg <- aggregate(. ~ scenario, res[, setdiff(names(res), "rep")], median, na.action = na.pass)
agg <- agg[order(as.integer(sub(" .*", "", agg$scenario))), ]
cat("\nsigma-hat / sigma_true (median of", reps, "simulated days); diagnostics\n")
print(format(cbind(agg[, 1, drop = FALSE], round(agg[, 2:11], 2)), width = 6), row.names = FALSE)
cat("\nDrift recovery (ppb/s): true vs shoulder vs baseline\n")
print(cbind(agg["scenario"], round(agg[, c("drift_true", "drift_shoulder", "drift_baseline")], 4)), row.names = FALSE)
write.csv(agg, "analysis/precision_bakeoff/synthetic_summary.csv", row.names = FALSE)
