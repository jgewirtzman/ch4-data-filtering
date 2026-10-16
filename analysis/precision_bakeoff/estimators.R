# Noise (sigma) and drift estimators + diagnostics for chamber concentration records.
# A record is a data.frame with t (s), x (concentration), id (closure id or NA outside
# windows), win (TRUE inside the fitted window), shoulder ("pre"/"post"/NA).

mad_s <- function(d) stats::mad(d[is.finite(d)], constant = 1.4826)

# ---- sigma estimators -------------------------------------------------------
sig_rec_d1 <- function(r) mad_s(diff(r$x)) / sqrt(2)                       # E1 whole record, 1st diffs (earlier campaign sigma)
sig_rec_d2 <- function(r) mad_s(diff(r$x, differences = 2)) / sqrt(6)      # E2 whole record, 2nd diffs
by_win <- function(r, f) lapply(split(r$x[r$win], r$id[r$win]), f)
sig_win_cent <- function(r) mad_s(unlist(by_win(r, function(v) { d <- diff(v); d - stats::median(d) }))) / sqrt(2)  # E3 per closure, centred 1st diffs
sig_win_d2   <- function(r) mad_s(unlist(by_win(r, function(v) diff(v, differences = 2)))) / sqrt(6)               # E4 pooled 2nd diffs
sig_i_d2     <- function(r) unlist(by_win(r, function(v) mad_s(diff(v, differences = 2)) / sqrt(6)))               # per closure
sig_med_d2   <- function(r) stats::median(sig_i_d2(r), na.rm = TRUE)                                               # E5 median of closures
sig_win_quad <- function(r) mad_s(unlist(lapply(split(r[r$win, ], r$id[r$win]), function(s) {                       # E6 quadratic residuals
  if (nrow(s) < 8) return(NULL); stats::residuals(stats::lm(x ~ poly(t, 2), data = s)) })))

# ---- drift estimators -------------------------------------------------------
# D1: slope of the ambient shoulders around each closure (pre and post fitted separately),
#     median across closures -> instrument drift rate (concentration units / s)
shoulder_slopes <- function(r) {
  s <- r[!is.na(r$shoulder), ]
  out <- lapply(split(s, paste(s$sid, s$shoulder)), function(z) if (nrow(z) >= 5) unname(stats::coef(stats::lm(x ~ t, data = z))[2]) else NA_real_)
  unlist(out) }
drift_shoulder <- function(r) stats::median(shoulder_slopes(r), na.rm = TRUE)
# D2: day-scale baseline drift: regress the median of each ambient gap on time
drift_baseline <- function(r) { a <- r[is.na(r$id) & !is.na(r$gap), ]; g <- stats::aggregate(cbind(x, t) ~ gap, a, stats::median)
  if (nrow(g) < 3) return(NA_real_); unname(stats::coef(stats::lm(x ~ t, data = g))[2]) }

# ---- diagnostics -------------------------------------------------------------
diag_white <- function(r) stats::median(unlist(by_win(r, function(v) { d <- diff(v, differences = 2); if (length(d) < 6) return(NULL)
  stats::cor(d[-1], d[-length(d)]) })), na.rm = TRUE)                       # white noise -> -2/3
diag_leak <- function(r) sig_win_cent(r) / sig_win_d2(r)                    # >1.2: trend/curvature leaking into 1st diffs
diag_conc <- function(r) { m <- unlist(by_win(r, mean)); s <- sig_i_d2(r)  # elasticity of sigma_i with mean concentration
  ok <- is.finite(m) & is.finite(s) & s > 0 & m > 0; if (sum(ok) < 8) return(NA_real_)
  unname(stats::coef(stats::lm(log(s[ok]) ~ log(m[ok])))[2]) }
diag_zero <- function(r) mean(diff(r$x) == 0, na.rm = TRUE)                 # alternating-gas / quantisation
