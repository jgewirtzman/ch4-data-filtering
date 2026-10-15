# ============================================================================
# 22_transform_inference.R — What the transformation does to inference
#
# The same questions asked of the same fluxes under each transformation:
#   raw; log10 (positives only, as usually done); log10(x + 1); arcsinh(x);
#   arcsinh(x / s) with s the dataset's median reference MDF; and ranks.
#
# 1. Real data (stems, noisier analyzer, 2023-2024):
#    (a) species and stand (upland vs wetland) differences, tree as a random effect;
#    (b) the temperature response, tree as a random effect.
# 2. Known-truth simulation (closures simulated as in 18, noisier analyzer):
#    (a) two groups with no difference (false-positive rate) and a 2-fold difference (power);
#    (b) the temperature response with true Q10 = 2.5: recovered Q10 per transformation.
#
# Loads:  outputs/10_seasonal_confounding.RData (df_gla)
# Writes: outputs/22_transform_inference.RData, outputs/22_transform_inference_tables.csv
# ============================================================================

source("scripts/00_setup.R")
suppressPackageStartupMessages(library(lme4))
set.seed(22)

e <- new.env(); load(file.path(output_dir, "10_seasonal_confounding.RData"), envir = e)
d <- e$df_gla %>% filter(is.finite(CH4_flux_nmolpm2ps), is.finite(airt)) %>%
  transmute(flux = CH4_flux_nmolpm2ps, species = species_label, stand = location, tree = factor(Tree), airt,
            mdf = CH4_MDF_ref)
s_mdf <- median(d$mdf, na.rm = TRUE)

TRANS <- list(
  "Raw"                 = list(f = function(x) x,                   keep = function(x) rep(TRUE, length(x)), q10 = NA),
  "log10, positives"    = list(f = function(x) log10(pmax(x, 1e-12)), keep = function(x) x > 0,           q10 = log(10)),
  "log10(x + 1)"        = list(f = function(x) log10(x + 1),          keep = function(x) x > -1,          q10 = log(10)),
  "arcsinh(x)"          = list(f = function(x) asinh(x),              keep = function(x) rep(TRUE, length(x)), q10 = 1),
  "arcsinh(x / MDF)"    = list(f = function(x) asinh(x / s_mdf),      keep = function(x) rep(TRUE, length(x)), q10 = 1)
)
# Implied Q10 from a slope b per degree on a log-like scale with base factor k: exp(10 * k * b)
q10_of <- function(b, k) if (is.na(k)) NA_real_ else exp(10 * k * b)
diag_of <- function(r) c(skew = mean(r^3) / sd(r)^3, kurt = mean(r^4) / sd(r)^4)

# ── 1. Real data ─────────────────────────────────────────────────────────────
lrt_p <- function(full, null) anova(null, full, refit = TRUE)[2, "Pr(>Chisq)"]
real <- list()
for (nm in names(TRANS)) {
  tr <- TRANS[[nm]]; k <- tr$keep(d$flux); dd <- d[k, ]; dd$y <- tr$f(dd$flux)
  m0 <- lmer(y ~ 1 + (1 | tree), dd, REML = FALSE)
  m_sp <- lmer(y ~ species + (1 | tree), dd, REML = FALSE)
  m_st <- lmer(y ~ stand + (1 | tree), dd, REML = FALSE)
  m_t  <- lmer(y ~ airt + (1 | tree), dd, REML = FALSE)
  b_t <- fixef(m_t)[["airt"]]
  dg <- diag_of(resid(m_t))
  real[[nm]] <- data.frame(method = nm, n = nrow(dd), pct_used = 100 * nrow(dd) / nrow(d),
                           p_species = lrt_p(m_sp, m0), p_stand = lrt_p(m_st, m0), p_temp = lrt_p(m_t, m0),
                           q10 = q10_of(b_t, tr$q10), resid_skew = dg[["skew"]], resid_kurt = dg[["kurt"]])
}
# Ranks: rank-transform the fluxes, then the same mixed models (Conover and Iman 1981)
dr <- d; dr$y <- rank(dr$flux)
r0 <- lmer(y ~ 1 + (1 | tree), dr, REML = FALSE)
real[["Ranks"]] <- data.frame(method = "Ranks", n = nrow(d), pct_used = 100,
                              p_species = lrt_p(lmer(y ~ species + (1 | tree), dr, REML = FALSE), r0),
                              p_stand = lrt_p(lmer(y ~ stand + (1 | tree), dr, REML = FALSE), r0),
                              p_temp = lrt_p(lmer(y ~ airt + (1 | tree), dr, REML = FALSE), r0),
                              q10 = NA, resid_skew = NA, resid_kurt = NA)
real <- do.call(rbind, real); rownames(real) <- NULL
# Raw-scale exponential model, flux = a * Q10^((T - 15) / 10), fitted to the untransformed fluxes
q10_nls <- function(y, temp) { m <- try(nls(y ~ a * q^((temp - 15) / 10), start = list(a = max(mean(y), 0.01), q = 2)), silent = TRUE)
  if (inherits(m, "try-error")) NA_real_ else coef(m)[["q"]] }
real$q10[real$method == "Raw"] <- q10_nls(d$flux, d$airt)
# Reference Q10 on the raw scale: log of 5 C-bin arithmetic means against temperature
bins <- d %>% mutate(b = cut(airt, seq(-10, 35, 5))) %>% group_by(b) %>%
  summarise(t = mean(airt), m = mean(flux), n = n(), .groups = "drop") %>% filter(n >= 20, m > 0)
q10_ref <- exp(10 * coef(lm(log(m) ~ t, bins, weights = n))[["t"]])
# Species ordering on each scale vs arithmetic means
ord_raw <- names(sort(tapply(d$flux, d$species, mean)))
ord <- sapply(names(TRANS), function(nm) { tr <- TRANS[[nm]]; k <- tr$keep(d$flux)
  paste(names(sort(tapply(tr$f(d$flux[k]), d$species[k], mean))), collapse = " < ") })
# "How much": each species' arithmetic mean flux and the back-transformed means on each scale
sp_means <- d %>% group_by(species) %>%
  summarise(`Arithmetic mean` = mean(flux),
            `Back-transformed, log10 (positives)` = 10^mean(log10(flux[flux > 0])),
            `Back-transformed, arcsinh(flux / MDF)` = s_mdf * sinh(mean(asinh(flux / s_mdf))),
            n = n(), .groups = "drop")
print(as.data.frame(sp_means), digits = 3)
print(real, digits = 3); message("Reference Q10 (arithmetic bin means): ", round(q10_ref, 2)); print(ord)

# ── 1b. The same workflows on 20 everyday comparisons ────────────────────────
dd0 <- d %>% mutate(season = ifelse(grepl("Growing", e$df_gla$season[match(seq_len(nrow(d)), seq_len(nrow(d)))]), "Growing", "Dormant"))
dd0$season <- ifelse(grepl("Growing", (e$df_gla %>% filter(is.finite(CH4_flux_nmolpm2ps), is.finite(airt)))$season), "Growing", "Dormant")
TRK <- list("Raw" = function(x) x, "log10, positives" = function(x) log10(pmax(x, 1e-12)), "log10(x + 1)" = function(x) log10(x + 1),
            "arcsinh(x)" = function(x) asinh(x), "arcsinh(x / MDF)" = function(x) asinh(x / s_mdf), "Ranks" = function(x) rank(x))
test_all <- function(dd, term) sapply(names(TRK), function(nm) {
  k <- if (nm == "log10, positives") dd$flux > 0 else rep(TRUE, nrow(dd)); x <- dd[k, ]; x$y <- TRK[[nm]](x$flux)
  m1 <- try(suppressMessages(lmer(as.formula(paste("y ~", term, "+ (1 | tree)")), x, REML = FALSE)), silent = TRUE)
  m0 <- try(suppressMessages(lmer(y ~ 1 + (1 | tree), x, REML = FALSE)), silent = TRUE)
  if (inherits(m1, "try-error") || inherits(m0, "try-error")) NA else suppressMessages(anova(m0, m1))[2, "Pr(>Chisq)"] })
cmp <- list(); spp <- sort(unique(dd0$species))
for (i in 1:3) for (j in (i + 1):4) cmp[[paste(spp[i], "vs", spp[j])]] <- test_all(dd0 %>% filter(species %in% spp[c(i, j)]), "species")
for (x in spp) { z <- dd0 %>% filter(species == x); if (length(unique(z$stand)) > 1) cmp[[paste0(x, ": upland vs wetland")]] <- test_all(z, "stand") }
for (x in spp) cmp[[paste0(x, ": growing vs dormant")]] <- test_all(dd0 %>% filter(species == x), "season")
for (x in c("Upland", "Wetland")) cmp[[paste(x, "stems: growing vs dormant")]] <- test_all(dd0 %>% filter(stand == x), "season")
for (x in spp) cmp[[paste0(x, ": temperature")]] <- test_all(dd0 %>% filter(species == x), "airt")
for (x in c("Upland", "Wetland")) cmp[[paste(x, "stems: temperature")]] <- test_all(dd0 %>% filter(stand == x), "airt")
cmp_p <- do.call(rbind, cmp)
print(signif(cmp_p, 2))
# Worked example: A. rubrum vs T. canadensis, end to end
ex <- dd0 %>% filter(species %in% c("A. rubrum", "T. canadensis")) %>% mutate(species = factor(species, levels = c("T. canadensis", "A. rubrum")))
back <- list("Raw" = function(m) m, "log10, positives" = function(m) 10^m, "log10(x + 1)" = function(m) 10^m - 1,
             "arcsinh(x)" = function(m) sinh(m), "arcsinh(x / MDF)" = function(m) s_mdf * sinh(m), "Ranks" = function(m) NA)
wf <- do.call(rbind, lapply(names(TRK), function(nm) {
  k <- if (nm == "log10, positives") ex$flux > 0 else rep(TRUE, nrow(ex)); x <- ex[k, ]; x$y <- TRK[[nm]](x$flux)
  m1 <- lmer(y ~ species + (1 | tree), x, REML = FALSE); m0 <- lmer(y ~ 1 + (1 | tree), x, REML = FALSE)
  a <- fixef(m1)[["(Intercept)"]]; b <- fixef(m1)[["speciesA. rubrum"]]
  data.frame(method = nm, used = sprintf("%d of %d", sum(k), nrow(ex)), p = anova(m0, m1)[2, "Pr(>Chisq)"],
             ratio = if (nm == "Raw") mean(ex$flux[ex$species == "A. rubrum"]) / mean(ex$flux[ex$species == "T. canadensis"]) else back[[nm]](a + b) / back[[nm]](a),
             resid = I(list(resid(m1))))
}))
print(wf[, c("method", "used", "p", "ratio")], digits = 3)

# ── 2. Known-truth simulation ────────────────────────────────────────────────
FLUX_TERM <- 2.2; N_OBS <- 180; SIG <- 3.6; Z95 <- qnorm(0.975)
tt <- seq_len(N_OBS) - 1; tc <- tt - mean(tt); Sxx <- sum(tc^2)
temp_m <- 11 + 14 * sin(2 * pi * ((1:12) - 4.5) / 12)
true_flux <- function(n, temp, mult = 1) {
  f <- rlnorm(n, log(0.12 * 2.5^((temp - 15) / 10)), 1.3) * mult
  up <- runif(n) < 0.10; f[up] <- -abs(rnorm(sum(up), 0.02, 0.01)); f
}
measure <- function(f) {                      # fitted flux of each closure under white noise
  x <- outer(f / FLUX_TERM, tt) + matrix(rnorm(length(f) * N_OBS, 0, SIG), length(f))
  as.vector(x %*% tc) / Sxx * FLUX_TERM
}
s_sim <- Z95 * SIG / N_OBS * FLUX_TERM        # the reference MDF of these closures
TRANS_SIM <- TRANS; TRANS_SIM[["arcsinh(x / MDF)"]]$f <- function(x) asinh(x / s_sim)
true_flux_up <- function(n, temp, p_up) {     # same emission distribution, a larger share of small uptake
  f <- rlnorm(n, log(0.12 * 2.5^((temp - 15) / 10)), 1.3)
  up <- runif(n) < p_up; f[up] <- -abs(rnorm(sum(up), 0.02, 0.01)); f
}
one_group_test <- function(mult, n_month = 10, p_up2 = NULL) {
  g1 <- unlist(lapply(1:12, function(m) measure(true_flux(n_month, temp_m[m]))))
  g2 <- unlist(lapply(1:12, function(m) measure(if (is.null(p_up2)) true_flux(n_month, temp_m[m], mult) else true_flux_up(n_month, temp_m[m], p_up2))))
  out <- sapply(names(TRANS_SIM), function(nm) { tr <- TRANS_SIM[[nm]]
    a <- g1[tr$keep(g1)]; b <- g2[tr$keep(g2)]; t.test(tr$f(a), tr$f(b))$p.value })
  c(out, Ranks = wilcox.test(g1, g2)$p.value)
}
N_REP <- 500
p_null <- t(replicate(N_REP, one_group_test(1)))
p_alt  <- t(replicate(N_REP, one_group_test(2)))
p_up   <- t(replicate(N_REP, one_group_test(1, p_up2 = 0.40)))   # 40 % vs 10 % uptake: true mean 33 % lower
grp_sim <- data.frame(method = colnames(p_null), false_pos = 100 * colMeans(p_null < 0.05),
                      power_2fold = 100 * colMeans(p_alt < 0.05), power_uptake_share = 100 * colMeans(p_up < 0.05))
# Whether and how much: two groups, true means 2-fold apart; each workflow's test and reported ratio
BACK <- list("Raw" = function(m) m, "log10, positives" = function(m) 10^m, "log10(x + 1)" = function(m) 10^m - 1,
             "arcsinh(x)" = function(m) sinh(m), "arcsinh(x / MDF)" = function(m) s_sim * sinh(m))
true_flux_hot <- function(n, temp) {           # same median, wider spread: mean 2-fold, typical closure unchanged
  sd2 <- sqrt(1.3^2 + 2 * log(2))
  f <- rlnorm(n, log(0.12 * 2.5^((temp - 15) / 10)), sd2)
  up <- runif(n) < 0.10; f[up] <- -abs(rnorm(sum(up), 0.02, 0.01)); f
}
two_groups <- function(n_month = 10, hot = FALSE) {
  g1 <- unlist(lapply(1:12, function(m) measure(true_flux(n_month, temp_m[m]))))
  g2 <- unlist(lapply(1:12, function(m) measure(if (hot) true_flux_hot(n_month, temp_m[m]) else true_flux(n_month, temp_m[m], 2))))
  out <- lapply(names(TRANS_SIM), function(nm) { tr <- TRANS_SIM[[nm]]; a <- g1[tr$keep(g1)]; b <- g2[tr$keep(g2)]
    ya <- tr$f(a); yb <- tr$f(b)
    c(p = t.test(ya, yb)$p.value, ratio = if (nm == "Raw") mean(b) / mean(a) else BACK[[nm]](mean(yb)) / BACK[[nm]](mean(ya))) })
  names(out) <- names(TRANS_SIM)
  out[["Ranks"]] <- c(p = wilcox.test(g1, g2)$p.value, ratio = NA)
  out[["Recommended"]] <- c(p = out[["arcsinh(x / MDF)"]][["p"]], ratio = mean(g2) / mean(g1))
  out
}
summ_wh <- function(tg, scen) do.call(rbind, lapply(names(tg[[1]]), function(nm) {
  pp <- sapply(tg, function(z) z[[nm]][["p"]]); rr <- sapply(tg, function(z) z[[nm]][["ratio"]])
  data.frame(scenario = scen, method = nm, detected = 100 * mean(pp < 0.05), ratio_med = median(rr),
             ratio_lo = quantile(rr, 0.05, na.rm = TRUE), ratio_hi = quantile(rr, 0.95, na.rm = TRUE)) }))
wh <- rbind(summ_wh(replicate(N_REP, two_groups(), simplify = FALSE), "Every flux 2\u00d7 larger"),
            summ_wh(replicate(N_REP, two_groups(hot = TRUE), simplify = FALSE), "More hot spots (mean 2\u00d7, median unchanged)"))
rownames(wh) <- NULL; print(wh, digits = 3)
q10_sim <- t(replicate(200, {
  temp <- rep(temp_m, each = 40); f <- measure(unlist(lapply(temp_m, function(tm) true_flux(40, tm))))
  out <- sapply(names(TRANS_SIM), function(nm) { tr <- TRANS_SIM[[nm]]; k <- tr$keep(f)
    q10_of(coef(lm(tr$f(f[k]) ~ temp[k]))[[2]], tr$q10) })
  out[["Raw"]] <- q10_nls(f, temp); out
}))
q10_sum <- data.frame(method = colnames(q10_sim), q10_median = apply(q10_sim, 2, median, na.rm = TRUE),
                      q10_lo = apply(q10_sim, 2, quantile, 0.05, na.rm = TRUE), q10_hi = apply(q10_sim, 2, quantile, 0.95, na.rm = TRUE))
print(grp_sim, digits = 3); print(q10_sum, digits = 3)

tab <- real %>% left_join(grp_sim, by = "method") %>% left_join(q10_sum, by = "method")
write.csv(tab, file.path(output_dir, "22_transform_inference_tables.csv"), row.names = FALSE)
save(real, sp_means, cmp_p, wf, wh, q10_ref, ord, ord_raw, grp_sim, q10_sum, tab, s_mdf, s_sim, file = file.path(output_dir, "22_transform_inference.RData"))
message("Saved 22_transform_inference")
