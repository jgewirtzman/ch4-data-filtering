# ============================================================================
# 21_logging_interval.R — What the logging interval does to realized precision
#
# The 1 Hz stem records are re-logged at 2, 5 and 10 s in two ways an analyzer
# can log: averaging the readings in each interval, or keeping one snapshot per
# interval. sigma is re-estimated from each re-logged closure (same estimator,
# campaign median per analyzer). Averaging should scale sigma by 1/sqrt(dt), so a
# 1 s datasheet value must be rescaled before comparison; snapshot logging should
# not change sigma, so rescaling would then be wrong.
#
# Loads: outputs/02b_manID.RData (closure windows with 1 Hz readings)
# Writes: outputs/21_logging_interval.RData
# ============================================================================

source("scripts/00_setup.R")
set.seed(21)
load(file.path(output_dir, "02b_manID.RData"))      # manID

win <- manID %>% filter(flag == 1, is.finite(CH4dry_ppb)) %>% arrange(UniqueID, POSIX.time)
ids <- win %>% group_by(UniqueID, analyzer) %>% summarise(n = n(), .groups = "drop") %>%
  filter(n >= 150) %>% group_by(analyzer) %>% slice_sample(n = 300) %>% ungroup()

relog <- function(x, dt, how) {
  if (dt == 1) return(x)
  if (how == "average") {
    k <- floor(length(x) / dt); colMeans(matrix(x[seq_len(k * dt)], nrow = dt))
  } else x[seq(1, length(x), by = dt)]
}
res <- list()
for (i in seq_len(nrow(ids))) {
  x <- win$CH4dry_ppb[win$UniqueID == ids$UniqueID[i]]
  for (dt in c(1, 2, 5, 10)) for (how in c("average", "snapshot")) {
    if (dt == 1 && how == "snapshot") next
    res[[length(res) + 1]] <- data.frame(analyzer = ids$analyzer[i], id = ids$UniqueID[i], dt = dt,
                                         how = if (dt == 1) "1 s record" else how, sigma = sigma_closure(relog(x, dt, how)))   # evenly spaced by construction: no timestamps needed
  }
}
res <- do.call(rbind, res)
res$analyzer <- ifelse(grepl("7810", res$analyzer), "Precise analyzer", "Noisier analyzer")
logsum <- res %>% group_by(analyzer, dt, how) %>% summarise(sigma = median(sigma), .groups = "drop")
base <- logsum %>% filter(dt == 1) %>% select(analyzer, sigma1 = sigma)
logsum <- logsum %>% left_join(base, by = "analyzer") %>%
  mutate(expected_avg = sigma1 / sqrt(dt), ratio_to_avg_expectation = sigma / expected_avg)
print(as.data.frame(logsum), digits = 3)
spec <- c("Noisier analyzer" = PREC_CH4_LGR, "Precise analyzer" = PREC_CH4_7810)
save(res, logsum, spec, file = file.path(output_dir, "21_logging_interval.RData"))
message("Saved 21_logging_interval")
