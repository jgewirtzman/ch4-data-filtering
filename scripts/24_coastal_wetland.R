# ============================================================================
# 24_coastal_wetland.R — Coastal wetland closures (CO2, CH4, N2O) from raw records
#
# Mangrove-saltmarsh creek, Woogoompah Island, southern Moreton Bay, Australia;
# winter campaign FW5, 15-19 August 2024 (J. A. Rosentreter and colleagues, Southern
# Cross University). 36 chamber positions (6 habitats x 3 replicates x light/dark);
# CO2 and CH4 with an LI-7810 and N2O with an LI-7820, in separate closures on the same
# collars, at 1 Hz. Closure windows, chamber volume and temperature from the field sheet;
# timestamps as corrected for analyzer clock offset in the source workbook.
#
# Processing as for the stems: goFlux LM + HM fits over the window picked by hand in
# 24a_click_coastal_windows.R (fallback: field-sheet window after a dead band), with prec = campaign sigma and conf = 0.95; per-closure second-difference sigma
# (goFlux::empirical.prec) and campaign sigma = median per analyzer x field session;
# reference MDF = 1.96 sigma / t x flux term; datasheet MDF = datasheet precision / t x
# flux term.
#
# Loads:  data/local_full_sources/woogoompah/ (raw records and field sheet; local only)
# Writes: data/from_coastal_wetland_study/coastal_closures.csv (one row per closure x gas)
# ============================================================================

source("scripts/00_setup.R")
suppressMessages({library(readr); library(goFlux)})

src_dir <- file.path(data_dir, "local_full_sources", "woogoompah")
out_dir <- file.path(data_dir, "from_coastal_wetland_study")
if (!dir.exists(src_dir)) {
  message("Coastal wetland raw records not found (", src_dir, "); skipping")
} else {

dir.create(out_dir, showWarnings = FALSE)
P <- function(x) as.POSIXct(x, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
# Dead band after the field-sheet start (s): the traces show 2-5 min of headspace mixing after
# placement before the concentration change becomes steady
DEAD_BAND <- as.numeric(Sys.getenv("COASTAL_DEAD_BAND", 300))
# Windows picked by hand (24a_click_coastal_windows.R) replace the field-sheet window and
# dead band when present; closures without a picked window are then left out
click_file <- file.path(out_dir, "coastal_windows_clicked.csv")
clicked <- if (file.exists(click_file)) {
  read.csv(click_file) %>% mutate(start = P(start.time_corr), end = P(end.time_corr))
} else NULL
message(if (is.null(clicked)) paste0("Windows: field sheet, dead band ", DEAD_BAND, " s")
        else paste0("Windows: picked by hand (", nrow(clicked), " closures)"))

# --- field sheet: two blocks side by side (N2O; CO2 and CH4) -----------------
fs <- read_csv(file.path(src_dir, "field_sheet.csv"), col_names = FALSE, skip = 4,
               col_types = cols(.default = "c"))
block <- function(cols, analyzer) {
  x <- fs[, cols]
  names(x) <- c("id", "camp", "date", "start", "end", "chamber", "season", "tide", "habitat",
                "light", "area_m2", "vol_m3", "vsys_m3", "tempK")
  x %>% filter(!is.na(id)) %>%
    mutate(across(c(area_m2, vol_m3, vsys_m3, tempK), as.numeric),
           analyzer = analyzer, session = substr(date, 1, 10),
           t0 = P(paste(session, start)), t1 = P(paste(session, end)),
           habitat = tolower(sub("^carb", "crab", habitat)), light = tolower(light))
}
fld <- bind_rows(block(10:23, "LI-7820"), block(27:40, "LI-7810"))
# a closure that ends past midnight belongs to the next calendar day
fld <- fld %>% mutate(t1 = if_else(t1 < t0, t1 + 86400, t1))

# --- raw records (corrected timestamps) ---------------------------------------
rd <- function(f) { r <- read_csv(file.path(src_dir, f), col_names = FALSE, skip = 7,
                                  col_types = cols(.default = "c"))
                    r[!is.na(r$X1) & r$X1 == "DATA", ] }
a <- rd("licor_n2o.csv"); b <- rd("licor_co2_ch4.csv")
raw <- list(
  "LI-7820" = tibble(POSIX.time = P(a$X10), diag = a$X5, H2O_ppm = as.numeric(a$X11),
                     N2Odry_ppb = as.numeric(a$X12)),
  "LI-7810" = tibble(POSIX.time = P(b$X9), diag = b$X5, H2O_ppm = as.numeric(b$X10),
                     CO2dry_ppm = as.numeric(b$X11), CH4dry_ppb = as.numeric(b$X12)))
# readings with a non-zero analyzer diagnostic are dropped (10 inside the windows, all spikes)
raw <- lapply(raw, function(r) r[r$diag == "0", ])
rm(a, b)

# --- one window per closure ---------------------------------------------------
win <- bind_rows(lapply(seq_len(nrow(fld)), function(i) {
  f <- fld[i, ]; r <- raw[[f$analyzer]]
  lo <- f$t0 + DEAD_BAND; hi <- f$t1
  if (!is.null(clicked)) { k <- match(f$id, clicked$UniqueID)
    if (is.na(k)) return(NULL)
    lo <- clicked$start[k]; hi <- clicked$end[k] }
  w <- r[r$POSIX.time >= lo & r$POSIX.time <= hi, ]
  if (nrow(w) == 0) return(NULL)
  w %>% mutate(UniqueID = f$id, analyzer = f$analyzer, session = f$session, flag = 1L,
               Etime = as.numeric(POSIX.time - POSIX.time[1], units = "secs"),
               Area = f$area_m2 * 1e4, Vtot = (f$vol_m3 + f$vsys_m3) * 1e3,
               Tcham = f$tempK - 273.15, Pcham = 101.325)
}))
message("Closures with readings: ", length(unique(win$UniqueID)), " of ", nrow(fld))

# --- precision: per closure, then campaign sigma per analyzer x session -------
gases <- list("LI-7810" = c("CO2dry_ppm", "CH4dry_ppb"), "LI-7820" = "N2Odry_ppb")
spec  <- c(CO2dry_ppm = PREC_CO2_7810, CH4dry_ppb = PREC_CH4_7810, N2Odry_ppb = PREC_N2O_7820)
res <- list()
for (an in names(gases)) for (g in gases[[an]]) {
  d <- win %>% filter(analyzer == an, is.finite(.data[[g]]))
  sc <- d %>% group_by(UniqueID, session) %>%
    summarise(sigma_i = sigma_closure(.data[[g]], POSIX.time), .groups = "drop") %>%
    group_by(session) %>% mutate(sigma_c = median(sigma_i, na.rm = TRUE)) %>% ungroup()
  prec_col <- sub("(dry)?_(ppm|ppb)$", "_prec", g)
  d[[prec_col]] <- sc$sigma_c[match(d$UniqueID, sc$UniqueID)]
  fit <- suppressWarnings(goFlux(as.data.frame(d), g, H2O_col = "H2O_ppm", warn.length = 60, conf = 0.95))
  best <- suppressWarnings(best.flux(fit, criteria = c("MAE", "RMSE", "AICc", "SE", "g.factor", "kappa",
                                                       "MDF", "nb.obs", "intercept", "p-value"),
                                     g.limit = 2, k.ratio = 1, p.val = 0.05, warn.length = 60))
  # per closure: closure length; OLS slope p-value, residual lag-1 autocorrelation and the
  # regression interval half-width (Cowan et al. 2025) with the campaign sigma; and the
  # lag-1 correlation of consecutive second differences (-2/3 for white noise; Box S1)
  tl <- d %>% group_by(UniqueID) %>%
    summarise(n = n(), t_s = closure.time(Etime),
              lm_p = slope_fit_stats(.data[[g]], Etime)[["lm_p"]],
              resid_ac1 = slope_fit_stats(.data[[g]], Etime)[["ac1"]],
              sxx = sum((Etime - mean(Etime))^2),
              d2_ac1 = { d2 <- diff(.data[[g]], differences = 2); cor(d2[-1], d2[-length(d2)]) },
              .groups = "drop")
  out <- best %>% select(UniqueID, model, best.flux, LM.flux, LM.r2, flux.term) %>%
    left_join(sc, by = "UniqueID") %>% left_join(tl, by = "UniqueID") %>%
    mutate(gas = sub("dry_ppm|dry_ppb", "", g), analyzer = an,
           MDF_ref = qnorm(0.975) * sigma_c / t_s * flux.term,
           MDF_ds  = spec[[g]] / t_s * flux.term,
           reg_half_95 = qt(0.975, n - 2) * sigma_c / sqrt(sxx) * flux.term) %>% select(-sxx)
  res[[g]] <- out
}
cl <- bind_rows(res) %>%
  left_join(fld %>% select(UniqueID = id, habitat, light, chamber), by = "UniqueID") %>%
  mutate(unit = if_else(gas == "CO2", "umol m-2 s-1", "nmol m-2 s-1")) %>%
  relocate(UniqueID, gas, analyzer, session, habitat, light, chamber)
cl$window_source <- if (is.null(clicked)) paste0("field sheet + ", DEAD_BAND, " s dead band") else "picked by hand"
write.csv(cl, file.path(out_dir, "coastal_closures.csv"), row.names = FALSE)
message("Wrote ", nrow(cl), " closure x gas rows")
print(cl %>% group_by(gas) %>% summarise(n = n(), HM = mean(model == "HM"), sigma = median(sigma_c),
        spec = spec[paste0(gas[1], ifelse(gas[1] == "CO2", "dry_ppm", "dry_ppb"))],
        t_med = median(t_s), below_ref = mean(abs(best.flux) < MDF_ref), r2_le_05 = mean(LM.r2 <= 0.5),
        mean_flux = mean(best.flux), median_flux = median(best.flux), neg = mean(best.flux < 0)))
}
