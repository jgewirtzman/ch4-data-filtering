# ============================================================================
# 20_generality.R — The same processing chain across gases, surfaces and analyzers
#
# Datasets (rows of the generality table):
#   stems, noisier analyzer (CH4) and precise analyzer (CH4)   — 03 output
#   soils along an upland-wetland gradient (CH4): upland + intermediate positions,
#     and wetland positions; sigma and MDF computed here from the raw closure windows
#     (Gewirtzman (2026), doi:10.5281/zenodo.18779715) with the paper's convention
#   agricultural soil (CH4 and N2O, precise analyzers)           — data/from_agsoil_study (de-identified)
# Metrics per dataset: realized sigma vs datasheet; % below the reference and the
# datasheet MDF; % removed by R2 <= 0.5; change in the magnitude of the mean and the
# median when below-MDF or low-R2 closures are removed; % of values <= 0 (lost to log).
#
# Loads:  outputs/03_mdf_results.RData, outputs/10_seasonal_confounding.RData,
#         data/from_ymf_soil_zenodo/ (soil closure windows, flagged readings only, or the full
#         source file named in SOIL_WINDOWS; goFlux fits), data/from_agsoil_study/
# Writes: outputs/20_generality.RData, outputs/20_generality_table.csv
# ============================================================================

source("scripts/00_setup.R")
suppressPackageStartupMessages(library(data.table))

# Second-difference sigma of one closure, over consecutive readings (logging gaps not
# split, unlike sigma_closure() in 00_setup.R; switching changes the wetland row of
# Table S4 slightly, so the published values keep this form)
sig2 <- function(x) { d2 <- diff(x, differences = 2); 1.4826 * median(abs(d2 - median(d2))) / sqrt(6) }
mag_change <- function(v, keep) {           # % change in |statistic| when only `keep` is retained
  c(mean = 100 * (abs(mean(v[keep])) / abs(mean(v)) - 1),
    median = 100 * (abs(median(v[keep])) / abs(median(v)) - 1))
}
summ <- function(name, gas, surface, analyzer, flux, r2, sigma, spec, mdf_ref, mdf_ds) {
  ok <- is.finite(flux)
  flux <- flux[ok]; r2 <- r2[ok]; mdf_ref <- mdf_ref[ok]; mdf_ds <- mdf_ds[ok]
  below <- abs(flux) < mdf_ref
  rm_mdf <- mag_change(flux, !below); rm_r2 <- mag_change(flux, is.finite(r2) & r2 > 0.5)
  data.frame(dataset = name, gas = gas, surface = surface, analyzer = analyzer, n = length(flux),
             sigma = sigma, spec = spec, sigma_ratio = sigma / spec,
             pct_below_ref = 100 * mean(below), pct_below_ds = 100 * mean(abs(flux) < mdf_ds, na.rm = TRUE),
             pct_r2_removed = 100 * mean(!(is.finite(r2) & r2 > 0.5)),
             mean_change_mdf = rm_mdf[["mean"]], median_change_mdf = rm_mdf[["median"]],
             mean_change_r2 = rm_r2[["mean"]], median_change_r2 = rm_r2[["median"]],
             pct_nonpositive = 100 * mean(flux <= 0), pct_negative = 100 * mean(flux < 0))
}

rows <- list()
# ── Stems ────────────────────────────────────────────────────────────────────
st <- local({ e <- new.env(); load(file.path(output_dir, "03_mdf_results.RData"), envir = e); e$df })
for (a in c("Noisier analyzer", "Precise analyzer")) {
  d <- st[st$inst_label == a & st$t_src == "measured", ]
  rows[[a]] <- summ(paste("Stems,", tolower(a)), "CH4", "Forest stems", a, d$CH4_flux_nmolpm2ps, d$CH4_r2,
                    median(d$sigma_camp_CH4, na.rm = TRUE), ifelse(a == "Noisier analyzer", PREC_CH4_LGR, PREC_CH4_7810),
                    d$CH4_MDF_ref, d$CH4_MDF_manufacturer)
}

# ── Soils along an upland–wetland gradient ───────────────────────────────────
# Soil closure windows of the published soil-gradient dataset (Gewirtzman (2026),
# doi:10.5281/zenodo.18779715): by default the copy in data/from_ymf_soil_zenodo/, reduced to the
# readings inside the fitting windows (flag == 1) and the four columns used here. SOIL_WINDOWS
# can name the complete source file (lgr_manual_identification_results_soil.csv, 75 MB) instead;
# the results are the same.
win_file <- Sys.getenv("SOIL_WINDOWS", file.path(soil_dir, "lgr_manual_identification_results_soil_flag1.csv.gz"))
if (!file.exists(win_file)) stop("Soil closure windows not found: ", win_file, " (see data/from_ymf_soil_zenodo/README.md)")
win <- fread(win_file)
win <- win[flag == 1 & is.finite(CH4dry_ppb)]
win[, POSIX.time := as.POSIXct(POSIX.time, tz = "UTC")]
setorder(win, UniqueID, POSIX.time)
cl <- win[, .(sigma_i = sig2(CH4dry_ppb), t = as.numeric(diff(range(POSIX.time)), units = "secs") + 1,
              day = as.character(as.Date(POSIX.time[1]))), by = UniqueID]
cl[, sigma_c := median(sigma_i, na.rm = TRUE), by = day]
gf <- read.csv(file.path(soil_dir, "CH4_best_flux_lgr_results_soil.csv"))
cl <- merge(cl, gf[, c("UniqueID", "flux.term", "best.flux")], by = "UniqueID")
cl[, `:=`(mdf_ref = Z95 * sigma_c / t * flux.term, mdf_ds = PREC_CH4_LGR / t * flux.term)]
soil <- local({ e <- new.env(); load(file.path(output_dir, "10_seasonal_confounding.RData"), envir = e); e$soil })
soil$UniqueID <- soil$unique_id
sm <- merge(as.data.frame(cl), soil[, c("UniqueID", "flux", "r2", "landscape_position")], by = "UniqueID")
message("Soil closures with sigma and MDF: ", nrow(sm), " of ", nrow(soil))
for (pos in list(c("U", "I"), c("WD", "WS"))) {
  d <- sm[sm$landscape_position %in% pos, ]
  nm <- if ("U" %in% pos) "Soils, upland" else "Soils, wetland"
  rows[[nm]] <- summ(nm, "CH4", if ("U" %in% pos) "Upland soil" else "Wetland soil", "Noisier analyzer",
                     d$flux, d$r2, median(d$sigma_c), PREC_CH4_LGR, d$mdf_ref, d$mdf_ds)
}
soil_sigma <- sm

# ── Agricultural soil (de-identified) ────────────────────────────────────────
ag <- read.csv(file.path(agsoil_dir, "agsoil_closures_deidentified.csv"))
ag <- ag[!duplicated(ag[, c("CH4_flux", "N2O_flux", "CO2_flux")]), ]
spec_ag <- c(CH4 = PREC_CH4_7810, N2O = PREC_N2O_7820)   # LI-7810 CH4, LI-7820 N2O, 1 s
for (g in c("CH4", "N2O")) {
  f <- ag[[paste0(g, "_flux")]]; s <- ag[[paste0(g, "_sigma")]]; m <- ag[[paste0(g, "_MDF")]]
  ds <- m * spec_ag[[g]] / (Z95 * s)          # datasheet MDF, k = 1: same t and flux term
  rows[[paste("Ag", g)]] <- summ(paste("Agricultural soil,", g), g, "Agricultural soil", "Precise analyzer",
                                f, ag[[paste0(g, "_r2")]], median(s, na.rm = TRUE), spec_ag[[g]], m, ds)
}

gen <- do.call(rbind, rows); rownames(gen) <- NULL
print(gen %>% mutate(across(where(is.numeric), ~ round(.x, 2))), width = 200)
write.csv(gen, file.path(output_dir, "20_generality_table.csv"), row.names = FALSE)
save(gen, soil_sigma, file = file.path(output_dir, "20_generality.RData"))
message("Saved 20_generality")
