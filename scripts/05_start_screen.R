# ============================================================================
# 05_start_screen.R — Validity screen: apparent uptake from a headspace that
#                     starts above the ambient air
#
# A closure that begins with its headspace far above the surrounding air (carry-over,
# breath, a bubble or disturbance at placement) relaxes back toward ambient and registers
# as uptake (main-text Fig. 4e). Such a closure is invalid. Criterion, the same for every
# dataset: a detected uptake (flux <= -MDF) whose concentration at the start of the fitted
# window is more than 1.5 times the air before the chamber was placed. (A large emission
# can legitimately raise the headspace that far during the dead band, so the screen is
# applied to uptake only.) Subtler elevated starts, such as Fig. 4e itself, fall below this
# margin; the trace screens of Sect. S3 flag them for review.
# The air before placement is the median of the first 10 readings of the extracted record,
# or the day's median of that value where the record itself starts elevated (carry-over from
# the previous closure) (stems, soils), the 2 min before the field-sheet start (coastal wetland), or the reading
# at chamber closure (agricultural soil, whose record starts at closure).
#
# Loads:  outputs/02b_manID.RData, 03 inputs (stem flux, MDF, sigma: outputs/02b_goflux_fluxes.RData
#         via 03), data/external/lgr_manual_identification_results_soil.csv (soil records, local),
#         data/from_ymf_soil_zenodo/ (soil fits), data/local_full_sources/woogoompah/ (coastal,
#         local), the agricultural raw traces (local), and the closure tables of 20 and 24
# Writes: data/invalid_closures.csv (dataset, closure, air before placement, window start,
#         flux, MDF; tracked so that downstream scripts run without the raw records)
# Run after 03 (stems) and 24 (coastal); 03, 10, 20 and 25 drop the closures listed.
# ============================================================================

source("scripts/00_setup.R")
suppressMessages({library(data.table); library(readr)})
out_file <- file.path(data_dir, "invalid_closures.csv")
res <- list()
flag <- function(d) d %>% filter(flux <= -mdf, start > 1.5 * pre)

# ── Stems ─────────────────────────────────────────────────────────────────────
if (file.exists(file.path(output_dir, "02b_manID.RData"))) {
  m <- new.env(); load(file.path(output_dir, "02b_manID.RData"), envir = m)
  e <- new.env(); load(file.path(output_dir, "03_mdf_results.RData"), envir = e)
  s <- m$manID %>% filter(is.finite(CH4dry_ppb)) %>% arrange(UniqueID, POSIX.time) %>% group_by(UniqueID) %>%
    summarise(pre = median(head(CH4dry_ppb, 10)), start = CH4dry_ppb[flag == 1][1], day = as.Date(min(POSIX.time)),
              .groups = "drop") %>%
    group_by(day) %>% mutate(pre = pmin(pre, median(pre))) %>% ungroup() %>%
    inner_join(e$df_unscreened %>% transmute(UniqueID = stem_id, flux = CH4_flux_nmolpm2ps, mdf = CH4_MDF_ref,
                                  sigma = sigma_camp_CH4), by = "UniqueID")
  res$stems <- flag(s) %>% mutate(dataset = "Stems", gas = "CH4")
}

# ── Soils along the upland-wetland gradient ───────────────────────────────────
soil_full <- Sys.getenv("SOIL_FULL", file.path(data_dir, "external", "lgr_manual_identification_results_soil.csv"))
if (file.exists(soil_full)) {
  w <- fread(soil_full)[is.finite(CH4dry_ppb)]
  w[, POSIX.time := as.POSIXct(POSIX.time, tz = "UTC")]; setorder(w, UniqueID, POSIX.time)
  st <- w[, .(pre = median(head(CH4dry_ppb, 10)), start = CH4dry_ppb[flag == 1][1], day = as.Date(min(POSIX.time))), by = UniqueID]
  st[, pre := pmin(pre, median(pre)), by = day]; st[, day := NULL]
  win <- w[flag == 1]
  sg <- win[, .(sigma_i = sigma_closure(CH4dry_ppb, POSIX.time),           # as in 20_generality.R
                t = goFlux::closure.time(as.numeric(POSIX.time) - as.numeric(min(POSIX.time))),
                day = as.character(as.Date(POSIX.time[1]))), by = UniqueID]
  sg[, sigma := median(sigma_i, na.rm = TRUE), by = day]
  gf <- read.csv(file.path(soil_dir, "CH4_best_flux_lgr_results_soil.csv"))
  s <- merge(merge(st, sg, by = "UniqueID"), gf[, c("UniqueID", "flux.term", "best.flux")], by = "UniqueID")
  s <- as.data.frame(s) %>% transmute(UniqueID, pre, start, sigma, flux = best.flux, mdf = Z95 * sigma / t * flux.term)
  res$soils <- flag(s) %>% mutate(dataset = "Soils", gas = "CH4")
}

# ── Coastal wetland ───────────────────────────────────────────────────────────
cw_dir <- file.path(data_dir, "local_full_sources", "woogoompah")
cw_file <- file.path(data_dir, "from_coastal_wetland_study", "coastal_closures.csv")
if (dir.exists(cw_dir) && file.exists(cw_file)) {
  P <- function(x) as.POSIXct(x, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
  fs <- read_csv(file.path(cw_dir, "field_sheet.csv"), col_names = FALSE, skip = 4, col_types = cols(.default = "c"))
  blk <- function(cols, an) { x <- fs[, cols]; names(x)[c(1, 3, 4)] <- c("id", "date", "start")
    x %>% filter(!is.na(id)) %>% transmute(UniqueID = id, analyzer = an, t0 = P(paste(substr(date, 1, 10), start))) }
  aux <- bind_rows(blk(10:23, "LI-7820"), blk(27:40, "LI-7810"))
  rd <- function(f) { r <- read_csv(file.path(cw_dir, f), col_names = FALSE, skip = 7, col_types = cols(.default = "c"))
    r[!is.na(r$X1) & r$X1 == "DATA" & r$X5 == "0", ] }
  a <- rd("licor_n2o.csv"); b <- rd("licor_co2_ch4.csv")
  raw <- list("LI-7820" = list(N2O = data.frame(t = P(a$X10), g = as.numeric(a$X12))),
              "LI-7810" = list(CH4 = data.frame(t = P(b$X9), g = as.numeric(b$X12))))
  clk <- read.csv(file.path(data_dir, "from_coastal_wetland_study", "coastal_windows_clicked.csv"))
  cw <- read.csv(cw_file) %>% filter(gas %in% c("CH4", "N2O"))
  s <- cw %>% rowwise() %>% mutate(
    r = list(raw[[analyzer]][[gas]]), t0 = aux$t0[match(UniqueID, aux$UniqueID)],
    ws = P(clk$start.time_corr[match(UniqueID, clk$UniqueID)]),
    pre = median(r$g[r$t >= t0 - 120 & r$t < t0]), start = r$g[which.min(abs(r$t - ws))]) %>%
    ungroup() %>% transmute(UniqueID, gas, pre, start, sigma = sigma_c, flux = best.flux, mdf = MDF_ref)
  res$coastal <- flag(s) %>% mutate(dataset = "Coastal wetland")
}

# ── Agricultural soil (record starts at chamber closure; CH4 and N2O) ───────────
ag_raw <- Sys.getenv("AGSOIL_TRACES",
  "/Users/jongewirtzman/My Drive/Research/Churchtown Dairy/churchtown-dairy/data/intermediate/flux_raw_traces.rds")
ag_fit <- file.path(dirname(ag_raw), "flux_goflux.csv")
if (file.exists(ag_raw) && file.exists(ag_fit)) {
  x <- readRDS(ag_raw)$ch; fit <- read.csv(ag_fit)
  res$agsoil <- bind_rows(lapply(c("CH4", "N2O"), function(g) {
    col <- paste0(g, "dry_ppb")
    st <- x %>% filter(is.finite(.data[[col]])) %>% arrange(chamID, POSIX.time) %>% group_by(chamID) %>%
      summarise(pre = .data[[col]][1], start = .data[[col]][Etime >= 0][1], .groups = "drop")
    st %>% inner_join(fit %>% transmute(chamID, flux = .data[[paste0(g, "_flux")]], mdf = .data[[paste0(g, "_MDF")]]),
                      by = "chamID") %>%
      transmute(UniqueID = chamID, pre, start, flux, mdf) %>% flag() %>% mutate(gas = g)
  })) %>% mutate(dataset = "Agricultural soil")
}

# ── Write: replace only the datasets rescreened here, keep the rest of the saved list ──
cols <- c("dataset", "gas", "closure", "air_before_ppb", "window_start_ppb", "flux", "mdf")
saved <- if (file.exists(out_file)) read.csv(out_file, stringsAsFactors = FALSE) else
  setNames(data.frame(matrix(ncol = 7, nrow = 0)), cols)
screened <- c(stems = "Stems", soils = "Soils", coastal = "Coastal wetland", agsoil = "Agricultural soil")[names(res)]
new <- bind_rows(res)
if (nrow(new)) new <- new %>% transmute(dataset, gas, closure = UniqueID, air_before_ppb = round(pre, 1),
                                        window_start_ppb = round(start, 1), flux = signif(flux, 4), mdf = signif(mdf, 3))
inv <- bind_rows(saved %>% filter(!dataset %in% screened), new) %>% arrange(dataset, closure)
if (length(screened)) {
  write.csv(inv, out_file, row.names = FALSE)
  message("Start-above-ambient screen: rescreened ", paste(screened, collapse = ", "), "; ", nrow(inv), " invalid closures")
} else message("05: no raw records available; keeping ", out_file)
print(inv)
