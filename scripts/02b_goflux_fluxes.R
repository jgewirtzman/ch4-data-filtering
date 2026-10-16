# ============================================================================
# 02b_goflux_fluxes.R — Re-fit every stem closure with the goFlux fork (v0.5.0.9002)
#
# Replaces the legacy fluxes in HF_2023-2025_tree_flux_corrected.csv (the source
# study's OLS on the field-log window for 2023-24; a 60-s best-CO2-R2 sub-window of the
# REMARK for 2025) with the convention of Gewirtzman et al. (in prep.) and the goFlux fork:
#
#   closure window  = field-log comp_start..comp_end (LGR, analyzer clock) or the
#                     REMARK run (LI-7810); no dead band, no sub-window search
#   fit             = goFlux::goFlux() LM + HM, goFlux::best.flux() with its default
#                     criteria (MAE, RMSE, AICc, SE, g-factor < 2, kappa <= kappa.max
#                     [Hüppi et al. 2018], MDF, nb.obs, intercept, p-value)
#   detection       = goFlux::flux.class(prec = that analyzer x field day's campaign
#                     sigma from 02 [closure median of second-difference precision],
#                     conf = 0.95 -> z = 1.96 benchmark multiplier, t = closure.time)
#   QC              = goFlux::qc.flags (c0, convex, min obs, ambient at sealing, clock,
#                     noisy) + goFlux::co2.tracer: flags only, nothing deleted
#
# Legacy values are kept side by side (CH4_flux_legacy, CH4_r2_legacy,
# CH4_SE_legacy, CO2_*_legacy); closures without an archived trace keep the
# legacy flux (flux_src = "legacy (no trace)") and are dropped in 03.
#
# Loads:  outputs/02_allan_deviation.RData (df, field_logs, closures_7810, mad_days)
# Saves:  outputs/02b_goflux_fluxes.RData   (df with goFlux fluxes; ch4_best, co2_best)
#         outputs/02b_manID.RData           (the fitted concentration windows)
#         outputs/02b_stem_fluxes_goflux.csv, outputs/02b_legacy_vs_goflux.csv,
#         outputs/02b_qc_review.csv, outputs/02b_qc_review_plots.pdf
# ============================================================================

source("scripts/00_setup.R")
suppressPackageStartupMessages(library(goFlux))   # fork >= 0.5.0.9002 (checked in 00_setup)

message("\n=== Loading 02 results ===")
load(file.path(output_dir, "02_allan_deviation.RData"))
tz_local <- "UTC"   # all stem clock times are labelled UTC (00_setup loaders, 02 field logs)
message("goFlux ", packageVersion("goFlux"), "; time zone ", tz_local)

# ============================================================================
# PART 1: ONE WINDOW PER CLOSURE WITH A TRACE
# ============================================================================

# settings (best.flux defaults incl. the Hüppi et al. 2018 g-factor / kappa criteria)
best_args <- list(g.limit = 2, k.ratio = 1, p.val = 0.05, warn.length = 60)
qc_set <- list(c0.mult = 1.5, min.obs = 60, ambient.sigma = 3, ambient.secs = 10, ambient.pre = 60,
               noisy.mult = 1.5, co2.tracer.p = 0.05)
# RESUME_02B=1 skips window extraction and fitting (~minutes) and reloads them from
# outputs/02b_checkpoint.RData + 02b_manID.RData written by a previous full run (README).
resume <- nzchar(Sys.getenv("RESUME_02B")) && file.exists(file.path(output_dir, "02b_checkpoint.RData")) &&
          file.exists(file.path(output_dir, "02b_manID.RData"))
if (resume) {
  message("RESUME_02B set: loading the checkpointed fits instead of re-fitting")
  load(file.path(output_dir, "02b_checkpoint.RData")); load(file.path(output_dir, "02b_manID.RData"))
}
message("\n=== Building closure windows ===")

# LGR (2023-24): the field-log entry each row was matched to in 02
win_lgr <- df %>%
  filter(year < 2025, !is.na(fl_UniqueID)) %>%
  select(stem_id, fl_UniqueID) %>%
  inner_join(field_logs %>%
               select(fl_UniqueID = UniqueID, machine, date_str,
                      start = comp_start_posix, end = comp_end_posix),
             by = "fl_UniqueID") %>%
  transmute(stem_id, analyzer = machine, date_str, start, end)

# LI-7810 (2025): the REMARK run matched in 02 (k-th visit of the tag that day)
win_7810 <- df %>%
  filter(year == 2025, !is.na(t_sec)) %>%
  select(stem_id, seg_key_7810) %>%
  inner_join(closures_7810 %>% select(seg_key_7810 = seg_key, date_str = DATE, start, end),
             by = "seg_key_7810") %>%
  transmute(stem_id, analyzer = "LI-7810", date_str, start, end)

windows <- bind_rows(win_lgr, win_7810) %>%
  mutate(obs.length = as.numeric(difftime(end, start, units = "secs")))
message("Windows: ", nrow(win_lgr), " LGR + ", nrow(win_7810), " LI-7810 = ",
        nrow(windows), " of ", nrow(df), " closures (",
        sum(is.na(df$t_sec)), " have no archived trace)")

# Chamber geometry per closure, in goFlux units
geom <- df %>%
  transmute(stem_id,
            Area  = surfarea * 1e4,                       # m2 -> cm2
            Vtot  = vol_system,                           # L (chamber + tubing)
            Pcham = if_else(is.na(bar), 101.3, bar / 10),  # hPa -> kPa
            Tcham = if_else(is.na(airt), 15, airt))       # deg C
n_geom_na <- sum(is.na(geom$Area) | is.na(geom$Vtot))
if (n_geom_na > 0) message("  ", n_geom_na, " closures lack chamber geometry and are skipped")

# ============================================================================
# PART 2: obs.win + windows.from.table, one analyzer x day at a time
# ============================================================================

if (!resume) {
message("\n=== Extracting windows from the raw records ===")

day_list <- windows %>% distinct(analyzer, date_str) %>% arrange(date_str, analyzer)
raw_7810_dir <- file.path(raw_dir, "7810", "raw-7810-tree_flux_data")
manID_all <- list()

for (i in seq_len(nrow(day_list))) {
  a <- day_list$analyzer[i]; d <- day_list$date_str[i]
  w <- windows %>% filter(analyzer == a, date_str == d) %>%
    inner_join(geom, by = "stem_id") %>% filter(!is.na(Area), !is.na(Vtot))
  if (nrow(w) == 0) next

  rec <- NULL
  if (a == "LI-7810") {
    f <- list.files(raw_7810_dir, pattern = paste0("^TG10-.*-", d, "T.*\\.data$"), full.names = TRUE)
    if (length(f)) rec <- tryCatch(load_li7810_day(f[1], tz = tz_local), error = function(e) NULL)
    if (!is.null(rec)) rec <- rec[, c("POSIX.time", "DATE", "CO2dry_ppm", "CH4dry_ppb", "H2O_ppm",
                                      "CO2_prec", "CH4_prec", "H2O_prec")]
  } else {
    folders <- lgr_catalog$folder[lgr_catalog$machine == a & lgr_catalog$date_str == d]
    if (length(folders) == 0) folders <- lgr_catalog$folder[lgr_catalog$date_str == d]
    for (folder in folders) { dd <- load_lgr_day(folder); if (!is.null(dd)) rec <- rbind(rec, dd) }
    if (!is.null(rec)) { rec <- rec[order(rec$POSIX.time), ]; rec <- rec[!duplicated(rec$POSIX.time), ] }
  }
  if (is.null(rec) || nrow(rec) < 10) { message("  ", d, " ", a, ": no record"); next }
  rec <- rec[!is.na(rec$CH4dry_ppb) & !is.na(rec$CO2dry_ppm), ]

  aux <- w %>% transmute(UniqueID = stem_id, start.time = start, obs.length,
                         Area, Vtot, Pcham, Tcham)
  ow <- tryCatch(suppressMessages(obs.win(inputfile = rec, auxfile = aux, shoulder = 60)),
                 error = function(e) { message("  ", d, " ", a, ": obs.win failed: ",
                                               conditionMessage(e)); NULL })
  if (is.null(ow) || length(ow) == 0) next
  man <- tryCatch(suppressMessages(suppressWarnings(
           windows.from.table(ow, w %>% transmute(UniqueID = stem_id, start, end),
                              tz = tz_local, warn.length = 60))),
         error = function(e) { message("  ", d, " ", a, ": windows.from.table failed: ",
                                       conditionMessage(e)); NULL })
  if (is.null(man) || nrow(man) == 0) next
  man$analyzer <- a; man$date_str <- d
  manID_all[[length(manID_all) + 1]] <- man
  if (i %% 10 == 0 || i == nrow(day_list))
    message("  ", i, "/", nrow(day_list), " analyzer-days; ",
            length(unique(do.call(rbind, manID_all)$UniqueID)), " windows so far")
  rm(rec, ow)
}

manID <- do.call(rbind, manID_all); rm(manID_all)
manID$grp <- paste(manID$analyzer, manID$date_str)
n_win <- length(unique(manID$UniqueID))
message("Windows extracted: ", n_win, " closures, ", nrow(manID), " rows (",
        sum(manID$flag == 1), " inside windows)")
win_len <- closure_seconds(manID)
message("  window length (s): median ", round(median(win_len$t_sec)), ", range ",
        round(min(win_len$t_sec)), "-", round(max(win_len$t_sec)),
        "; logging interval ", paste(unique(round(win_len$dt_s, 1)), collapse = "/"), " s")

# ============================================================================
# PART 3: goFlux LM + HM, best.flux
# ============================================================================

message("\n=== goFlux fits (LM + HM) ===")
ch4_flux <- suppressWarnings(goFlux(manID, "CH4dry_ppb", H2O_col = "H2O_ppm", warn.length = 60))
co2_flux <- suppressWarnings(goFlux(manID, "CO2dry_ppm", H2O_col = "H2O_ppm", warn.length = 60))
# best.flux defaults: MAE, RMSE, AICc, SE, g.factor (< 2), kappa (k.ratio 1;
# kappa.max after Hüppi et al. 2018), MDF, nb.obs, intercept, p-value
ch4_best <- do.call(best.flux, c(list(ch4_flux), best_args))
co2_best <- do.call(best.flux, c(list(co2_flux), best_args))
message("CH4 model chosen: "); print(table(ch4_best$model))
message("CO2 model chosen: "); print(table(co2_best$model))

# ============================================================================
# PART 4: detection class (campaign sigma of that analyzer x day) + QC flags
# ============================================================================

message("\n=== goFlux::flux.class (campaign sigma of each analyzer x day, z = 1.96) ===")
key <- manID %>% distinct(UniqueID, grp, analyzer, date_str)
ch4_best <- ch4_best %>% left_join(key, by = "UniqueID")
co2_best <- co2_best %>% left_join(key, by = "UniqueID")

# sigma per group: the closure-median second-difference precision of that analyzer x
# day from 02; fall back to the analyzer's median day (mad_inst_median) if a day is missing
sig_day <- mad_days %>% filter(!duplicated(paste(instrument, date_str))) %>%
  transmute(grp = paste(instrument, date_str), sig_CH4 = mad_sd_CH4, sig_CO2 = mad_sd_CO2)
sig_grp <- key %>% distinct(grp, analyzer) %>%
  left_join(sig_day, by = "grp") %>%
  left_join(mad_inst_median %>% transmute(analyzer = instrument, med_CH4 = mad_sd_CH4,
                                          med_CO2 = mad_sd_CO2), by = "analyzer") %>%
  mutate(sig_src = if_else(is.na(sig_CH4), "analyzer median day", "day"),
         sig_CH4 = coalesce(sig_CH4, med_CH4), sig_CO2 = coalesce(sig_CO2, med_CO2))
message("  sigma source: "); print(table(sig_grp$sig_src))

classify <- function(best, gas, sig_col) {
  pr <- data.frame(UniqueID = best$UniqueID, prec = sig_grp[[sig_col]][match(best$grp, sig_grp$grp)])
  res <- suppressWarnings(flux.class(best, dataframe = manID, gastype = gas, prec = pr, conf = 0.95))
  res$sigma_camp <- pr$prec
  res
}
ch4_best <- classify(ch4_best, "CH4dry_ppb", "sig_CH4")
co2_best <- classify(co2_best, "CO2dry_ppm", "sig_CO2")
message("CH4 detection class (campaign sigma, z = 1.96): "); print(table(ch4_best$det.class))

message("\n=== goFlux::qc.flags + co2.tracer ===")
ch4_best <- suppressWarnings(qc.flags(ch4_best, dataframe = manID, gastype = "CH4dry_ppb", by = "grp",
                                      c0.mult = qc_set$c0.mult, min.obs = qc_set$min.obs,
                                      ambient.sigma = qc_set$ambient.sigma, ambient.secs = qc_set$ambient.secs,
                                      ambient.pre = qc_set$ambient.pre, noisy.mult = qc_set$noisy.mult))
ch4_best$qc.co2.tracer <- !co2.tracer(co2_best, ch4_best, p.val = qc_set$co2.tracer.p)   # TRUE = CO2 not clearly rising
qc_cols <- c(qc.c0 = "C0 > 1.5x group median", qc.co2.tracer = "CO2 not clearly rising",
             qc.convex = "convex (accelerating) trace", qc.min.obs = "window < 60 obs",
             qc.ambient = "headspace at sealing off ambient", qc.clock = "clock mismatch",
             qc.noisy = "closure noisier than campaign")
ch4_best$qc.any.all <- Reduce(`|`, lapply(names(qc_cols), function(v) ch4_best[[v]] %in% TRUE))
ch4_best$qc.note <- apply(sapply(names(qc_cols), function(v) ifelse(ch4_best[[v]] %in% TRUE, qc_cols[[v]], NA)), 1,
                          function(r) paste(na.omit(r), collapse = "; "))
message("QC flags fired: "); print(sapply(names(qc_cols), function(v) sum(ch4_best[[v]] %in% TRUE)))
message("Closures with any QC flag: ", sum(ch4_best$qc.any.all), " / ", nrow(ch4_best))

# ============================================================================
# PART 5: MERGE INTO THE STEM TABLE, LEGACY SIDE BY SIDE
# ============================================================================

message("\n=== Merging into the stem table ===")
# checkpoint: the fitted tables and windows, so the merge can be re-run without re-fitting
save(ch4_best, co2_best, sig_grp, windows, file = file.path(output_dir, "02b_checkpoint.RData"))
save(manID, file = file.path(output_dir, "02b_manID.RData"))
}  # end if (!resume)
# (parameter names avoid dplyr data-mask clashes with the `model` column and stats::lm)
pick_model <- function(mod, hm_v, lm_v) ifelse(mod == "HM", as.numeric(hm_v), as.numeric(lm_v))
num_cols <- function(b) { b <- as.data.frame(b)
                          for (v in names(b)) if (is.data.frame(b[[v]]) || is.list(b[[v]])) b[[v]] <- unlist(b[[v]][[1]])
                          for (v in grep("^(LM|HM)\\.|^(C0|Ct|MDF|prec|flux\\.term|nb\\.obs|g\\.fact|best\\.flux)$", names(b), value = TRUE))
                            b[[v]] <- suppressWarnings(as.numeric(b[[v]])); b }
ch4_best <- num_cols(ch4_best); co2_best <- num_cols(co2_best)
new_ch4 <- ch4_best %>%
  transmute(stem_id = UniqueID,
            CH4_flux_goflux = best.flux, CH4_model = model,
            CH4_LM_flux = LM.flux, CH4_HM_flux = HM.flux,
            CH4_r2_goflux = pick_model(model, HM.r2, LM.r2),
            CH4_SE_goflux = pick_model(model, HM.SE, LM.SE),
            CH4_LM_r2 = LM.r2, CH4_LM_p = LM.p.val, CH4_g_fact = g.fact,
            CH4_quality_check = quality.check,
            CH4_MDF_det = det.MDF, CH4_det_class = det.class,
            nb_obs_goflux = nb.obs, flux_term_goflux = flux.term,
            qc_c0 = qc.c0, qc_co2_tracer = qc.co2.tracer, qc_convex = qc.convex,
            qc_min_obs = qc.min.obs, qc_ambient = qc.ambient, qc_clock = qc.clock,
            qc_noisy = qc.noisy, qc_any = qc.any.all, qc_note = qc.note)
new_co2 <- co2_best %>%
  transmute(stem_id = UniqueID,
            CO2_flux_goflux = best.flux, CO2_model = model,
            CO2_r2_goflux = pick_model(model, HM.r2, LM.r2),
            CO2_SE_goflux = pick_model(model, HM.SE, LM.SE),
            CO2_LM_p = LM.p.val, CO2_det_class = det.class)

df <- df %>%
  select(-any_of(c(names(new_ch4)[-1], names(new_co2)[-1]))) %>%
  left_join(new_ch4, by = "stem_id") %>%
  left_join(new_co2, by = "stem_id") %>%
  mutate(
    CH4_flux_legacy = CH4_flux_nmolpm2ps, CH4_r2_legacy = CH4_r2, CH4_SE_legacy = CH4_SE_corr,
    CO2_flux_legacy = CO2_flux_umolpm2ps, CO2_r2_legacy = CO2_r2, CO2_SE_legacy = CO2_SE_corr,
    flux_src = if_else(!is.na(CH4_flux_goflux), "goFlux", "legacy (no trace)"),
    CH4_flux_nmolpm2ps = coalesce(CH4_flux_goflux, CH4_flux_legacy),
    CH4_r2             = coalesce(CH4_r2_goflux, CH4_r2_legacy),
    CH4_SE_corr        = coalesce(CH4_SE_goflux, CH4_SE_legacy),
    CO2_flux_umolpm2ps = coalesce(CO2_flux_goflux, CO2_flux_legacy),
    CO2_r2             = coalesce(CO2_r2_goflux, CO2_r2_legacy),
    CO2_SE_corr        = coalesce(CO2_SE_goflux, CO2_SE_legacy),
    lm_p_CH4           = coalesce(CH4_LM_p, lm_p_CH4),
    CH4_snr_se         = abs(CH4_flux_nmolpm2ps) / CH4_SE_corr,
    CO2_snr_se         = abs(CO2_flux_umolpm2ps) / CO2_SE_corr
  )
message("Flux source: "); print(table(df$flux_src, df$year))

# --- Legacy vs goFlux comparison (closures with both) ---
cmp <- df %>% filter(flux_src == "goFlux") %>%
  mutate(analyzer = if_else(year == 2025, "B (LI-7810)", "A (LGR)"),
         d_best = CH4_flux_goflux - CH4_flux_legacy,
         d_LM   = CH4_LM_flux - CH4_flux_legacy,
         sign_change = sign(CH4_flux_goflux) != sign(CH4_flux_legacy),
         d_gt_mdf = abs(d_best) > CH4_MDF_det)
legacy_cmp <- bind_rows(cmp, cmp %>% mutate(analyzer = "all")) %>%
  group_by(analyzer) %>%
  summarise(n = n(),
            pct_HM = 100 * mean(CH4_model == "HM"),
            r_pearson  = cor(CH4_flux_goflux, CH4_flux_legacy),
            r_spearman = cor(CH4_flux_goflux, CH4_flux_legacy, method = "spearman"),
            median_legacy = median(CH4_flux_legacy), median_goflux = median(CH4_flux_goflux),
            mean_legacy = mean(CH4_flux_legacy), mean_goflux = mean(CH4_flux_goflux),
            median_ratio_LM = median(CH4_LM_flux / CH4_flux_legacy),
            median_ratio_best = median(CH4_flux_goflux / CH4_flux_legacy),
            median_abs_diff_LM = median(abs(d_LM)), median_abs_diff_best = median(abs(d_best)),
            pct_diff_gt_MDF = 100 * mean(d_gt_mdf, na.rm = TRUE),
            pct_sign_change = 100 * mean(sign_change),
            pct_neg_legacy = 100 * mean(CH4_flux_legacy < 0), pct_neg_goflux = 100 * mean(CH4_flux_goflux < 0),
            median_r2_legacy = median(CH4_r2_legacy, na.rm = TRUE), median_r2_goflux = median(CH4_r2_goflux),
            median_SE_legacy = median(CH4_SE_legacy, na.rm = TRUE), median_SE_goflux = median(CH4_SE_goflux),
            .groups = "drop")
message("\n=== Legacy vs goFlux (CH4, closures with a trace) ===")
print(as.data.frame(legacy_cmp), digits = 3)
write.csv(legacy_cmp, file.path(output_dir, "02b_legacy_vs_goflux.csv"), row.names = FALSE)

# ============================================================================
# PART 6: OUTPUTS
# ============================================================================

out_cols <- c("stem_id", "UniqueID", "fl_UniqueID", "date", "Tree", "PLOT", "SPECIES", "year",
              "instrument", "flux_src", "t_sec", "dt_s", "nb_obs_goflux",
              "CH4_flux_legacy", "CH4_r2_legacy", "CH4_SE_legacy",
              "CH4_flux_goflux", "CH4_model", "CH4_LM_flux", "CH4_HM_flux", "CH4_r2_goflux",
              "CH4_SE_goflux", "CH4_LM_p", "CH4_g_fact", "CH4_quality_check",
              "sigma_camp_CH4", "CH4_MDF_det", "CH4_det_class",
              "CO2_flux_legacy", "CO2_flux_goflux", "CO2_model", "CO2_r2_goflux", "CO2_det_class",
              "qc_c0", "qc_co2_tracer", "qc_convex", "qc_min_obs", "qc_ambient", "qc_clock",
              "qc_noisy", "qc_any", "qc_note")
write.csv(df[order(df$date, df$Tree), intersect(out_cols, names(df))],
          file.path(output_dir, "02b_stem_fluxes_goflux.csv"), row.names = FALSE)

# Closures the QC screens or goFlux's own checks caught: for review, not removal
review <- df %>%
  filter(flux_src == "goFlux", qc_any | grepl("g.factor|kappa|intercept", CH4_quality_check)) %>%
  arrange(date, Tree) %>%
  select(any_of(out_cols))
write.csv(review, file.path(output_dir, "02b_qc_review.csv"), row.names = FALSE)
message("\nQC review list: ", nrow(review), " closures -> outputs/02b_qc_review.csv")
message("  screens fired: ")
print(colSums(review[, c("qc_c0", "qc_co2_tracer", "qc_convex", "qc_min_obs",
                         "qc_ambient", "qc_clock", "qc_noisy")], na.rm = TRUE))

# one page per flagged closure (CH4 and CO2), for the review
pdf_ok <- tryCatch({
  rv_ids <- review$stem_id
  if (length(rv_ids) > 0) {
    p_ch4 <- flux.plot(ch4_best[ch4_best$UniqueID %in% rv_ids, ], manID[manID$UniqueID %in% rv_ids, ],
                       "CH4dry_ppb", shoulder = 60,
                       plot.display = c("MDF", "nb.obs", "flux.term"), quality.check = TRUE)
    p_co2 <- flux.plot(co2_best[co2_best$UniqueID %in% rv_ids, ], manID[manID$UniqueID %in% rv_ids, ],
                       "CO2dry_ppm", shoulder = 60, quality.check = TRUE)
    flux2pdf(c(p_ch4, p_co2), outfile = file.path(output_dir, "02b_qc_review_plots.pdf"))
  }
  TRUE
}, error = function(e) { message("  (review plots not written: ", conditionMessage(e), ")"); FALSE })

settings <- list(goFlux = as.character(packageVersion("goFlux")),
                 goFlux_source = "https://github.com/jgewirtzman/goFlux v0.5.0.9002 (doi:10.5281/zenodo.23256675)",
                 window = "field-log comp_start..comp_end (LGR) / REMARK run (LI-7810); no dead band",
                 shoulder_s = 60, best.flux = best_args,
                 detection = "flux.class: prec = analyzer x day campaign sigma (closure median of second-difference precision, 02), conf = 0.95 (z = 1.96 benchmark), t = closure.time",
                 qc = qc_set, date = as.character(Sys.Date()))
save(df, ch4_best, co2_best, legacy_cmp, sig_grp, settings, windows,
     file = file.path(output_dir, "02b_goflux_fluxes.RData"))
save(manID, file = file.path(output_dir, "02b_manID.RData"))
message("\nSaved: outputs/02b_goflux_fluxes.RData, 02b_manID.RData, 02b_stem_fluxes_goflux.csv,",
        " 02b_legacy_vs_goflux.csv, 02b_qc_review.csv")
message("=== 02b complete ===")
