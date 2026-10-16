# ============================================================================
# 14_figure_closures.R — Readings of the closures drawn in Fig. 2, Fig. 5 and the
#                        worked example (Sect. S1, Fig. S1), from the raw records
#
# Writes small, tracked CSVs to data/figure_closures/ so that 15_co2_validity_figure.R
# and 17_sigma_figures.R run without the raw archive (data/stem_raw/, ~650 MB):
#   fig5_co2_tracer.csv          Fig. 5: the four CO2-tracer example closures (LGR1)
#   fig2_sigma_day.csv           Fig. 2: every closure of the precise analyzer's example
#                                day (LI-7810, 2025-06-03)
#   worked_example_closures.csv  Sect. S1 / Fig. S1: the two noisier analyzer (LGR1)
#                                closures of the worked example
# Without data/stem_raw/ the script does nothing and the CSVs already in the
# repository are used. Numbers are written with as many digits as needed to read
# back bit-identical (checked below), so the figures do not change.
#
# Loads: data/stem_raw/LGR/LGR1/, outputs/02_allan_deviation.RData,
#        outputs/02b_manID.RData, outputs/02b_goflux_fluxes.RData
# Run after 02b.
# ============================================================================

source("scripts/00_setup.R")

local({
  if (!dir.exists(file.path(data_dir, "stem_raw", "LGR", "LGR1"))) {
    message("14: raw records (data/stem_raw/) not found; skipping. 15 and 17 use the ",
            "readings already in ", fig_closures_dir, "/")
    return(invisible(NULL))
  }
  dir.create(fig_closures_dir, recursive = TRUE, showWarnings = FALSE)

  # Shortest decimal text (15-17 significant digits) that parses back to the same double.
  # R's own number parser is not correctly rounded on every platform, so some doubles have
  # no decimal text that reads back exactly; LGR methane is therefore stored in ppm, the
  # unit the analyzer logs (6 significant digits), and converted to ppb by 15 and 17
  # exactly as when they read the raw files.
  num_txt <- function(x) {
    s <- sprintf("%.15g", x)
    for (d in 16:17) {
      bad <- which(!is.na(x) & suppressWarnings(as.numeric(s)) != x)
      if (length(bad) == 0) break
      s[bad] <- sprintf(paste0("%.", d, "g"), x[bad])
    }
    s[is.na(x)] <- NA
    s
  }
  # Write, then read back exactly as 15 and 17 do and require identical values
  write_checked <- function(df, file) {
    path <- file.path(fig_closures_dir, file)
    out <- as.data.frame(df)
    is_dbl <- vapply(out, is.double, logical(1))
    out[is_dbl] <- lapply(out[is_dbl], num_txt)
    write.csv(out, path, row.names = FALSE, na = "NA",
              quote = which(vapply(df, is.character, logical(1))))
    back <- read.csv(path, stringsAsFactors = FALSE,       # 15 and 17 read numbers as doubles
                     colClasses = setNames(rep("numeric", sum(is_dbl)), names(df)[is_dbl]))
    if (!identical(names(back), names(df)) || !all(mapply(identical, back, as.list(as.data.frame(df)))))
      stop("14: ", file, " does not read back identically")
    message("14: wrote ", path, " (", nrow(df), " readings, ",
            round(file.size(path) / 1024), " KB)")
  }

  # --------------------------------------------------------------------------
  # Fig. 5 (15_co2_validity_figure.R): the analyzer's "Time" column and dry CO2/CH4
  # inside each closure's computation window (same closures and windows as the
  # `measurements` table of 15; keep the two in step)
  # --------------------------------------------------------------------------
  fig5 <- data.frame(
    panel      = 1:4,
    UniqueID   = c("45509_414", "45489_153", "45580_439", "45610_414"),
    date_str   = c("2024-08-05", "2024-07-16", "2024-10-15", "2024-11-14"),
    comp_start = c("11:28:30", "12:04:13", "13:02:55", "11:31:00"),
    comp_end   = c("11:31:30", "12:07:13", "13:05:55", "11:34:00"),
    stringsAsFactors = FALSE)
  lgr1 <- file.path(data_dir, "stem_raw", "LGR", "LGR1")
  read_window <- function(date_str, comp_start, comp_end) {
    t_start <- as.POSIXct(paste(as.Date(date_str), comp_start), format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    t_end   <- as.POSIXct(paste(as.Date(date_str), comp_end),   format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    read_one <- function(fpath) {
      raw <- read.csv(fpath, skip = 1, stringsAsFactors = FALSE, check.names = FALSE)
      names(raw) <- trimws(names(raw))
      raw$Time_parsed <- as.POSIXct(trimws(raw$Time), format = "%m/%d/%Y %H:%M:%OS", tz = "UTC")
      raw[!is.na(raw$Time_parsed) & raw$Time_parsed >= t_start & raw$Time_parsed <= t_end, ]
    }
    dir_path <- file.path(lgr1, date_str)
    subsets <- lapply(list.files(dir_path, pattern = "^micro_.*_f[0-9]+\\.txt$", full.names = TRUE), read_one)
    subsets <- Filter(function(s) nrow(s) > 0, subsets)
    if (length(subsets) == 0) {           # zipped files only when the plain ones hold nothing
      for (z in list.files(dir_path, pattern = "^micro_.*_f[0-9]+\\.txt\\.zip$", full.names = TRUE)) {
        for (f in unzip(z, exdir = tempdir())) {
          s <- read_one(f); file.remove(f)
          if (nrow(s) > 0) subsets[[length(subsets) + 1]] <- s
        }
      }
    }
    if (length(subsets) == 0) stop("14: no readings for ", date_str, " ", comp_start, "-", comp_end)
    w <- do.call(rbind, subsets)
    w <- data.frame(time_s   = as.numeric(w$Time_parsed),
                    CO2d_ppm = as.numeric(trimws(w[["[CO2]d_ppm"]])),
                    CH4d_ppm = as.numeric(trimws(w[["[CH4]d_ppm"]])))   # 15 converts to ppb
    el <- w$time_s - as.numeric(t_start)
    w[el >= 0 & el <= 180, ]               # 15 shows a 180 s window
  }
  f5 <- do.call(rbind, lapply(seq_len(nrow(fig5)), function(i) {
    w <- read_window(fig5$date_str[i], fig5$comp_start[i], fig5$comp_end[i])
    data.frame(panel = fig5$panel[i], UniqueID = fig5$UniqueID[i], analyzer = "LGR1", w,
               stringsAsFactors = FALSE)
  }))
  write_checked(f5, "fig5_co2_tracer.csv")

  # --------------------------------------------------------------------------
  # Fig. 2 (17_sigma_figures.R): every closure window (flag == 1) of the precise
  # analyzer's example day, as extracted by 02b
  # --------------------------------------------------------------------------
  e <- new.env()
  load(file.path(output_dir, "02b_manID.RData"), envir = e)
  manID <- e$manID
  wk <- manID %>%
    filter(flag == 1, is.finite(CH4dry_ppb), analyzer == "LI-7810", date_str == "2025-06-03") %>%
    arrange(UniqueID, POSIX.time)
  write_checked(data.frame(UniqueID = wk$UniqueID, analyzer = wk$analyzer,
                           time_s = as.numeric(wk$POSIX.time), CH4dry_ppb = wk$CH4dry_ppb,
                           stringsAsFactors = FALSE),
                "fig2_sigma_day.csv")

  # --------------------------------------------------------------------------
  # Worked example (17_sigma_figures.R): the two closures chosen from 02b's windows,
  # re-read from the raw record over the field-log computation window used by 02
  # (02's segments and per-closure sigma are the source of every MDF in the paper)
  # --------------------------------------------------------------------------
  pick <- pick_worked_example(manID %>% filter(flag == 1) %>% arrange(UniqueID, POSIX.time))
  a <- new.env(); load(file.path(output_dir, "02_allan_deviation.RData"), envir = a)
  f <- new.env(); load(file.path(output_dir, "02b_goflux_fluxes.RData"), envir = f)
  stem_fl <- setNames(f$df$fl_UniqueID, f$df$stem_id)
  folders <- a$lgr_catalog$folder[a$lgr_catalog$machine == "LGR1" & a$lgr_catalog$date_str == pick$day]
  day_raw <- do.call(rbind, lapply(folders, load_lgr_day))
  we <- do.call(rbind, lapply(c(pick$low$UniqueID, pick$high$UniqueID), function(stem) {
    fl <- a$field_logs[a$field_logs$UniqueID == stem_fl[[stem]], ][1, ]
    s <- day_raw[day_raw$POSIX.time >= fl$comp_start_posix & day_raw$POSIX.time <= fl$comp_end_posix, ]
    ppm <- as.numeric(sprintf("%.6g", s$CH4dry_ppb / 1000))   # as logged (6 significant digits)
    stopifnot(identical(ppm * 1000, s$CH4dry_ppb))
    data.frame(UniqueID = stem, fl_UniqueID = stem_fl[[stem]], analyzer = "LGR1",
               time_s = as.numeric(s$POSIX.time), CH4d_ppm = ppm, stringsAsFactors = FALSE)
  }))
  write_checked(we, "worked_example_closures.csv")
})
