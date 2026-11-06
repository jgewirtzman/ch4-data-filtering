# ============================================================================
# 24a_click_coastal_windows.R — Pick the coastal wetland closure windows by hand
#
# INTERACTIVE: run in RStudio or R (an on-screen graphics device is needed), from
# the repository root:  source("scripts/24a_click_coastal_windows.R")
#
# For each closure goFlux::click.peak2 shows the record from 3 min before the
# field-sheet start to 3 min after its end. Click the START, then the END of the
# steady part of the closure (after the headspace has mixed, before the chamber
# comes off or any step). CO2/CH4 closures are picked on CH4 with CO2 shown above;
# N2O closures on N2O with water vapour shown above.
#
# Resumable: every pick is appended to the output file at once, closures already
# in it are skipped, and to redo one, delete its row and run again.
#
# Loads:  data/local_full_sources/woogoompah/ (raw records and field sheet; local only)
# Writes: data/from_coastal_wetland_study/coastal_windows_clicked.csv
#         (UniqueID, analyzer, start.time_corr, end.time_corr, obs.length_corr, picked_at)
# ============================================================================

source("scripts/00_setup.R")
suppressMessages({library(readr); library(goFlux)})

src_dir  <- file.path(data_dir, "local_full_sources", "woogoompah")
out_file <- file.path(data_dir, "from_coastal_wetland_study", "coastal_windows_clicked.csv")
dir.create(dirname(out_file), showWarnings = FALSE)
SHOULDER <- 180   # s shown before the field start and after the field end
P <- function(x) as.POSIXct(x, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")

# --- field sheet ---------------------------------------------------------------
fs <- read_csv(file.path(src_dir, "field_sheet.csv"), col_names = FALSE, skip = 4,
               col_types = cols(.default = "c"))
block <- function(cols, analyzer) {
  x <- fs[, cols]; names(x)[c(1, 3, 4, 5)] <- c("UniqueID", "date", "start", "end")
  x %>% filter(!is.na(UniqueID)) %>%
    transmute(UniqueID, analyzer = analyzer,
              start.time = P(paste(substr(date, 1, 10), start)),
              end.time   = P(paste(substr(date, 1, 10), end)),
              end.time   = if_else(end.time < start.time, end.time + 86400, end.time),
              obs.length = as.numeric(end.time - start.time, units = "secs"))
}
aux <- bind_rows(block(10:23, "LI-7820"), block(27:40, "LI-7810"))

# --- raw records (corrected timestamps; non-zero diagnostics dropped) -----------
rd <- function(f) { r <- read_csv(file.path(src_dir, f), col_names = FALSE, skip = 7,
                                  col_types = cols(.default = "c"))
                    r[!is.na(r$X1) & r$X1 == "DATA" & r$X5 == "0", ] }
a <- rd("licor_n2o.csv"); b <- rd("licor_co2_ch4.csv")
raw <- list(
  "LI-7820" = data.frame(POSIX.time = P(a$X10), H2O_ppm = as.numeric(a$X11), N2Odry_ppb = as.numeric(a$X12)),
  "LI-7810" = data.frame(POSIX.time = P(b$X9), H2O_ppm = as.numeric(b$X10), CO2dry_ppm = as.numeric(b$X11),
                         CH4dry_ppb = as.numeric(b$X12)))
rm(a, b)

settings <- list(
  "LI-7810" = list(gastype = "CH4dry_ppb", gases = "CO2dry_ppm", plot.lim = c(1800, 3500)),
  "LI-7820" = list(gastype = "N2Odry_ppb", gases = "H2O_ppm",   plot.lim = c(280, 370)))

done <- if (file.exists(out_file)) read.csv(out_file)$UniqueID else character(0)
todo <- aux %>% filter(!UniqueID %in% done)
message(nrow(todo), " of ", nrow(aux), " closures still to pick (", length(done), " done)")

for (i in seq_len(nrow(todo))) {
  x <- todo[i, ]; st <- settings[[x$analyzer]]
  ow <- obs.win(inputfile = raw[[x$analyzer]], auxfile = as.data.frame(x[, c("UniqueID", "start.time", "obs.length")]),
                shoulder = SHOULDER)
  message("\n[", i, "/", nrow(todo), "] ", x$UniqueID, ": click START, then END")
  pk <- click.peak2(ow, gastype = st$gastype, gases = st$gases, plot.lim = st$plot.lim, sleep = 2)
  if (is.null(pk) || !nrow(pk)) { message("  no window picked; skipped"); next }
  row <- pk %>% distinct(UniqueID, start.time_corr, end.time_corr, obs.length_corr) %>%
    mutate(analyzer = x$analyzer, picked_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
           start.time_corr = format(start.time_corr, "%Y-%m-%d %H:%M:%S"),
           end.time_corr = format(end.time_corr, "%Y-%m-%d %H:%M:%S"))
  write.table(row, out_file, sep = ",", row.names = FALSE, append = file.exists(out_file),
              col.names = !file.exists(out_file))
}
message("\nPicked windows: ", out_file)
