# ============================================================================
# 02_allan_deviation.R — Per-closure precision for all flux datasets
#
# Computes per-closure instrument noise: the second-difference precision
# sigma_i (sc_*; sigma_closure() in 00_setup.R) and, for the x3t comparison, the
# Allan deviation at one logging interval (allan_sd_*); and the campaign
# precision (median of sigma_i over the closures of an analyzer x field day;
# the sigma that enters the reference MDF; stored as mad_sd_* for historical
# reasons, with the whole-record first-difference MAD kept as wr_sd_*) from:
#   1. HF canopy manID timeseries (LGR1/2/3)
#   2. YMF canopy manID timeseries
#   3. 2025 LI-7810 stem raw 1-Hz CSVs
#   4. 2023-24 LGR stem raw recordings (segmented by field log timestamps)
#
# Also runs rolling window scan on continuous canopy data to estimate
# global per-instrument precision (kept for comparison), records the logging
# interval of every record (LGR2 = 10 s, YMF = 5 s, others 1 s) and the closure
# length in seconds (t_sec; NOT nb.obs).
#
# Merges sigma_i, Allan deviation, campaign sigma and t_sec into the stem flux data frame (df).
# Saves all results to outputs/.
# ============================================================================

source("scripts/00_setup.R")

message("\n=== Loading data from 01_load_data.R ===")
load(file.path(output_dir, "01_loaded_data.RData"))

# Load YMF imp object (saved separately due to size)
ymf_imp_file <- file.path(output_dir, "01_ymf_imp.RData")
if (file.exists(ymf_imp_file)) {
  load(ymf_imp_file)
}

# ============================================================================
# PART 1: ALLAN DEVIATION — HF CANOPY (manID data)
# ============================================================================

message("\n=== Computing Allan deviation per measurement (HF canopy) ===")

allan_hf  <- compute_allan_per_uid(manID_all_hf)
allan_ymf <- compute_allan_per_uid(manID.YMF)

message("HF Allan deviation computed for ", nrow(allan_hf), " measurements")
message("  CO2 Allan SD range: ", round(min(allan_hf$allan_sd_CO2, na.rm = TRUE), 4),
        " - ", round(max(allan_hf$allan_sd_CO2, na.rm = TRUE), 4), " ppm")
message("  CH4 Allan SD range: ", round(min(allan_hf$allan_sd_CH4, na.rm = TRUE), 4),
        " - ", round(max(allan_hf$allan_sd_CH4, na.rm = TRUE), 4), " ppb")

message("YMF Allan deviation computed for ", nrow(allan_ymf), " measurements")
message("  CO2 Allan SD range: ", round(min(allan_ymf$allan_sd_CO2, na.rm = TRUE), 4),
        " - ", round(max(allan_ymf$allan_sd_CO2, na.rm = TRUE), 4), " ppm")
message("  CH4 Allan SD range: ", round(min(allan_ymf$allan_sd_CH4, na.rm = TRUE), 4),
        " - ", round(max(allan_ymf$allan_sd_CH4, na.rm = TRUE), 4), " ppb")

# ============================================================================
# PART 1b: CANOPY — CLOSURE LENGTH (s), LOGGING INTERVAL, CAMPAIGN SIGMA
# ============================================================================

message("\n=== Closure length and campaign precision (canopy) ===")

# Closure length in seconds per measurement (span of flagged rows + 1 interval)
tsec_hf  <- closure_seconds(manID_all_hf)
tsec_ymf <- closure_seconds(manID.YMF)

# Campaign sigma (reference): median over the record's closures of the per-closure
# second-difference precision (sigma_closure), from the flagged window rows. The whole-
# record MAD of first differences is kept as wr_sd_* for comparison.
# The YMF record is the LGR3 unit on a different day; treated as its own campaign.
canopy_records <- list(LGR1 = imp.LGR1, LGR2 = imp.LGR2, LGR3 = imp.LGR3,
                       YMF = ymf_imp)
# MAD is computed per run of constant logging interval (precision_mad_runs);
# a record with a mid-record settings change yields one row per interval
mad_results <- do.call(rbind, lapply(names(canopy_records), function(nm) {
  r <- canopy_records[[nm]]
  co2 <- precision_mad_runs(r$CO2dry_ppm, r$POSIX.time)
  ch4 <- precision_mad_runs(r$CH4dry_ppb, r$POSIX.time)
  data.frame(instrument = nm, dt_s = ch4$dt_s, mad_sd_CO2 = co2$mad_sigma[match(ch4$dt_s, co2$dt_s)],
             mad_sd_CH4 = ch4$mad_sigma, n_rows = ch4$n, stringsAsFactors = FALSE)
}))
if (any(duplicated(mad_results$instrument)))
  warning("a canopy record has more than one logging interval; downstream joins use the first")
mad_results <- mad_results[!duplicated(mad_results$instrument), ]
canopy_windows <- list(LGR1 = manID.LGR1, LGR2 = manID.LGR2, LGR3 = manID.LGR3, YMF = manID.YMF)
mad_results <- mad_results %>%
  rename(wr_sd_CO2 = mad_sd_CO2, wr_sd_CH4 = mad_sd_CH4) %>%
  rowwise() %>%
  mutate(w = list(canopy_windows[[instrument]] %>% filter(flag == 1) %>% arrange(UniqueID, POSIX.time)),
         mad_sd_CH4 = sigma_campaign(w$CH4dry_ppb, w$UniqueID, w$POSIX.time)[["sigma"]],
         mad_sd_CO2 = sigma_campaign(w$CO2dry_ppm, w$UniqueID, w$POSIX.time)[["sigma"]],
         n_closures = sigma_campaign(w$CH4dry_ppb, w$UniqueID, w$POSIX.time)[["n_closures"]]) %>%
  select(-w) %>% ungroup() %>% as.data.frame()
# per-closure sigma for every canopy closure (figures)
sigma_closures_canopy <- bind_rows(lapply(names(canopy_windows), function(nm) {
  w <- canopy_windows[[nm]] %>% filter(flag == 1) %>% arrange(UniqueID, POSIX.time)
  w %>% group_by(UniqueID) %>% summarise(sc_CH4 = sigma_closure(CH4dry_ppb, POSIX.time), sc_CO2 = sigma_closure(CO2dry_ppm, POSIX.time),
                                          .groups = "drop") %>% mutate(instrument = nm) }))
message("Campaign sigma (closure median, 2nd differences) vs whole-record MAD, per canopy record:")
print(mad_results)

# ============================================================================
# PART 2: LOAD FIELD LOGS (for LGR measurement timestamps)
# ============================================================================

message("\n=== Loading field logs for LGR measurement timestamps ===")

# Closure times of the 2023-2024 field logs: the four logs (Summer 2023, Sept 2023 - Sept 2024,
# Summer 2024 and the Sept-Dec 2024 sheets) reduced to the columns used here, as text
# (data/from_hf_stem_study/README.md): closure ID, date, tree tag, computation window
# (QC-updated times where the Sept 2023 - Sept 2024 log has them), analyzer start time,
# analyzer unit as logged, and the log each row came from
field_logs <- read.csv(file.path(stem_dir, "field_log_closure_times.csv"), colClasses = "character")
message("Combined field logs: ", nrow(field_logs), " entries")

# Standardize machine names
field_logs$machine <- gsub("LGR #", "LGR", field_logs$machine)
field_logs$machine <- gsub("\\s+", "", field_logs$machine)
field_logs$machine[is.na(field_logs$machine) |
                    field_logs$machine == ""] <- "LGR1"  # default

# Parse dates — handle multiple formats:
#   "9/21/2023" or "01/05/2024" (%m/%d/%Y)
#   "2/27/24" or "3/6/24"       (%m/%d/%y — 2-digit year!)
#   "6-29-2023"                  (%m-%d-%Y)
#   "24-02-27"                   (%y-%m-%d)
field_logs$date_parsed <- as.Date(field_logs$date_raw,
                                   format = "%m/%d/%Y")
# Fix 2-digit years that %Y parsed as year 24 AD instead of 2024
bad_year <- !is.na(field_logs$date_parsed) &
            as.integer(format(field_logs$date_parsed, "%Y")) < 100
if (any(bad_year)) {
  field_logs$date_parsed[bad_year] <- as.Date(
    field_logs$date_raw[bad_year], format = "%m/%d/%y")
}
# Try alternate format for dates like "6-29-2023"
na_dates <- is.na(field_logs$date_parsed)
if (any(na_dates)) {
  field_logs$date_parsed[na_dates] <- as.Date(
    field_logs$date_raw[na_dates], format = "%m-%d-%Y")
}
# Handle 2-digit years like "24-02-27"
na_dates <- is.na(field_logs$date_parsed)
if (any(na_dates)) {
  field_logs$date_parsed[na_dates] <- as.Date(
    field_logs$date_raw[na_dates], format = "%y-%m-%d")
}
# Handle ISO format "2024-09-17" (Sept-Dec 2024 sheets: readxl dates converted to text)
na_dates <- is.na(field_logs$date_parsed)
if (any(na_dates)) {
  field_logs$date_parsed[na_dates] <- as.Date(
    field_logs$date_raw[na_dates], format = "%Y-%m-%d")
}

message("Field log date range: ",
        min(field_logs$date_parsed, na.rm = TRUE), " to ",
        max(field_logs$date_parsed, na.rm = TRUE))
message("Dates with NAs: ", sum(is.na(field_logs$date_parsed)))

# Build comp_start/comp_end as POSIXct (date + time)
field_logs <- field_logs %>%
  filter(!is.na(date_parsed)) %>%
  mutate(
    # analyzer clock times, labelled UTC like the raw records (00_setup loaders)
    comp_start_posix = as.POSIXct(paste(date_parsed, comp_start),
                                    format = "%Y-%m-%d %H:%M:%S", tz = "UTC"),
    comp_end_posix   = as.POSIXct(paste(date_parsed, comp_end),
                                    format = "%Y-%m-%d %H:%M:%S", tz = "UTC"),
    t_sec_fl         = as.numeric(difftime(comp_end_posix,
                                            comp_start_posix,
                                            units = "secs")),
    date_str         = format(date_parsed, "%Y-%m-%d")
  )

# Remove entries with unparseable times or unreasonable durations
field_logs <- field_logs %>%
  filter(!is.na(comp_start_posix), !is.na(comp_end_posix),
         t_sec_fl > 30, t_sec_fl < 1800)

message("After cleaning: ", nrow(field_logs), " field log entries with valid timestamps")
message("  Measurement duration range: ",
        round(min(field_logs$t_sec_fl)), " - ",
        round(max(field_logs$t_sec_fl)), " sec")

# De-duplicate: prefer updated2 over summer2023 (has QC corrections), and
# prefer summer2024 over updated2 (more complete for those dates)
field_logs <- field_logs %>%
  mutate(priority = case_when(
    source == "xlsx_late2024" ~ 4,
    source == "summer2024" ~ 3,
    source == "updated2"   ~ 2,
    source == "summer2023" ~ 1
  )) %>%
  group_by(UniqueID) %>%
  slice_max(priority, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(-priority)

message("After de-duplication: ", nrow(field_logs), " unique field log entries")

# ============================================================================
# PART 3: ALLAN DEVIATION — 2025 LI-7810 STEM DATA
# ============================================================================

message("\n=== Allan deviation: 2025 LI-7810 ===")

# Closures are the REMARK runs of the instrument's own continuous record (the
# operator types the tree tag at closure, clears it at opening): the whole
# closure, not the 60-s best-CO2-R2 sub-window the legacy tables archived. The
# whole-day first-difference MAD computed here is replaced below by the closure
# median of sigma_i (the campaign sigma), as for the LGR days.
raw_dir_7810   <- file.path(raw_dir, "7810", "raw-7810-tree_flux_data")
raw_files_7810 <- list.files(raw_dir_7810, pattern = "\\.data$", full.names = TRUE)
message("Found ", length(raw_files_7810), " raw LI-7810 daily .data files")

allan_7810    <- list()
mad_7810_days <- list()
for (f in raw_files_7810) {
  rec <- tryCatch(load_li7810_day(f), error = function(e) {
    message("  could not import ", basename(f), ": ", conditionMessage(e)); NULL })
  if (is.null(rec) || nrow(rec) < 10) next
  for (d in sort(unique(rec$DATE))) {
    day <- rec[rec$DATE == d, ]
    co2_runs <- precision_mad_runs(day$CO2dry_ppm, day$POSIX.time)
    ch4_runs <- precision_mad_runs(day$CH4dry_ppb, day$POSIX.time)
    mad_7810_days[[length(mad_7810_days) + 1]] <- data.frame(
      instrument = "LI-7810", date_str = d, dt_s = ch4_runs$dt_s,
      mad_sd_CO2 = co2_runs$mad_sigma[match(ch4_runs$dt_s, co2_runs$dt_s)],
      mad_sd_CH4 = ch4_runs$mad_sigma, n_rows = ch4_runs$n, stringsAsFactors = FALSE)
    segs <- segment_remarks(day)
    if (is.null(segs)) next
    for (j in seq_len(nrow(segs))) {
      sub <- day[segs$i0[j]:segs$i1[j], ]
      if (nrow(sub) < 5) next
      dt_seg <- logging_interval(sub$POSIX.time)
      sf <- slope_fit_stats(sub$CH4dry_ppb,
                            as.numeric(sub$POSIX.time - sub$POSIX.time[1], units = "secs"))
      allan_7810[[length(allan_7810) + 1]] <- data.frame(
        seg_key      = segs$seg_key[j],
        REMARK       = segs$REMARK[j],
        DATE         = d,
        k            = segs$k[j],
        start        = segs$start[j],
        end          = segs$end[j],
        allan_sd_CO2 = allan_sd(sub$CO2dry_ppm),
        allan_sd_CH4 = allan_sd(sub$CH4dry_ppb),
        sc_CO2       = sigma_closure(sub$CO2dry_ppm, sub$POSIX.time),
        sc_CH4       = sigma_closure(sub$CH4dry_ppb, sub$POSIX.time),
        n_pts        = nrow(sub),
        t_sec        = as.numeric(diff(range(sub$POSIX.time)), units = "secs") + dt_seg,
        dt_s         = dt_seg,
        instrument   = "LI-7810",
        lm_p_CH4     = sf[["lm_p"]],
        lm_ac1_CH4   = sf[["ac1"]],
        stringsAsFactors = FALSE)
    }
    message("  ", d, ": ", nrow(segs), " REMARK closures, ", nrow(day), " rows")
  }
  rm(rec)
}
allan_df_7810 <- do.call(rbind, allan_7810)
mad_7810_days <- do.call(rbind, mad_7810_days)
if (!is.null(mad_7810_days)) {
  message("LI-7810 campaign sigma: ", nrow(mad_7810_days), " days; CH4 MAD median ",
          round(median(mad_7810_days$mad_sd_CH4), 4), " ppb")
}

if (is.null(allan_df_7810) || nrow(allan_df_7810) == 0) {
  warning("No LI-7810 Allan results. Continuing with LGR only.")
  allan_df_7810 <- data.frame(
    seg_key = character(), REMARK = character(), DATE = character(), k = integer(),
    start = as.POSIXct(character()), end = as.POSIXct(character()),
    allan_sd_CO2 = numeric(), allan_sd_CH4 = numeric(), sc_CO2 = numeric(), sc_CH4 = numeric(),
    n_pts = integer(), t_sec = numeric(), dt_s = numeric(), instrument = character(),
    lm_p_CH4 = numeric(), lm_ac1_CH4 = numeric(),
    stringsAsFactors = FALSE)
} else {
  message("LI-7810 closures (REMARK runs): ", nrow(allan_df_7810), " on ",
          length(unique(allan_df_7810$DATE)), " days; repeat visits of a tag on a day: ",
          sum(allan_df_7810$k > 1))
  message("  CH4 Allan SD: median ", round(median(allan_df_7810$allan_sd_CH4, na.rm = TRUE), 4),
          " ppb, range ", round(min(allan_df_7810$allan_sd_CH4, na.rm = TRUE), 4),
          " - ", round(max(allan_df_7810$allan_sd_CH4, na.rm = TRUE), 4))
  message("  closure length (s): median ", round(median(allan_df_7810$t_sec)),
          ", range ", round(min(allan_df_7810$t_sec)), " - ", round(max(allan_df_7810$t_sec)))
}
# closure windows for 02b (instrument time)
closures_7810 <- allan_df_7810[, c("seg_key", "REMARK", "DATE", "k", "start", "end", "n_pts", "t_sec")]

# ============================================================================
# PART 4: ALLAN DEVIATION — 2023-24 LGR/UGGA STEM DATA
# ============================================================================

message("\n=== Allan deviation: 2023-24 LGR/UGGA ===")

lgr_base <- file.path(raw_dir, "LGR")
message("LGR raw data directory: ", lgr_base)
message("  Exists: ", dir.exists(lgr_base))

# Catalog all available LGR date folders
lgr_catalog <- data.frame(machine = character(), date_str = character(),
                           folder = character(), stringsAsFactors = FALSE)
for (lgr_name in c("LGR1", "LGR2", "LGR3")) {
  lgr_dir <- file.path(lgr_base, lgr_name)
  if (!dir.exists(lgr_dir)) next
  # Folder names like "2023-06-29", "2023-06-06 (1)", "BGS2025-01-09"; some
  # archive several days as sub-folders ("2024-12-03/2024-12-02"), so every
  # directory (any depth) whose name carries a date and holds _f data files is
  # catalogued under that date.
  all_dirs <- list.dirs(lgr_dir, recursive = TRUE, full.names = TRUE)
  for (folder in all_dirs) {
    date_part <- regmatches(basename(folder), regexpr("\\d{4}-\\d{2}-\\d{2}", basename(folder)))
    if (length(date_part) == 0) next
    has_data <- any(grepl("_f\\d+\\.txt(\\.zip)?$", list.files(folder)))
    if (!has_data) next
    lgr_catalog <- rbind(lgr_catalog, data.frame(
      machine = lgr_name, date_str = date_part, folder = folder,
      stringsAsFactors = FALSE
    ))
  }
}
message("LGR catalog: ", nrow(lgr_catalog), " date-folders across ",
        length(unique(lgr_catalog$machine)), " instruments")

# parse_lgr_file() / load_lgr_day() are defined in 00_setup.R (shared with 02b)

# Process each field log entry: find LGR data, segment, compute Allan deviation
message("\nProcessing LGR measurements...")

# Group field log entries by date + machine for efficiency
# (load LGR file once per date, segment all measurements)
field_logs_lgr <- field_logs %>%
  filter(date_parsed < as.Date("2025-01-01"))

# Track unique dates to load
unique_day_machine <- field_logs_lgr %>%
  distinct(date_str, machine) %>%
  arrange(date_str)

message("  Unique LGR date-machine combos: ", nrow(unique_day_machine))

allan_lgr_all <- list()
mad_lgr_days  <- list()   # campaign sigma per analyzer x field day
n_matched <- 0
n_missing <- 0

for (i in seq_len(nrow(unique_day_machine))) {
  d <- unique_day_machine$date_str[i]
  m <- unique_day_machine$machine[i]

  # Find matching folder in catalog
  folders <- lgr_catalog$folder[lgr_catalog$machine == m &
                                  lgr_catalog$date_str == d]

  if (length(folders) == 0) {
    # Try other LGR instruments as fallback
    folders <- lgr_catalog$folder[lgr_catalog$date_str == d]
  }

  if (length(folders) == 0) {
    n_missing <- n_missing + 1
    next
  }

  # Load LGR data for this day (concatenate all folders if multiple)
  day_data <- NULL
  for (folder in folders) {
    dd <- load_lgr_day(folder)
    if (!is.null(dd)) {
      day_data <- rbind(day_data, dd)
    }
  }

  if (is.null(day_data) || nrow(day_data) < 10) {
    n_missing <- n_missing + 1
    next
  }
  day_data <- day_data[order(day_data$POSIX.time), ]
  day_data <- day_data[!duplicated(day_data$POSIX.time), ]

  # Logging interval and campaign sigma for this analyzer x field day
  # (LGR2 logs every 10 s; its sigma is at that interval)
  co2_runs <- precision_mad_runs(day_data$CO2dry_ppm, day_data$POSIX.time)
  ch4_runs <- precision_mad_runs(day_data$CH4dry_ppb, day_data$POSIX.time)
  mad_lgr_days[[length(mad_lgr_days) + 1]] <- data.frame(
    instrument = m, date_str = d, dt_s = ch4_runs$dt_s,
    mad_sd_CO2 = co2_runs$mad_sigma[match(ch4_runs$dt_s, co2_runs$dt_s)],
    mad_sd_CH4 = ch4_runs$mad_sigma,          # ppb (loader converts ppm -> ppb)
    n_rows = ch4_runs$n, stringsAsFactors = FALSE)

  # Get all field log entries for this date+machine
  entries <- field_logs_lgr %>%
    filter(date_str == d, machine == m)

  # Segment each measurement and compute Allan deviation
  for (j in seq_len(nrow(entries))) {
    start_t <- entries$comp_start_posix[j]
    end_t   <- entries$comp_end_posix[j]

    segment <- day_data[day_data$POSIX.time >= start_t &
                          day_data$POSIX.time <= end_t, ]

    if (nrow(segment) < 5) next

    n_matched <- n_matched + 1
    dt_seg <- logging_interval(segment$POSIX.time)   # this closure's own spacing
    sf <- slope_fit_stats(segment$CH4dry_ppb,
                          as.numeric(segment$POSIX.time - segment$POSIX.time[1], units = "secs"))
    allan_lgr_all[[n_matched]] <- data.frame(
      UniqueID     = entries$UniqueID[j],
      date_str     = d,
      allan_sd_CO2 = allan_sd(segment$CO2dry_ppm),   # ppm (same as LI-7810)
      allan_sd_CH4 = allan_sd(segment$CH4dry_ppb),   # ppb (same as LI-7810)
      sc_CO2       = sigma_closure(segment$CO2dry_ppm, segment$POSIX.time),
      sc_CH4       = sigma_closure(segment$CH4dry_ppb, segment$POSIX.time),
      n_pts        = nrow(segment),
      # closure length in seconds: span + one logging interval (= nrow only at 1 Hz)
      t_sec        = as.numeric(diff(range(segment$POSIX.time)), units = "secs") + dt_seg,
      dt_s         = dt_seg,
      instrument   = m,
      lm_p_CH4     = sf[["lm_p"]],
      lm_ac1_CH4   = sf[["ac1"]],
      stringsAsFactors = FALSE
    )
  }

  # Progress message every 10 dates
  if (i %% 10 == 0 || i == nrow(unique_day_machine)) {
    message("  Processed ", i, "/", nrow(unique_day_machine),
            " dates (", n_matched, " measurements matched so far)")
  }
}

allan_df_lgr <- do.call(rbind, allan_lgr_all)
mad_lgr_days <- do.call(rbind, mad_lgr_days)
if (!is.null(mad_lgr_days)) {
  message("\nLGR campaign sigma (MAD, per analyzer x day):")
  print(mad_lgr_days %>% group_by(instrument) %>%
          summarise(days = n(), dt_s = median(dt_s),
                    mad_CH4_median = median(mad_sd_CH4), mad_CO2_median = median(mad_sd_CO2),
                    .groups = "drop"))
}

if (is.null(allan_df_lgr) || nrow(allan_df_lgr) == 0) {
  warning("No LGR Allan results computed.")
  allan_df_lgr <- data.frame(
    UniqueID = character(), date_str = character(),
    allan_sd_CO2 = numeric(), allan_sd_CH4 = numeric(), sc_CO2 = numeric(), sc_CH4 = numeric(),
    n_pts = integer(), t_sec = numeric(), dt_s = numeric(), instrument = character(),
    lm_p_CH4 = numeric(), lm_ac1_CH4 = numeric(),
    stringsAsFactors = FALSE
  )
} else {
  message("\nLGR Allan deviation: ", nrow(allan_df_lgr), " measurements")
  message("  CH4 Allan SD: median ", round(median(allan_df_lgr$allan_sd_CH4, na.rm = TRUE), 4),
          " ppb, range ", round(min(allan_df_lgr$allan_sd_CH4, na.rm = TRUE), 4),
          " - ", round(max(allan_df_lgr$allan_sd_CH4, na.rm = TRUE), 4))
  message("  Missing date-folders: ", n_missing)
}

# ============================================================================
# PART 5: MERGE ALLAN RESULTS INTO STEM FLUX DATA
# ============================================================================

message("\n=== Merging Allan deviation into corrected flux data ===")

# Stable row id for the stem table (UniqueID is NA for 2025 and garbled for
# some 2023 rows); 02b keys the goFlux results on it.
df$stem_id <- sprintf("s%04d", seq_len(nrow(df)))

# --- Match 2025 LI-7810 by Tree_Date key: the k-th visit of that tag that day
# (6 trees were closed twice on a day; rows are ordered by their timestamp) ---
df <- df %>%
  group_by(Tree, date) %>%
  mutate(k_day = rank(coalesce(as.numeric(as.POSIXct(datetime_posx,
                                format = "%Y-%m-%d %H:%M:%S", tz = "UTC")), 0),
                      ties.method = "first")) %>%
  ungroup()
df$seg_key_7810 <- paste(df$Tree, df$date, df$k_day, sep = "_")

allan_7810_unique <- allan_df_7810 %>%
  select(seg_key, allan_sd_CO2, allan_sd_CH4, sc_CO2, sc_CH4, n_pts, t_sec, dt_s, instrument,
         lm_p_CH4, lm_ac1_CH4)

df <- merge(df, allan_7810_unique,
            by.x = "seg_key_7810", by.y = "seg_key",
            all.x = TRUE, suffixes = c("", "_7810"))
# Fallback: a row whose k-th REMARK run is missing (a run too short to keep, or
# a label typo) takes the one unused run of that tag on that day, if unique
fb_cols <- c("allan_sd_CO2", "allan_sd_CH4", "sc_CO2", "sc_CH4", "n_pts", "t_sec", "dt_s",
             "instrument", "lm_p_CH4", "lm_ac1_CH4")
n_fb <- 0
for (idx in which(df$year == 2025 & is.na(df$t_sec))) {
  cand <- allan_df_7810[allan_df_7810$REMARK == as.character(df$Tree[idx]) &
                          allan_df_7810$DATE == df$date[idx] &
                          !(allan_df_7810$seg_key %in% df$seg_key_7810), ]
  if (nrow(cand) != 1) next
  df[idx, fb_cols] <- cand[1, fb_cols]
  df$seg_key_7810[idx] <- cand$seg_key[1]
  n_fb <- n_fb + 1
}
message("LI-7810 rows matched through the unused-run fallback: ", n_fb)

# --- Match 2023-24 LGR by UniqueID ---
if (nrow(allan_df_lgr) > 0) {
  allan_lgr_unique <- allan_df_lgr %>%
    group_by(UniqueID) %>%
    slice(1) %>%
    ungroup() %>%
    select(UniqueID, allan_sd_CO2, allan_sd_CH4, sc_CO2, sc_CH4, n_pts, t_sec, dt_s, instrument,
           lm_p_CH4, lm_ac1_CH4)

  df <- merge(df, allan_lgr_unique,
              by = "UniqueID", all.x = TRUE,
              suffixes = c("", "_lgr"))

  # Combine: prefer existing (7810) values, fill in LGR where missing
  df <- df %>%
    mutate(
      fl_UniqueID  = if_else(!is.na(allan_sd_CH4_lgr), UniqueID, NA_character_),
      allan_sd_CO2 = if_else(is.na(allan_sd_CO2), allan_sd_CO2_lgr, allan_sd_CO2),
      allan_sd_CH4 = if_else(is.na(allan_sd_CH4), allan_sd_CH4_lgr, allan_sd_CH4),
      sc_CO2       = if_else(is.na(sc_CO2), sc_CO2_lgr, sc_CO2),
      sc_CH4       = if_else(is.na(sc_CH4), sc_CH4_lgr, sc_CH4),
      n_pts        = if_else(is.na(n_pts), n_pts_lgr, n_pts),
      t_sec        = if_else(is.na(t_sec), t_sec_lgr, t_sec),
      dt_s         = if_else(is.na(dt_s), dt_s_lgr, dt_s),
      lm_p_CH4     = if_else(is.na(lm_p_CH4), lm_p_CH4_lgr, lm_p_CH4),
      lm_ac1_CH4   = if_else(is.na(lm_ac1_CH4), lm_ac1_CH4_lgr, lm_ac1_CH4),
      instrument   = if_else(is.na(instrument) | instrument == "",
                              if_else(!is.na(instrument_lgr), instrument_lgr, NA_character_),
                              instrument)
    ) %>%
    select(-ends_with("_lgr"))
}

# --- For 2023 garbled UniqueIDs, try matching by datetime ---
# The corrected flux file datetime_posx = "Real start" from field logs
# Match 2023-24 unmatched rows by (date + tree + datetime proximity)
unmatched_idx <- which(is.na(df$allan_sd_CH4) & df$year < 2025)
message("Unmatched 2023-24 rows before datetime fallback: ", length(unmatched_idx))

if (length(unmatched_idx) > 0 && nrow(field_logs_lgr) > 0 && nrow(allan_df_lgr) > 0) {
  # Build lookup from field logs: real_start datetime + tree -> UniqueID
  fl_lookup <- field_logs %>%
    filter(!is.na(date_parsed)) %>%
    mutate(
      real_start_posix = as.POSIXct(paste(date_parsed, real_start),
                                     format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    ) %>%
    filter(!is.na(real_start_posix))

  # For each unmatched row, try to find a field log entry with matching
  # datetime_posx (within 2 minutes tolerance)
  for (idx in unmatched_idx) {
    dt_str <- df$datetime_posx[idx]
    if (is.na(dt_str) || dt_str == "") next

    # Parse datetime from corrected flux file
    dt <- as.POSIXct(dt_str, format = "%Y-%m-%dT%H:%M:%OSZ", tz = "UTC")
    if (is.na(dt)) dt <- as.POSIXct(dt_str, format = "%Y-%m-%dT%H:%M:%S", tz = "UTC")
    if (is.na(dt)) next

    # Compare with field log real_start times on same date
    same_date <- fl_lookup$date_str == df$date[idx]
    if (!any(same_date, na.rm = TRUE)) next

    fl_sub <- fl_lookup[same_date & !is.na(same_date), ]
    if (nrow(fl_sub) == 0) next

    # Find closest match within 2 minutes
    time_diffs <- abs(as.numeric(difftime(fl_sub$real_start_posix, dt, units = "secs")))
    best <- which.min(time_diffs)
    if (time_diffs[best] > 120) next  # no match within 2 min

    # Get UniqueID from field log and look up Allan result
    matched_uid <- fl_sub$UniqueID[best]
    allan_row <- allan_df_lgr[allan_df_lgr$UniqueID == matched_uid, ]
    if (nrow(allan_row) == 0) next

    df$allan_sd_CO2[idx] <- allan_row$allan_sd_CO2[1]
    df$allan_sd_CH4[idx] <- allan_row$allan_sd_CH4[1]
    df$sc_CO2[idx]       <- allan_row$sc_CO2[1]
    df$sc_CH4[idx]       <- allan_row$sc_CH4[1]
    df$n_pts[idx]        <- allan_row$n_pts[1]
    df$t_sec[idx]        <- allan_row$t_sec[1]
    df$dt_s[idx]         <- allan_row$dt_s[1]
    df$lm_p_CH4[idx]     <- allan_row$lm_p_CH4[1]
    df$lm_ac1_CH4[idx]   <- allan_row$lm_ac1_CH4[1]
    df$instrument[idx]   <- allan_row$instrument[1]
    df$fl_UniqueID[idx]  <- matched_uid
  }
}

unmatched_final <- sum(is.na(df$allan_sd_CH4) & df$year < 2025)
message("Unmatched 2023-24 after datetime fallback: ", unmatched_final)

# Summary
message("\n=== Allan deviation coverage (stem) ===")
message("  2025 LI-7810: ", sum(!is.na(df$allan_sd_CH4) & df$year == 2025),
        " / ", sum(df$year == 2025))
message("  2023-24 LGR:  ", sum(!is.na(df$allan_sd_CH4) & df$year < 2025),
        " / ", sum(df$year < 2025))
message("  Total:        ", sum(!is.na(df$allan_sd_CH4)), " / ", n_total)

# --- Campaign sigma per stem measurement ---
# Reference campaign sigma = median of the analyzer x day's per-closure second-difference
# sigma (swap_sigma below). Rows whose raw record was not found get the analyzer-wide
# median over days (machine from the field log when known, else the pooled LGR median).
# The whole-record first-difference MAD computed above is kept as wr_sd_* (comparison).
closure_median_days <- bind_rows(
  allan_df_7810 %>% transmute(instrument, date_str = DATE, sc_CH4, sc_CO2),
  allan_df_lgr  %>% transmute(instrument, date_str, sc_CH4, sc_CO2)) %>%
  group_by(instrument, date_str) %>%
  summarise(cm_CH4 = median(sc_CH4, na.rm = TRUE), cm_CO2 = median(sc_CO2, na.rm = TRUE),
            n_closures = sum(is.finite(sc_CH4)), .groups = "drop")
swap_sigma <- function(d) d %>%
  rename(wr_sd_CO2 = mad_sd_CO2, wr_sd_CH4 = mad_sd_CH4) %>%
  left_join(closure_median_days, by = c("instrument", "date_str")) %>%
  mutate(mad_sd_CH4 = cm_CH4, mad_sd_CO2 = cm_CO2) %>% select(-cm_CH4, -cm_CO2)
mad_7810_days <- swap_sigma(mad_7810_days)
mad_lgr_days  <- swap_sigma(mad_lgr_days)
message("Campaign sigma per analyzer x day: closure median (reference) vs whole-record MAD (comparison)")
print(bind_rows(mad_7810_days, mad_lgr_days) %>% group_by(instrument) %>%
        summarise(days = n(), closure_median_CH4 = median(mad_sd_CH4, na.rm = TRUE),
                  whole_record_CH4 = median(wr_sd_CH4, na.rm = TRUE), .groups = "drop"))
mad_days <- bind_rows(mad_7810_days, mad_lgr_days) %>% filter(!is.na(mad_sd_CH4))
mad_inst_median <- mad_days %>%
  group_by(instrument) %>%
  summarise(mad_sd_CO2 = median(mad_sd_CO2, na.rm = TRUE),
            mad_sd_CH4 = median(mad_sd_CH4, na.rm = TRUE),
            dt_s = median(dt_s), .groups = "drop")
lgr_pooled <- mad_days %>% filter(instrument != "LI-7810") %>%
  summarise(mad_sd_CO2 = median(mad_sd_CO2, na.rm = TRUE),
            mad_sd_CH4 = median(mad_sd_CH4, na.rm = TRUE))

fl_machine <- field_logs %>% distinct(UniqueID, machine)
df <- df %>%
  mutate(instrument = if_else(is.na(instrument) & year == 2025, "LI-7810", instrument)) %>%
  left_join(fl_machine, by = "UniqueID") %>%
  mutate(instrument = if_else(is.na(instrument) & !is.na(machine), machine, instrument)) %>%
  select(-machine) %>%
  left_join(mad_days %>% select(instrument, date_str, sig_day_CO2 = mad_sd_CO2,
                                sig_day_CH4 = mad_sd_CH4, dt_day = dt_s) %>%
              # a day with two intervals has two rows; keep the one matching the closure
              group_by(instrument, date_str) %>% mutate(n_int = n()) %>% ungroup(),
            by = c("instrument", "date" = "date_str"), relationship = "many-to-many") %>%
  filter(is.na(n_int) | n_int == 1 | is.na(dt_s) | abs(dt_s - dt_day) < 0.2 * dt_day) %>%
  select(-n_int) %>%
  left_join(mad_inst_median %>% select(instrument, sig_inst_CO2 = mad_sd_CO2,
                                       sig_inst_CH4 = mad_sd_CH4, dt_inst = dt_s),
            by = "instrument") %>%
  mutate(
    sigma_camp_CO2 = coalesce(sig_day_CO2, sig_inst_CO2, lgr_pooled$mad_sd_CO2),
    sigma_camp_CH4 = coalesce(sig_day_CH4, sig_inst_CH4, lgr_pooled$mad_sd_CH4),
    sigma_camp_src = case_when(!is.na(sig_day_CH4) ~ "day",
                               !is.na(sig_inst_CH4) ~ "instrument median",
                               TRUE ~ "pooled LGR median"),
    # closures with no archived trace: interval and length are NOT measured;
    # keep them flagged so downstream can exclude them if the interval matters
    t_src = if_else(!is.na(t_sec), "measured", "assumed 300 s / 1 Hz"),
    dt_s = coalesce(dt_s, dt_day, dt_inst, 1)
  ) %>%
  select(-sig_day_CO2, -sig_day_CH4, -sig_inst_CO2, -sig_inst_CH4, -dt_day, -dt_inst)

message("\n=== Campaign sigma source (stem) ===")
print(table(df$instrument, df$sigma_camp_src, useNA = "ifany"))
message("=== Closure length / interval source (stem) ===")
print(table(df$t_src))
message("=== Logging intervals present, per analyzer x day (stem) ===")
print(mad_days %>% count(instrument, dt_s))

# ============================================================================
# PART 6: ROLLING WINDOW SCAN — CANOPY CONTINUOUS DATA
# ============================================================================

message("\n=== Rolling window scan for stable periods ===")

# Function: scan a continuous timeseries with rolling windows, scoring each
# by |slope| + residual_SD (both normalized), and return the SD from the
# quietest/flattest bottom 5% of windows.
rolling_window_precision <- function(df, col_name, window_sec = 30,
                                     quantile_cutoff = 0.05) {
  x <- df[[col_name]]
  t_sec <- as.numeric(difftime(df$POSIX.time, df$POSIX.time[1], units = "secs"))
  n <- length(x)

  if (n < window_sec + 1) {
    warning("Timeseries too short for rolling window (", n, " points)")
    return(list(sd = NA_real_, n_windows = 0, n_selected = 0))
  }

  # Pre-allocate
  n_windows <- n - window_sec + 1
  slopes    <- numeric(n_windows)
  resid_sds <- numeric(n_windows)

  for (i in seq_len(n_windows)) {
    idx <- i:(i + window_sec - 1)
    t_win <- t_sec[idx] - t_sec[idx[1]]
    x_win <- x[idx]

    # Quick linear regression via lm.fit for speed
    X <- cbind(1, t_win)
    fit <- .lm.fit(X, x_win)
    slopes[i] <- abs(fit$coefficients[2])
    resid_sds[i] <- sd(fit$residuals)
  }

  # Normalize both to [0, 1]
  slope_range <- range(slopes, na.rm = TRUE)
  sd_range    <- range(resid_sds, na.rm = TRUE)

  if (slope_range[2] - slope_range[1] == 0) {
    norm_slopes <- rep(0, n_windows)
  } else {
    norm_slopes <- (slopes - slope_range[1]) / (slope_range[2] - slope_range[1])
  }
  if (sd_range[2] - sd_range[1] == 0) {
    norm_sds <- rep(0, n_windows)
  } else {
    norm_sds <- (resid_sds - sd_range[1]) / (sd_range[2] - sd_range[1])
  }

  # Combined score: low = flat AND quiet
  combined_score <- norm_slopes + norm_sds

  # Select bottom 5%
  threshold <- quantile(combined_score, probs = quantile_cutoff, na.rm = TRUE)
  selected  <- which(combined_score <= threshold)
  n_selected <- length(selected)

  # Use the median residual SD of the selected windows as the precision estimate
  precision_sd <- median(resid_sds[selected])

  list(sd = precision_sd, n_windows = n_windows, n_selected = n_selected)
}

# Run for each HF instrument
message("Processing LGR1...")
rw_lgr1_co2 <- rolling_window_precision(imp.LGR1, "CO2dry_ppm")
rw_lgr1_ch4 <- rolling_window_precision(imp.LGR1, "CH4dry_ppb")

message("Processing LGR2...")
rw_lgr2_co2 <- rolling_window_precision(imp.LGR2, "CO2dry_ppm")
rw_lgr2_ch4 <- rolling_window_precision(imp.LGR2, "CH4dry_ppb")

message("Processing LGR3...")
rw_lgr3_co2 <- rolling_window_precision(imp.LGR3, "CO2dry_ppm")
rw_lgr3_ch4 <- rolling_window_precision(imp.LGR3, "CH4dry_ppb")

# YMF
if (exists("ymf_imp")) {
  message("Processing YMF...")
  rw_ymf_co2 <- rolling_window_precision(ymf_imp, "CO2dry_ppm")
  rw_ymf_ch4 <- rolling_window_precision(ymf_imp, "CH4dry_ppb")
} else {
  # Try to find YMF imp object by name pattern
  ymf_imp_name <- ls(pattern = "^imp\\.")
  ymf_imp_name <- ymf_imp_name[grepl("YMF", ymf_imp_name, ignore.case = TRUE)]
  if (length(ymf_imp_name) == 0) {
    ymf_imp_name <- ls(pattern = "^imp_YMF|^imp\\.YMF")
  }
  if (length(ymf_imp_name) > 0) {
    ymf_imp <- get(ymf_imp_name[1])
    message("Processing YMF (", ymf_imp_name[1], ")...")
    rw_ymf_co2 <- rolling_window_precision(ymf_imp, "CO2dry_ppm")
    rw_ymf_ch4 <- rolling_window_precision(ymf_imp, "CH4dry_ppb")
  } else {
    warning("YMF imported data not found; skipping YMF rolling window")
    rw_ymf_co2 <- list(sd = NA_real_, n_windows = 0, n_selected = 0)
    rw_ymf_ch4 <- list(sd = NA_real_, n_windows = 0, n_selected = 0)
  }
}

# Collect rolling window results
rw_results <- data.frame(
  instrument = c("LGR1", "LGR2", "LGR3", "YMF"),
  flat_sd_CO2 = c(rw_lgr1_co2$sd, rw_lgr2_co2$sd, rw_lgr3_co2$sd, rw_ymf_co2$sd),
  flat_sd_CH4 = c(rw_lgr1_ch4$sd, rw_lgr2_ch4$sd, rw_lgr3_ch4$sd, rw_ymf_ch4$sd),
  n_windows   = c(rw_lgr1_co2$n_windows, rw_lgr2_co2$n_windows,
                  rw_lgr3_co2$n_windows, rw_ymf_co2$n_windows),
  n_selected  = c(rw_lgr1_co2$n_selected, rw_lgr2_co2$n_selected,
                  rw_lgr3_co2$n_selected, rw_ymf_co2$n_selected),
  stringsAsFactors = FALSE
)

message("\n--- Rolling window precision estimates (bottom 5% quietest windows) ---")
print(rw_results)

message("\n--- Comparison: Manufacturer vs Empirical precision ---")
message("Manufacturer (GLA131): CO2 = ", PREC_CO2_LGR, " ppm, CH4 = ", PREC_CH4_LGR, " ppb")
for (i in seq_len(nrow(rw_results))) {
  message(rw_results$instrument[i], ": CO2 = ",
          round(rw_results$flat_sd_CO2[i], 4), " ppm, CH4 = ",
          round(rw_results$flat_sd_CH4[i], 4), " ppb")
}

# ============================================================================
# PART 7: SAVE ALL ALLAN DEVIATION RESULTS
# ============================================================================

message("\n=== Saving Allan deviation results ===")

save(df, n_total,
     allan_hf, allan_ymf,
     tsec_hf, tsec_ymf, mad_results, sigma_closures_canopy, closure_median_days,
     allan_df_7810, allan_df_lgr, closures_7810, lgr_catalog,
     mad_7810_days, mad_lgr_days, mad_days, mad_inst_median,
     field_logs, field_logs_lgr,
     rw_results,
     rw_lgr1_co2, rw_lgr1_ch4,
     rw_lgr2_co2, rw_lgr2_ch4,
     rw_lgr3_co2, rw_lgr3_ch4,
     rw_ymf_co2, rw_ymf_ch4,
     file = file.path(output_dir, "02_allan_deviation.RData"))

message("Saved: ", file.path(output_dir, "02_allan_deviation.RData"))

message("\n=== Allan deviation computation complete ===")
message("  HF canopy Allan:    ", nrow(allan_hf), " measurements")
message("  YMF canopy Allan:   ", nrow(allan_ymf), " measurements")
message("  LI-7810 stem Allan: ", nrow(allan_df_7810), " measurements")
message("  LGR stem Allan:     ", nrow(allan_df_lgr), " measurements")
message("  Rolling window:     ", nrow(rw_results), " instruments")
