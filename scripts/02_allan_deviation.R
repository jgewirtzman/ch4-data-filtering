# ============================================================================
# 02_allan_deviation.R — Compute Allan deviation for all flux datasets
#
# Computes per-measurement instrument noise (Allan deviation at tau=1s) from:
#   1. HF canopy manID timeseries (LGR1/2/3)
#   2. YMF canopy manID timeseries
#   3. 2025 LI-7810 stem raw 1-Hz CSVs
#   4. 2023-24 LGR stem raw recordings (segmented by field log timestamps)
#
# Also runs rolling window scan on continuous canopy data to estimate
# global per-instrument precision.
#
# Merges Allan results into stem flux data frame (df).
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
# PART 2: LOAD FIELD LOGS (for LGR measurement timestamps)
# ============================================================================

message("\n=== Loading field logs for LGR measurement timestamps ===")

field_data_dir <- file.path(data_dir, "stem", "field_logs")

# --- Field log 1: Summer 2023 ---
fl1 <- read.csv(file.path(field_data_dir, "Field_Data_Monthly_Summer2023.csv"),
                stringsAsFactors = FALSE)
fl1_std <- fl1 %>%
  transmute(
    UniqueID       = UniqueID,
    date_raw       = Date,
    tree_raw       = as.character(Tree_Tag),
    comp_start     = comp_start_time,
    comp_end       = comp_end_time,
    real_start     = Real.start,
    machine        = trimws(Analyzer),
    source         = "summer2023"
  )

# --- Field log 2: Updated2 (Sep 2023 - Sep 2024) ---
fl2 <- read.csv(file.path(field_data_dir, "Field_Data_Monthly_updated2.csv"),
                stringsAsFactors = FALSE)
# Use the "Updated Start Time" / "Updated End Time" columns if available
# (these incorporate QC corrections), else fall back to comp_start_time
fl2_std <- fl2 %>%
  transmute(
    UniqueID       = UniqueID,
    date_raw       = format_Date,
    tree_raw       = as.character(Tree.Tag),
    comp_start     = if_else(nchar(trimws(Updated.Start.Time)) > 0,
                              Updated.Start.Time, comp_start_time),
    comp_end       = if_else(nchar(trimws(Updated.End.Time)) > 0,
                              Updated.End.Time, comp_end_time),
    real_start     = Real.start,
    machine        = trimws(Machine),
    source         = "updated2"
  )

# --- Field log 3: Summer 2024 ---
fl3 <- read.csv(file.path(field_data_dir, "summer_2024.csv"),
                stringsAsFactors = FALSE)
fl3_std <- fl3 %>%
  transmute(
    UniqueID       = UniqueID,
    date_raw       = format_Date,
    tree_raw       = as.character(Tree),
    comp_start     = comp_start_time,
    comp_end       = comp_end_time,
    real_start     = Real.start,
    machine        = "LGR1",   # no Machine column; all 2024 used LGR1
    source         = "summer2024"
  )

# --- Field log 4: Sept-Dec 2024 Excel files ---
xlsx_dir <- file.path(field_data_dir, "Sept2024_onwards")
fl4_std <- data.frame()
if (requireNamespace("readxl", quietly = TRUE) && dir.exists(xlsx_dir)) {
  xlsx_files <- list.files(xlsx_dir, pattern = "[.]xlsx$", full.names = TRUE)
  xlsx_files <- xlsx_files[!grepl("Template", xlsx_files)]
  xlsx_files <- xlsx_files[!grepl("2025", xlsx_files)]  # skip 2025 (LI-7810)

  for (xf in xlsx_files) {
    tryCatch({
      xd <- as.data.frame(readxl::read_excel(xf))
      # Normalize column names (spaces -> dots)
      names(xd) <- make.names(names(xd))
      xd_std <- xd %>%
        transmute(
          UniqueID   = as.character(UniqueID),
          date_raw   = as.character(format_Date),
          tree_raw   = as.character(Tree.Tag),
          comp_start = as.character(comp_start_time),
          comp_end   = as.character(comp_end_time),
          real_start = as.character(Real.start),
          machine    = if ("Machine" %in% names(xd))
                         as.character(Machine) else "LGR1",
          source     = "xlsx_late2024"
        )
      fl4_std <- bind_rows(fl4_std, xd_std)
    }, error = function(e) {
      message("  Warning: could not read ", basename(xf), ": ", e$message)
    })
  }
  # Clean up comp_start/comp_end that may have Excel datetime artifacts
  # (e.g., "1899-12-31 12:18:52" -> extract just "12:18:52")
  fl4_std$comp_start <- sub("^.*\\s+", "", fl4_std$comp_start)
  fl4_std$comp_end   <- sub("^.*\\s+", "", fl4_std$comp_end)
  fl4_std$real_start <- sub("^.*\\s+", "", fl4_std$real_start)

  message("Excel field logs (Sept-Dec 2024): ", nrow(fl4_std), " entries")
} else {
  message("readxl not available or xlsx dir missing; skipping Excel field logs")
}

# Combine all field logs
field_logs <- bind_rows(fl1_std, fl2_std, fl3_std, fl4_std)
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
# Handle ISO format "2024-09-17" (from readxl Date objects converted to char)
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
    comp_start_posix = as.POSIXct(paste(date_parsed, comp_start),
                                    format = "%Y-%m-%d %H:%M:%S"),
    comp_end_posix   = as.POSIXct(paste(date_parsed, comp_end),
                                    format = "%Y-%m-%d %H:%M:%S"),
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

raw_dir_7810 <- file.path(data_dir, "stem_raw", "7810")
raw_files_7810 <- list.files(raw_dir_7810, pattern = "\\.csv$",
                              recursive = TRUE, full.names = TRUE)
raw_files_7810 <- raw_files_7810[grepl("/tables/", raw_files_7810)]
message("Found ", length(raw_files_7810), " raw 1-Hz LI-7810 files")

allan_7810 <- lapply(raw_files_7810, function(f) {
  tryCatch({
    raw <- read.csv(f, stringsAsFactors = FALSE)
    if (nrow(raw) < 5) return(NULL)
    data.frame(
      REMARK       = raw$REMARK[1],
      DATE         = raw$DATE[1],
      allan_sd_CO2 = allan_sd(raw$CO2),
      allan_sd_CH4 = allan_sd(raw$CH4),
      n_pts        = nrow(raw),
      t_sec        = nrow(raw),
      instrument   = "LI-7810",
      stringsAsFactors = FALSE
    )
  }, error = function(e) NULL)
})

allan_df_7810 <- do.call(rbind, Filter(Negate(is.null), allan_7810))

if (is.null(allan_df_7810) || nrow(allan_df_7810) == 0) {
  warning("No LI-7810 Allan results. Continuing with LGR only.")
  allan_df_7810 <- data.frame(
    REMARK = character(), DATE = character(),
    allan_sd_CO2 = numeric(), allan_sd_CH4 = numeric(),
    n_pts = integer(), t_sec = integer(), instrument = character(),
    stringsAsFactors = FALSE
  )
} else {
  message("LI-7810 Allan deviation: ", nrow(allan_df_7810), " measurements")
  message("  CH4 Allan SD: median ", round(median(allan_df_7810$allan_sd_CH4, na.rm = TRUE), 4),
          " ppb, range ", round(min(allan_df_7810$allan_sd_CH4, na.rm = TRUE), 4),
          " - ", round(max(allan_df_7810$allan_sd_CH4, na.rm = TRUE), 4))
}

# Build match key for 7810: REMARK_DATE
allan_df_7810$match_key <- paste(allan_df_7810$REMARK, allan_df_7810$DATE, sep = "_")

# ============================================================================
# PART 4: ALLAN DEVIATION — 2023-24 LGR/UGGA STEM DATA
# ============================================================================

message("\n=== Allan deviation: 2023-24 LGR/UGGA ===")

lgr_base <- file.path(data_dir, "stem_raw", "LGR")
message("LGR raw data directory: ", lgr_base)
message("  Exists: ", dir.exists(lgr_base))

# Catalog all available LGR date folders
lgr_catalog <- data.frame(machine = character(), date_str = character(),
                           folder = character(), stringsAsFactors = FALSE)
for (lgr_name in c("LGR1", "LGR2", "LGR3")) {
  lgr_dir <- file.path(lgr_base, lgr_name)
  if (!dir.exists(lgr_dir)) next
  date_folders <- list.dirs(lgr_dir, recursive = FALSE, full.names = TRUE)
  for (folder in date_folders) {
    folder_name <- basename(folder)
    # Folder names like "2023-06-29" or "2023-06-06 (1)"
    date_part <- substr(folder_name, 1, 10)
    lgr_catalog <- rbind(lgr_catalog, data.frame(
      machine = lgr_name, date_str = date_part, folder = folder,
      stringsAsFactors = FALSE
    ))
  }
}
message("LGR catalog: ", nrow(lgr_catalog), " date-folders across ",
        length(unique(lgr_catalog$machine)), " instruments")

# Parse one LGR raw file (txt or txt.zip) and return 1-Hz data frame
parse_lgr_file <- function(filepath) {
  tryCatch({
    if (grepl("\\.zip$", filepath)) {
      # Read from zip without extracting
      txt_name <- sub("\\.zip$", "", basename(filepath))
      con <- unz(filepath, txt_name)
      lines <- readLines(con, warn = FALSE)
      close(con)
    } else {
      lines <- readLines(filepath, warn = FALSE)
    }

    if (length(lines) < 3) return(NULL)

    # Skip line 1 (serial/metadata), read from line 2 (headers) onwards
    dat <- tryCatch(
      read.csv(text = paste(lines[-1], collapse = "\n"),
               stringsAsFactors = FALSE, strip.white = TRUE),
      error = function(e) NULL
    )
    if (is.null(dat) || nrow(dat) == 0) return(NULL)

    # Parse timestamp from first column (SysTime)
    # Format: "06/29/2023 10:06:11.087"
    ts_col <- names(dat)[1]  # usually "SysTime" or similar
    dat$timestamp <- as.POSIXct(dat[[ts_col]], format = "%m/%d/%Y %H:%M:%OS")

    # Get CH4 and CO2 dry mole fractions
    ch4_col <- grep("CH4.*d_ppm", names(dat), value = TRUE)[1]
    co2_col <- grep("CO2.*d_ppm", names(dat), value = TRUE)[1]
    if (is.na(ch4_col) || is.na(co2_col)) return(NULL)

    data.frame(
      timestamp = dat$timestamp,
      CH4_ppm   = as.numeric(dat[[ch4_col]]),
      CO2_ppm   = as.numeric(dat[[co2_col]]),
      stringsAsFactors = FALSE
    ) %>% filter(!is.na(timestamp))
  }, error = function(e) {
    NULL
  })
}

# Load and concatenate all LGR files for a given folder
load_lgr_day <- function(folder_path) {
  # Find all txt and txt.zip files with _f prefix (data files, not _b, _l, _p)
  all_files <- list.files(folder_path, full.names = TRUE)
  data_files <- all_files[grepl("_f\\d+\\.txt(\\.zip)?$", all_files)]

  if (length(data_files) == 0) return(NULL)

  # Parse each file and concatenate
  day_data <- lapply(data_files, parse_lgr_file)
  day_data <- do.call(rbind, Filter(Negate(is.null), day_data))

  if (is.null(day_data) || nrow(day_data) == 0) return(NULL)

  # Sort by timestamp
  day_data <- day_data[order(day_data$timestamp), ]
  day_data
}

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

  # Get all field log entries for this date+machine
  entries <- field_logs_lgr %>%
    filter(date_str == d, machine == m)

  # Segment each measurement and compute Allan deviation
  for (j in seq_len(nrow(entries))) {
    start_t <- entries$comp_start_posix[j]
    end_t   <- entries$comp_end_posix[j]

    segment <- day_data[day_data$timestamp >= start_t &
                          day_data$timestamp <= end_t, ]

    if (nrow(segment) < 5) next

    n_matched <- n_matched + 1
    allan_lgr_all[[n_matched]] <- data.frame(
      UniqueID     = entries$UniqueID[j],
      date_str     = d,
      allan_sd_CO2 = allan_sd(segment$CO2_ppm),          # ppm (same as LI-7810)
      allan_sd_CH4 = allan_sd(segment$CH4_ppm) * 1000,  # convert ppm -> ppb to match LI-7810
      n_pts        = nrow(segment),
      t_sec        = nrow(segment),
      instrument   = m,
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

if (is.null(allan_df_lgr) || nrow(allan_df_lgr) == 0) {
  warning("No LGR Allan results computed.")
  allan_df_lgr <- data.frame(
    UniqueID = character(), date_str = character(),
    allan_sd_CO2 = numeric(), allan_sd_CH4 = numeric(),
    n_pts = integer(), t_sec = integer(), instrument = character(),
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

# --- Match 2025 LI-7810 by Tree_Date key ---
df$match_key_7810 <- paste(df$Tree, df$date, sep = "_")

allan_7810_unique <- allan_df_7810 %>%
  group_by(match_key) %>%
  slice(1) %>%
  ungroup() %>%
  select(match_key, allan_sd_CO2, allan_sd_CH4, n_pts, t_sec, instrument)

df <- merge(df, allan_7810_unique,
            by.x = "match_key_7810", by.y = "match_key",
            all.x = TRUE, suffixes = c("", "_7810"))

# --- Match 2023-24 LGR by UniqueID ---
if (nrow(allan_df_lgr) > 0) {
  allan_lgr_unique <- allan_df_lgr %>%
    group_by(UniqueID) %>%
    slice(1) %>%
    ungroup() %>%
    select(UniqueID, allan_sd_CO2, allan_sd_CH4, n_pts, t_sec, instrument)

  df <- merge(df, allan_lgr_unique,
              by = "UniqueID", all.x = TRUE,
              suffixes = c("", "_lgr"))

  # Combine: prefer existing (7810) values, fill in LGR where missing
  df <- df %>%
    mutate(
      allan_sd_CO2 = if_else(is.na(allan_sd_CO2), allan_sd_CO2_lgr, allan_sd_CO2),
      allan_sd_CH4 = if_else(is.na(allan_sd_CH4), allan_sd_CH4_lgr, allan_sd_CH4),
      n_pts        = if_else(is.na(n_pts), n_pts_lgr, n_pts),
      t_sec        = if_else(is.na(t_sec), t_sec_lgr, t_sec),
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
                                     format = "%Y-%m-%d %H:%M:%S")
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
    df$n_pts[idx]        <- allan_row$n_pts[1]
    df$t_sec[idx]        <- allan_row$t_sec[1]
    df$instrument[idx]   <- allan_row$instrument[1]
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
     allan_df_7810, allan_df_lgr,
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
