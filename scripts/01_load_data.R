# ============================================================================
# 01_load_data.R — Load and prepare all flux datasets
#
# Loads:
#   1. Stem flux data (HF 2023-2025, corrected) + chamber geometry
#   2. Canopy flux data (HF goFlux compiled)
#   3. YMF canopy flux data (Yale-Myers Forest black oak)
#   4. manID concentration timeseries (per-measurement, for Allan deviation)
#   5. Continuous imported data (for rolling window precision)
#   6. Soil flux data (semirigid chamber)
#
# Saves key objects to outputs/01_loaded_data.RData
# ============================================================================

source("scripts/00_setup.R")

# ============================================================================
# PART 1: LOAD CORRECTED STEM FLUX DATA
# ============================================================================

message("\n=== Loading corrected stem flux data ===")

fluxes <- read.csv(file.path(data_dir, "stem",
                              "HF_2023-2025_tree_flux_corrected.csv"),
                   stringsAsFactors = FALSE)
message("Loaded corrected flux file: ", nrow(fluxes), " rows")

df <- fluxes %>%
  filter(
    !is.na(year),
    !is.na(CO2_flux_umolpm2ps), !is.na(CO2_r2), !is.na(CO2_SE),
    !is.na(CH4_flux_nmolpm2ps), !is.na(CH4_r2), !is.na(CH4_SE)
  ) %>%
  mutate(
    # Pre-2025 CH4_SE needs x1000 correction (flux already corrected in file)
    CH4_SE_corr = if_else(year < 2025, CH4_SE * 1000, CH4_SE),
    # SE-based SNR
    CO2_snr_se = abs(CO2_flux_umolpm2ps) / CO2_SE,
    CH4_snr_se = abs(CH4_flux_nmolpm2ps) / CH4_SE_corr
  )

n_total <- nrow(df)
message("After NA removal: ", n_total, " measurements")
message("  2023-24 (LGR/UGGA): ", sum(df$year < 2025),
        " | 2025 (LI-7810): ", sum(df$year == 2025))

# ============================================================================
# PART 2: CHAMBER GEOMETRY FROM tree_volumes.csv
# ============================================================================

message("\n=== Loading chamber geometry ===")

tv <- read.csv(file.path(data_dir, "stem", "tree_volumes.csv"),
               stringsAsFactors = FALSE)
tv$Tag_int <- as.integer(round(tv$Tag))
# vol_system = Volume (L) + 0.028 L extra connecting tubing
tv$vol_system_tv <- tv$Volume + EXTRA_TUBE_VOL

message("tree_volumes.csv: ", nrow(tv), " trees")

# Join tree_volumes to corrected flux by Tree = Tag_int
df <- df %>%
  left_join(tv[, c("Tag_int", "vol_system_tv")],
            by = c("Tree" = "Tag_int"))

# Fill vol_system for ALL rows using tree_volumes.csv
# (2025 data already has vol_system; 2023-24 is NA)
df <- df %>%
  mutate(
    vol_system = if_else(is.na(vol_system), vol_system_tv, vol_system),
    surfarea   = if_else(is.na(surfarea), SURFAREA_M2, surfarea),
    nmol       = if_else(is.na(nmol) & !is.na(bar) & !is.na(airt) & !is.na(vol_system),
                          bar * vol_system / (R_hPa_L * (airt + 273.15)),
                          nmol)
  ) %>%
  select(-vol_system_tv)

message("Chamber geometry filled:")
message("  vol_system NAs: ", sum(is.na(df$vol_system)), " / ", n_total)
message("  nmol NAs:       ", sum(is.na(df$nmol)), " / ", n_total)
message("  surfarea NAs:   ", sum(is.na(df$surfarea)), " / ", n_total)

# ============================================================================
# PART 3: LOAD CANOPY FLUX DATA (Harvard Forest)
# ============================================================================

message("\n=== Loading canopy flux data ===")

hf_compiled <- read.csv(
  file.path(data_dir, "canopy", "canopy_flux_goFlux_compiled.csv"),
  stringsAsFactors = FALSE)
hf_compiled$UniqueID <- trimws(hf_compiled$UniqueID)

message("HF compiled: ", nrow(hf_compiled), " rows")

# ============================================================================
# PART 4: LOAD YMF CANOPY DATA (Yale-Myers Forest)
# ============================================================================

message("\n=== Loading YMF canopy flux data ===")

ymf_compiled <- read.csv(
  file.path(data_dir, "ymf_canopy", "ymf_black_oak_flux_compiled.csv"),
  stringsAsFactors = FALSE)
ymf_compiled$UniqueID <- trimws(ymf_compiled$UniqueID)

message("YMF compiled: ", nrow(ymf_compiled), " rows")

# ============================================================================
# PART 5: LOAD manID DATA (per-measurement concentration timeseries)
# ============================================================================

message("\n=== Loading manID data ===")

load(file.path(data_dir, "canopy", "manID_LGR1.RData"))
load(file.path(data_dir, "canopy", "manID_LGR2.RData"))
load(file.path(data_dir, "canopy", "manID_LGR3.RData"))
load(file.path(data_dir, "ymf_canopy", "manID_YMF.RData"))

manID_all_hf <- rbind(manID.LGR1, manID.LGR2, manID.LGR3)
message("HF manID rows: ", nrow(manID_all_hf),
        " (", length(unique(manID_all_hf$UniqueID)), " measurements)")
message("YMF manID rows: ", nrow(manID.YMF),
        " (", length(unique(manID.YMF$UniqueID)), " measurements)")

# Determine which instrument each HF measurement came from
# UniqueIDs in manID data tell us the mapping
lgr1_uids <- unique(manID.LGR1$UniqueID)
lgr2_uids <- unique(manID.LGR2$UniqueID)
lgr3_uids <- unique(manID.LGR3$UniqueID)

hf_compiled$instrument <- NA_character_
hf_compiled$instrument[hf_compiled$UniqueID %in% lgr1_uids] <- "LGR1"
hf_compiled$instrument[hf_compiled$UniqueID %in% lgr2_uids] <- "LGR2"
hf_compiled$instrument[hf_compiled$UniqueID %in% lgr3_uids] <- "LGR3"
ymf_compiled$instrument <- "YMF"

# ============================================================================
# PART 6: LOAD CONTINUOUS IMPORTED DATA (for rolling window precision)
# ============================================================================

message("\n=== Loading continuous imported data ===")

load(file.path(data_dir, "canopy", "imp_LGR1_combined.RData"))
load(file.path(data_dir, "canopy", "imp_LGR2_combined.RData"))
load(file.path(data_dir, "canopy", "imp_LGR3_combined.RData"))
load(file.path(data_dir, "ymf_canopy", "imp_YMF_combined.RData"))

# Tag each with instrument name
imp.LGR1$instrument <- "LGR1"
imp.LGR2$instrument <- "LGR2"
imp.LGR3$instrument <- "LGR3"

message("  imp.LGR1: ", nrow(imp.LGR1), " rows")
message("  imp.LGR2: ", nrow(imp.LGR2), " rows")
message("  imp.LGR3: ", nrow(imp.LGR3), " rows")

# Detect YMF imported object name
ymf_imp_name <- ls(pattern = "^imp\\.")
ymf_imp_name <- ymf_imp_name[grepl("YMF", ymf_imp_name, ignore.case = TRUE)]
if (length(ymf_imp_name) == 0) {
  ymf_imp_name <- ls(pattern = "^imp_YMF|^imp\\.YMF")
}
if (length(ymf_imp_name) > 0) {
  ymf_imp <- get(ymf_imp_name[1])
  message("  ", ymf_imp_name[1], ": ", nrow(ymf_imp), " rows")
} else {
  message("  Warning: YMF imported data object not found")
}

# ============================================================================
# PART 7: LOAD SOIL FLUX DATA
# ============================================================================

message("\n=== Loading soil flux data ===")

soil_flux <- read.csv(
  file.path(data_dir, "soil", "semirigid_chamber_flux.csv"),
  stringsAsFactors = FALSE)

message("Soil flux: ", nrow(soil_flux), " rows")

# ============================================================================
# PART 8: SAVE ALL LOADED DATA
# ============================================================================

message("\n=== Saving loaded data to outputs/ ===")

save(df, n_total,
     hf_compiled, ymf_compiled,
     manID_all_hf, manID.LGR1, manID.LGR2, manID.LGR3, manID.YMF,
     lgr1_uids, lgr2_uids, lgr3_uids,
     imp.LGR1, imp.LGR2, imp.LGR3,
     soil_flux,
     file = file.path(output_dir, "01_loaded_data.RData"))

# Also save YMF imp object if it exists
if (exists("ymf_imp")) {
  save(ymf_imp, file = file.path(output_dir, "01_ymf_imp.RData"))
}

message("Saved: ", file.path(output_dir, "01_loaded_data.RData"))

# ============================================================================
# PART 9: SUMMARY
# ============================================================================

message("\n=== Data loading summary ===")
message("  Stem flux (df):        ", nrow(df), " measurements (",
        sum(df$year < 2025), " LGR + ", sum(df$year == 2025), " LI-7810)")
message("  HF canopy compiled:    ", nrow(hf_compiled), " measurements")
message("  YMF canopy compiled:   ", nrow(ymf_compiled), " measurements")
message("  HF manID timeseries:   ", nrow(manID_all_hf), " rows (",
        length(unique(manID_all_hf$UniqueID)), " UIDs)")
message("  YMF manID timeseries:  ", nrow(manID.YMF), " rows (",
        length(unique(manID.YMF$UniqueID)), " UIDs)")
message("  Soil flux:             ", nrow(soil_flux), " measurements")
message("\nAll data loaded successfully.")
