# Copy of: Harvard Forest tree-stem CH₄ flux study (2023–2025)

**Source.** J. Gewirtzman et al., in preparation. The closure table, chamber volumes, field logs
and raw analyzer records will be archived in the Harvard Forest Data Archive (DOI to follow).
Cite that package, and the source paper when published, for any use of these data.

**Not new data.** These files are a minimal copy of the stem study's tables, reproduced here only
so that this paper's analysis can be rerun. They are not new observations. Obtain the data from
the source archive, and do not count these closures as additional or independent data in
compilations, syntheses or meta-analyses.

## What is included, and why

Only the columns and rows that the scripts read; every other column of the source tables
(e.g. relative humidity, PAR, concentrations, error notes, DBH, sampling-round bookkeeping,
soil moisture, plot notes, operator notes) is removed. Values are copied verbatim (no rounding).

| File | Rows | Columns kept | Read by |
|---|---|---|---|
| `HF_2023-2025_tree_flux_corrected.csv` | all 1640 closures of the source table | 20 of 38: `airt`, `bar`, `datetime_posx`, `Tree`, `PLOT`, `SPECIES`, `vol_system`, `nmol`, `surfarea`, `CO2_r2`, `CH4_r2`, `CO2_SE`, `CH4_SE`, `CO2_flux_umolpm2ps`, `CH4_flux_nmolpm2ps`, `year`, `UniqueID`, `date`, `location`, `species_label` | `01_load_data.R` |
| `tree_volumes.csv` | all 60 stem collars | `Tag` (tree tag), `Volume` (collar + cap + tubing volume, L) | `01_load_data.R` |
| `field_log_closure_times.csv` | 1437 closure entries, June 2023 – December 2024 (LGR/UGGA campaigns) | `UniqueID`, `date_raw`, `tree_raw`, `comp_start`, `comp_end`, `real_start`, `machine`, `source` | `02_allan_deviation.R` (only with the raw records) |
| `figure_closures/` | readings of the 36 closures drawn in Figs. 2, 5 and S1 | see `figure_closures/README.md` | `15_co2_validity_figure.R`, `17_sigma_figures.R` |

`field_log_closure_times.csv` combines the four field logs of the 2023–2024 campaigns
(`Field_Data_Monthly_Summer2023.csv`, `Field_Data_Monthly_updated2.csv`, `summer_2024.csv` and the
six Sept–Dec 2024 `Stem_flux_*.xlsx` sheets) into the columns that `02_allan_deviation.R` uses to
locate each closure in the raw records, as text exactly as entered (or, for the `.xlsx` sheets,
as read by `readxl` with the spurious Excel date stripped from the clock times): closure ID, date,
tree tag, start and end of the computation window and the analyzer's start time on the analyzer
clock, analyzer unit (`machine`: as logged, or `LGR1` for logs without an analyzer column; empty
or `NA` entries are read as LGR1 by 02, the only unit used for stems in 2023–2024), and which log the row came from (`source`: `summer2023`, `updated2`, `summer2024`,
`xlsx_late2024`; 02 keeps the later log when a closure appears in two). For the
September 2023 – September 2024 log (`updated2`), `comp_start`/`comp_end` are the QC-updated window
times where the log has them.

The raw 1 Hz records are not in this repository (`data/stem_raw/README.md`). The results that
02 and 02b compute from them are in `data/derived_from_stem_raw/`.

**Licence.** CC BY 4.0 (`LICENSE-data`), as for the source tables.
