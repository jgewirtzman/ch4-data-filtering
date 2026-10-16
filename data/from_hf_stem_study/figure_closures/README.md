# Readings of the closures drawn in the figures

Small extracts of the raw stem records (`data/stem_raw/`, not in this repository), written by
`scripts/14_figure_closures.R` so that `15_co2_validity_figure.R` and `17_sigma_figures.R` run
without the raw archive. Like the rest of `data/from_hf_stem_study/`, they are a copy of the stem
study's data (J. Gewirtzman et al., in prep.; Harvard Forest Data Archive, DOI to follow), not
new observations; see `../README.md`. Each file holds only the readings inside the closure windows the
figures show. Numbers are written with enough digits to read back identically, so the figures
built from these files match those built from the raw records.

`time_s` is the analyzer clock time in seconds since 1970-01-01 00:00:00, labelled UTC
(as everywhere in the stem pipeline, it is the analyzer's local clock, the clock the field logs
refer to; in R, `.POSIXct(time_s, tz = "UTC")`). Read numeric columns as doubles
(`colClasses = "numeric"`).

| File | Figure | Closures | Readings |
|---|---|---|---|
| `fig5_co2_tracer.csv` | Fig. 5 (CO₂ as a validity tracer) | 4 closures, ABB/LGR GLA131 unit LGR1, 2024 | 724 |
| `fig2_sigma_day.csv` | Fig. 2 (how σ is estimated) | all 30 closures of 3 June 2025, LI-COR LI-7810 | 6068 |
| `worked_example_closures.csv` | Sect. S1 worked example, Fig. S1, `figures/worked_example.json` | 2 closures, LGR1, 14 Nov 2024 | 362 |

## fig5_co2_tracer.csv

The 180 s computation window of each closure, from the analyzer's `Time` column and dry mole
fractions.

| Column | Unit | Description |
|---|---|---|
| `panel` | | Column of Fig. 5 (1 = a, ..., 4 = d) |
| `UniqueID` | | Field-log closure ID |
| `analyzer` | | Analyzer unit |
| `time_s` | s | Reading time (see above) |
| `CO2d_ppm` | ppm (µmol mol⁻¹) | Dry CO₂ mole fraction, as logged |
| `CH4d_ppm` | ppm (µmol mol⁻¹) | Dry CH₄ mole fraction, as logged (Fig. 5 shows ppb = 1000 × ppm) |

## fig2_sigma_day.csv

The closure windows (flag = 1) that `02b_goflux_fluxes.R` extracted for every closure of the
LI-7810 on 2025-06-03; panels a-c show three of them, panel d all 30.

| Column | Unit | Description |
|---|---|---|
| `UniqueID` | | Stem closure ID (`closure_id` in `results/stem_closures_processed.csv`) |
| `analyzer` | | Analyzer |
| `time_s` | s | Reading time (see above) |
| `CH4dry_ppb` | ppb (nmol mol⁻¹) | Dry CH₄ mole fraction |

## worked_example_closures.csv

The two closures of the worked example (a flux near the MDF and a large flux; chosen by
`pick_worked_example()` in `scripts/00_setup.R`), over the field-log computation window used by
`02_allan_deviation.R`, from the analyzer's `SysTime` column.

| Column | Unit | Description |
|---|---|---|
| `UniqueID` | | Stem closure ID |
| `fl_UniqueID` | | Field-log closure ID |
| `analyzer` | | Analyzer unit |
| `time_s` | s | Reading time (see above) |
| `CH4d_ppm` | ppm (µmol mol⁻¹) | Dry CH₄ mole fraction, as logged (17 uses ppb = 1000 × ppm) |

Licence: CC BY 4.0 (`LICENSE-data`).
