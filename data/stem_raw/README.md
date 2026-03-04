# Raw 1-Hz instrument data

These files contain the raw concentration timeseries from which per-measurement
Allan deviations are computed in `02_allan_deviation.R`.

## 7810/
LI-COR LI-7810 pre-segmented CSVs from 2025 Harvard Forest stem flux measurements.
Each CSV covers one measurement (~60 seconds at 1 Hz).
Source: `tree-flux-2025/data/raw/7810_Processed/Tree_Fluxes/`

## LGR/
ABB GLA131/UGGA (LGR1, LGR2, LGR3) continuous recordings from 2023-2024
Harvard Forest stem flux campaigns. Each subdirectory contains date-organized
folders with full-day continuous recordings.
Source: `Matthes_Lab/stem-CH4-flux/Raw LGR Data/`

## Note
These files (~650 MB total) are only needed to recompute Allan deviations from
scratch. The intermediate results saved by `02_allan_deviation.R` to `outputs/`
are sufficient for all downstream analyses (scripts 03-08).
