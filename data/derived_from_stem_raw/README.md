# Derived from the stem study's raw records (intermediate results of scripts 02 and 02b)

**What these are.** The saved results of `scripts/02_allan_deviation.R` and
`scripts/02b_goflux_fluxes.R`, computed from the raw 1 Hz analyzer records of the Harvard Forest
stem study (`data/stem_raw/`, ≈650 MB, not in this repository) together with
`data/from_hf_stem_study/` and `data/from_canopy_study/`. They are shipped so that scripts 03
onward run without the raw archive.

**Not new data.** They are a processing of the stem study's data (J. Gewirtzman et al., in prep.;
raw records to be archived in the Harvard Forest Data Archive, DOI to follow), not a separate
dataset and not new observations. Cite the stem study and this paper; do not count these closures
as additional or independent data in compilations, syntheses or meta-analyses.

| File | Written by | Main contents | Needed by |
|---|---|---|---|
| `02_allan_deviation.RData` | 02 | per-closure precision (second-difference σ, Allan deviation), closure length and logging interval, campaign σ per analyzer and day, slope p-values, the canopy precision results (`mad_results`, `rw_results`, `allan_hf`, `tsec_hf`, …), and the stem table `df` with these columns | 03, 17 |
| `02b_goflux_fluxes.RData` | 02b | goFlux refits of the 1547 stem closures with a raw record (`ch4_best`, `co2_best`), campaign σ per group (`sig_grp`), closure windows (`windows`), and `df` with the refitted fluxes, detection classes and trace screens | 03, 17 |
| `02b_manID.RData` | 02b | the raw readings inside every fitted stem closure window (`manID`, 477,255 rows) | 17, 21 |

**How the pipeline uses them.** When the raw records are absent (`data/stem_raw/LGR/LGR1` and
`data/stem_raw/7810` not found), `scripts/00_setup.R` copies any of these files that is missing
from `outputs/` into `outputs/`, and `scripts/run_all.R` skips 02, 02b and 14. With the raw records
present, 02 and 02b recompute them and `run_all.R` refreshes the copies here, so that these files
always match the last full run.

**Licence.** CC BY 4.0 (`LICENSE-data`).
