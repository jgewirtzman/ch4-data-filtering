# Copy of: canopy (whole-tree) CH₄ flux closures, Harvard Forest and Yale-Myers Forest

**Source.** Gewirtzman et al., in review (the whole-tree stem, branch and leaf flux study). The
data will be archived with that paper; cite the paper and its data package for any use.

**Not new data.** These files are a minimal copy of the source study's processed closures,
reproduced here only so that this paper's analysis can be rerun. This paper uses them only to
measure the precision of the analyzer units (Fig. 3a and the precision comparison). They are not
new observations: obtain the data from the source, and do not count these closures as additional
or independent data in compilations, syntheses or meta-analyses.

## What is included, and why

Columns not read by the scripts (site, researchers, tree and chamber descriptors, environmental
covariates, notes, field times, goFlux diagnostics other than those below) are removed; values
are copied verbatim. In the `.RData` files the objects keep their source names.

| File | Contents | Rows kept | Columns kept |
|---|---|---|---|
| `canopy_flux_goFlux_compiled.csv` | goFlux fits of the Harvard Forest canopy closures (units LGR1, LGR2, LGR3) | all 136 | 13 of 60: `UniqueID`, `CO2_best.flux`, `CO2_model`, `CO2_LM.r2`, `CO2_HM.r2`, `CO2_nb.obs`, `CO2_flux.term`, `CH4_best.flux`, `CH4_model`, `CH4_LM.r2`, `CH4_HM.r2`, `CH4_LM.p.val`, `CH4_flux.term` |
| `ymf_black_oak_flux_compiled.csv` | goFlux fits of the Yale-Myers black-oak closures (unit LGR3, another day) | all 9 | the same 13 of 44 |
| `manID_LGR1.RData`, `manID_LGR2.RData`, `manID_LGR3.RData`, `manID_YMF.RData` | readings inside each closure's fitting window (`manID.LGR1`, …) | the window readings (`flag == 1`): 6812, 1067, 15316 and 950 | `UniqueID`, `flag`, `POSIX.time`, `CO2dry_ppm`, `CH4dry_ppb` |
| `imp_LGR1_combined.RData`, `imp_LGR2_combined.RData`, `imp_LGR3_combined.RData`, `imp_YMF_combined.RData` | the continuous record of each unit on the sampling days (`imp.LGR1`, …, `imp_YMF`) | all (27505, 3608, 45534, 1673) | `POSIX.time`, `CO2dry_ppm`, `CH4dry_ppb` |

Read by `01_load_data.R`; used by `02_allan_deviation.R` (per-closure and campaign precision,
closure length, rolling-window and whole-record precision) and `03_mdf_computation.R` (canopy
MDFs and `outputs/precision_comparison.csv`), and through them by 09, 16 and 19.

The paper reports units LGR1 and LGR3 (95 closures at 1 Hz; Fig. 3a). LGR2 (10 s logging) and
the Yale-Myers oak record are kept because they enter the intermediate results of 02 and 03
(precision table, canopy detection rates); removing them would change those saved outputs.

**Licence.** CC BY 4.0 (`LICENSE-data`).
