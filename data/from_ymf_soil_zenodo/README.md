# Copy of: Yale-Myers Forest upland–wetland soil CH₄ fluxes (Zenodo)

**Source.** Gewirtzman, J. (2026). jgewirtzman/tree-methanogens: Tree Microbiomes and Methane
Exchange in Upland Forests — Code & Data (v0.1.0). Zenodo. https://doi.org/10.5281/zenodo.18779715.
Cite that record for any use of these data.

**Not new data.** These files are a minimal copy of three files of that record, reproduced here
only so that this paper's analysis can be rerun. They are not new observations: obtain the data
from the Zenodo record, and do not count these closures as additional or independent data in
compilations, syntheses or meta-analyses.

## What is included, and why

Values are copied verbatim (no rounding); only rows and columns are removed.

| File | Source file in the record | Rows kept | Columns kept | Read by |
|---|---|---|---|---|
| `semirigid_chamber_flux.csv` | `semirigid_chamber_flux.csv` | the 296 soil closures (`measurement_type == "soil"`) of 692 (the others are tree stems, not used) | 8 of 21: `unique_id`, `date`, `landscape_position`, `chamber_temp_C`, `CH4_best_flux_nmol_m2_s`, `CH4_LM_r2`, `CH4_LM_pval`, `measurement_type` | `01_load_data.R` → 08, 10 → 19, 20 |
| `CH4_best_flux_lgr_results_soil.csv` | `CH4_best_flux_lgr_results_soil.csv` (goFlux fits) | all 288 | 3 of 48: `UniqueID`, `flux.term`, `best.flux` | `20_generality.R` |
| `lgr_manual_identification_results_soil_flag1.csv.gz` | `lgr_manual_identification_results_soil.csv` (closure windows, 75 MB) | the 75,038 readings inside the fitting windows (`flag == 1`) of 257,100 | 4 of 29: `UniqueID`, `POSIX.time`, `CH4dry_ppb`, `flag` | `20_generality.R` (soil σ and MDF) |

The windows file is gzip-compressed (0.4 MB); `data.table::fread` reads it directly (needs the
R.utils package). It drops, among others, the sampling-team and file-name columns. To run 20
on the complete source file instead, download it from the record and set the environment variable
`SOIL_WINDOWS` to its path (for example in `data/external/`, which git ignores); the results are
identical.

**Licence.** CC BY 4.0, as published on Zenodo (the record's licence), with attribution to the
record above.
