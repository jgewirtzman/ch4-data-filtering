# Input data: copies of other datasets, not new observations

**Every dataset in this folder was collected and is published (or will be published) by
another project.** The files here are minimal copies, reduced to the rows and columns this
paper's scripts read, and reproduced only so that the analysis can be rerun from this
repository alone. They are not new observations and add no data to the source studies.

- To use the data, obtain them from the original source and cite that source (table below).
- Do not count these files as additional or independent data in compilations, syntheses or
  meta-analyses: each closure is already in its source dataset, and counting it here as well
  would count it twice.
- `derived_from_stem_raw/` holds intermediate results that this paper's scripts computed from
  the stem study's raw records: a processing of that dataset, not a separate one.

Each folder's README states the source, what was kept and why, and the licence.

| Folder | Source | Citation or archive | Used by | Licence |
|---|---|---|---|---|
| `from_hf_stem_study/` | Harvard Forest tree-stem CH₄ fluxes, 2023–2025: closure table, chamber volumes, field-log closure times; `figure_closures/` holds the raw readings of the few closures drawn in Figs. 2, 5 and S1 | J. Gewirtzman et al., in prep.; to be archived with the raw records in the Harvard Forest Data Archive (DOI to follow) | 01, 02, 15, 17 | CC BY 4.0 |
| `stem_raw/` | Raw 1 Hz analyzer records of the stem study (≈650 MB; **not in this repository**, local only) | same as above | 02, 02b, 14 | (source archive) |
| `derived_from_stem_raw/` | Results of scripts 02 and 02b computed from `stem_raw/`, so that 03 onward run without the raw records | this repository (processing of the stem study) | 03–23 | CC BY 4.0 |
| `from_canopy_study/` | Harvard Forest and Yale-Myers canopy closures (stems, branches, leaves): fluxes and processed concentration traces of three analyzer units | Gewirtzman et al., in review | 01 → 02, 03 | CC BY 4.0 |
| `from_ymf_soil_zenodo/` | Yale-Myers Forest upland–wetland soil CH₄ fluxes and soil closure windows | Gewirtzman, J. (2026), Zenodo, doi:10.5281/zenodo.18779715 (v0.1.0) | 01, 20 | CC BY 4.0 (as published) |
| `from_agsoil_study/` | De-identified agricultural-soil CH₄ and N₂O closures | J. Gewirtzman et al., in prep. | 20 | CC BY 4.0 |

`external/` (not in the repository) is an optional local folder for full source files, such as
the complete soil closure windows of the Zenodo record (`SOIL_WINDOWS`, see
`from_ymf_soil_zenodo/README.md`).
