# Raw 1 Hz analyzer records (not in this repository)

The raw concentration records for the stem closures (~650 MB) belong to the Harvard Forest stem
study (J. Gewirtzman et al., in prep.) and will be archived in its Harvard Forest Data Archive
package (DOI to follow). This folder is only the documented local location for them: git ignores
its contents. Place the records here to recompute precision and fluxes from scratch with
`02_allan_deviation.R` and `02b_goflux_fluxes.R` (`scripts/run_all.R` then runs them and refreshes
`data/derived_from_stem_raw/`). Without them, everything from 03 onward runs from the saved
results of 02 and 02b in `data/derived_from_stem_raw/`.

The readings of the few closures drawn in Figs. 2, 5 and S1 are copied from these
records to `data/from_hf_stem_study/figure_closures/` (tracked) by `14_figure_closures.R`, so the
figure scripts 15 and 17 never read this folder.

## 7810/
LI-COR LI-7810 records from the 2025 Harvard Forest stem campaigns.

- `raw-7810-tree_flux_data/TG10-01861-*.data`: the analyzer's continuous daily records
  (1 Hz), 16 field days, 9 Apr - 25 Aug 2025. Each closure is a run of the REMARK
  column (the tree tag typed at closure).

## LGR/
ABB/Los Gatos GLA131 (units LGR1, LGR2, LGR3) continuous daily records from the
2023-2024 Harvard Forest stem campaigns, in date-named folders. Some date folders hold
several days as sub-folders (e.g. `2024-12-03/2024-12-02`); the catalog in 02 descends
into them.

## Closures without a raw record
Some closures in `data/from_hf_stem_study/HF_2023-2025_tree_flux_corrected.csv` have no surviving raw
record (the 1 May 2024 files are empty, and no record exists for 1-2 Oct 2025). They
cannot be refitted or given their own precision, so `03_mdf_computation.R` drops them;
1547 stem closures are analysed.
