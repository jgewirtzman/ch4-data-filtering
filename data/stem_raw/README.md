# Raw 1 Hz analyzer records (not in this repository)

The raw concentration records for the stem closures (~650 MB) are archived separately,
in the Harvard Forest Data Archive package of the stem study (DOI to follow). Place them here to recompute precision and fluxes from
scratch with `02_allan_deviation.R` and `02b_goflux_fluxes.R`; everything downstream
(03 onward) runs from the intermediate results in `outputs/`.

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
Some closures in `data/stem/HF_2023-2025_tree_flux_corrected.csv` have no surviving raw
record (the 1 May 2024 files are empty, and no record exists for 1-2 Oct 2025). They
cannot be refitted or given their own precision, so `03_mdf_computation.R` drops them;
1547 stem closures are analysed.
