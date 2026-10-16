# Closure-level results (stems)

`stem_closures_processed.csv`: one row per stem closure analysed in the paper (1547), written by
`scripts/23_export_closures.R` from the pipeline outputs. Fluxes in nmol m⁻² s⁻¹ (CH₄) and
µmol m⁻² s⁻¹ (CO₂); concentrations in ppb.

| Column | Meaning |
|---|---|
| `closure_id` | Closure identifier used throughout the pipeline |
| `date`, `site`, `plot`, `tree`, `species` | Sampling date; upland or wetland site; plot code; tree tag; species code |
| `analyzer` | Noisier (LGR GLA131, 2023–2024) or precise (LI-COR LI-7810, 2025) |
| `air_temp_C` | Air temperature at the time of the closure (°C) |
| `n_readings`, `logging_interval_s`, `closure_s` | Readings in the fitted window; logging interval (s); closure length *t* (window span plus one interval, s) |
| `flux_term` | Converts a slope in ppb s⁻¹ to a flux in nmol m⁻² s⁻¹ (V P (1 − x_H2O) / (A R T)) |
| `CH4_flux`, `CH4_model`, `CH4_r2` | Best-fitting CH₄ flux, its model (LM linear, HM Hutchinson–Mosier) and R² |
| `CO2_flux` | Concurrent CO₂ flux |
| `sigma_closure_ppb` | The closure's own precision, 1.4826 MAD of second differences / √6 |
| `sigma_campaign_ppb` | Campaign precision: median of `sigma_closure_ppb` over the analyzer's closures on that day |
| `MDF_reference` | Reference MDF: 1.96 × campaign σ / *t* × flux term |
| `MDF_datasheet` | MDF from the manufacturer's 1 s precision with multiplier 1 (the goFlux default) |
| `detection_class` | *emission* (flux > MDF), *uptake* (flux < −MDF) or *below MDF*, against `MDF_reference` |
| `reg_half_95` | Half-width of the 95 % regression interval of the slope (Cowan et al. 2025) with the campaign σ, in flux units |
| `lm_p` | p-value of the straight-line slope fitted to the raw readings |
| `resid_ac1` | Lag-1 autocorrelation of that line's residuals (near 0 for independent noise) |
| `qc_*` | Closure screens (flags for review; nothing was removed): `c0` start far from ambient, `co2_tracer` CO₂ not rising, `convex` accelerating trace, `min_obs` too few readings, `ambient` headspace at sealing far from the air before, `noisy` noisier than the campaign; `qc_any` any of these |
