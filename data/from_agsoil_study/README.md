# Copy of: agricultural-soil CH₄ and N₂O closures (de-identified)

**Source.** J. Gewirtzman et al., in preparation (a manure-management study in a managed
agricultural grassland). The full data, with their site and experimental context, will be
released with that paper; cite it for any use.

**Not new data.** This file is a minimal, de-identified copy of the source study's closure-level
results, reproduced here only so that this paper's analysis (its agricultural-soil rows of Table S4
and the cross-dataset figures) can be rerun. It is not a new or independent dataset: obtain the data from the source study, and do
not count these closures as additional or independent data in compilations, syntheses or
meta-analyses. The de-identification also makes this copy unsuitable for such use (no site,
treatment, plot or date).

## What is included, and why

One growing season of soil closures, CO₂/CH₄ on a LI-COR LI-7810 and N₂O on a LI-COR LI-7820 in
the same closed loop (1 Hz). Fluxes, model choice, per-closure-group precision and MDF were computed
by the source study with goFlux 0.5.0.9002 and its additions, with the same detection convention
as this paper (z = 1.96, closure-median second-difference precision, t in seconds).

`agsoil_closures_deidentified.csv`: all 535 closures of the de-identified table, rows in shuffled
order; 9 of its 18 columns, the ones `20_generality.R` reads (values verbatim):

| Column | Unit | Description |
|---|---|---|
| `CO2_flux` | µmol m⁻² s⁻¹ | CO₂ flux (used only to recognise duplicate rows) |
| `CH4_flux`, `N2O_flux` | nmol m⁻² s⁻¹ | Best-fit flux |
| `CH4_r2`, `N2O_r2` | | R² of the best-fit model |
| `CH4_sigma`, `N2O_sigma` | ppb | Campaign precision (closure median of the second-difference σ) |
| `CH4_MDF`, `N2O_MDF` | nmol m⁻² s⁻¹ | Reference MDF, 1.96 σ / t × flux term |

De-identified for this paper: site, plots, collars, experimental treatments and dates are not
included; closure IDs, month, closure length, model names, standard errors and QC flags of the
earlier de-identified table are also removed here because the scripts do not use them.

**Licence.** CC BY 4.0 (`LICENSE-data`).
