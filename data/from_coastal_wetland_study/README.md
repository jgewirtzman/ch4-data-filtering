# Coastal wetland closures (copied from a source project; not new data)

Chamber closures of CO2, CH4 and N2O from a mangrove-saltmarsh tidal creek on Woogoompah
Island, southern Moreton Bay, Australia, winter campaign of 15-19 August 2024, measured by
J. A. Rosentreter and colleagues (Centre for Coastal Biogeochemistry, Southern Cross
University). Other results of the project: Vienna-Hallam et al. (2026),
doi:10.1029/2026JG009795; van Haastregt et al. (2026), doi:10.1029/2026JG009861.

These files are a minimal copy, reduced to what this analysis uses. They are not a release
of the source data: do not compile them into syntheses; contact the data owner instead.

| File | Content |
|---|---|
| `coastal_windows_clicked.csv` | Fitting window of each closure, picked by eye from the raw record (`scripts/24a_click_coastal_windows.R`) |
| `coastal_closures.csv` | One row per closure and gas: selected goFlux model and flux, linear R2 and p-value, per-closure and campaign second-difference sigma, closure length, reference and datasheet MDF, regression interval, residual and second-difference lag-1 autocorrelation, habitat, light/dark |

Units: CO2 umol m-2 s-1; CH4 and N2O nmol m-2 s-1; sigma in ppm (CO2) or ppb (CH4, N2O).
The raw analyzer records (LI-7810 and LI-7820, 1 Hz) are held by the data owner and are not
included.
