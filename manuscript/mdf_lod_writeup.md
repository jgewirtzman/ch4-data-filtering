# Minimal Detectable Flux and Quality Filtering: Methods and Results

## 1. Methods

### 1.1 Instrument precision

All flux measurements were made with ABB/Los Gatos Research GLA131-series Microportable Greenhouse Gas Analyzers (UGGA). The manufacturer specifies CH4 precision of 0.9 ppb and CO2 precision of 0.35 ppm (1-sigma at 1 s averaging). Because field conditions (vibration, temperature fluctuations, aging optics) can degrade precision relative to laboratory specifications, we estimated empirical precision using two complementary approaches.

**Allan deviation (per-measurement).** For each flux measurement, we estimated Allan deviation at tau = 1 s from the dry-mole-fraction timeseries within the measurement window. Because concentrations change monotonically during a chamber closure, first-differencing is used to remove the flux-driven trend before estimating noise, yielding a measurement-specific precision estimate that is independent of the flux model fit.

**Rolling window scan (per-instrument).** We scanned the full continuous timeseries for each instrument using a 30-second rolling window. For each window we computed the absolute slope (via linear regression) and the residual standard deviation around the fit. Windows were ranked by a combined normalized score of |slope| + residual SD; the bottom 5% (flattest and quietest) were selected as representative of true instrument noise floors. This yields a single global precision estimate per instrument.

### 1.2 Minimal Detectable Flux

We computed the Minimal Detectable Flux (MDF) using three published approaches, all sharing the general form MDF = precision_term / t x flux.term, where t is the measurement duration in seconds and flux.term encodes the chamber volume-to-area ratio and temperature via the ideal gas law (computed per measurement by goFlux):

1. **Manufacturer MDF (goFlux).** Uses the GLA131 datasheet precision (0.9 ppb CH4): MDF = precision_manufacturer / t x flux.term.

2. **Wassmann et al. (2018).** Uses the global empirical precision from the rolling window method, scaled by a factor of 3 for the desired confidence level: MDF = (3 x SD_global) / t x flux.term. Computed at 99%, 95%, and 90% confidence (using z-scores of 2.576, 1.960, and 1.645 respectively in place of the factor 3).

3. **Christiansen et al. (2015) MQL variant.** Uses the per-measurement Allan deviation, the factor 3, and the t-distribution critical value: MDF = (SD_allan x 3 x t_alpha) / t x flux.term, where t_alpha is the two-tailed t critical value at the specified confidence level with df = t - 2. Computed at 99%, 95%, and 90% confidence.

### 1.3 Additional quality criteria

Beyond MDF, we evaluated several additional quality filters commonly applied to chamber flux measurements:

- **CH4 R-squared thresholds** (> 0.7, > 0.9): The coefficient of determination from the best-fit model (linear or Hutchinson-Mosier, as selected by goFlux).
- **CO2 R-squared threshold** (> 0.7): Applied to the concurrent CO2 flux from the same measurement as an independent check on measurement quality.
- **Signal-to-noise ratio (SNR)** thresholds (> 2, > 3): Defined as |CH4 flux| / noise_floor, where noise_floor = σ_Allan / t × flux.term converts the per-measurement Allan deviation (1-sigma, concentration units) to flux units via the measurement duration and chamber geometry. This empirical, per-measurement definition captures instrument-specific and condition-specific noise rather than relying on manufacturer specifications.

### 1.4 Sensitivity analysis

To assess how the choice of quality filter propagates through whole-tree flux scaling, we applied each of the 13 criteria independently to the Harvard Forest dataset (n = 136 measurements). For each filter, measurements failing the criterion had their CH4 flux set to zero; passing measurements retained their original flux values. We then computed mean per-tissue flux rates for four compartments matching the truncation analysis (stem < 2 m, stem >= 2 m, branches, leaves) and scaled these by Whittaker & Woodwell (1967) surface area indices (stem bark = 0.45, branch bark = 1.70, LAI = 4.5 m2 m-2 ground) to obtain integrated stand-level budgets. The cone taper model (26 m tree, 40 cm DBH) was used to partition stem bark area into below- and above-2 m fractions.

The analysis was repeated with outliers removed (measurements exceeding Q3 + 3 x IQR per compartment) to assess the combined effect of quality filtering and outlier sensitivity.

---

## 2. Figures

### Figure 1. Empirical vs. manufacturer CH4 precision

![precision_comparison_CH4](plots/precision/precision_comparison_CH4.png)

**Caption.** Comparison of CH4 analyzer precision estimates across three ABB GLA131 instruments deployed at Harvard Forest (LGR1, LGR2, LGR3) and one at Yale-Myers Forest (LGR3-YMF). Red bar: manufacturer specification (0.9 ppb, 1-sigma at 1 s). Light blue bars: empirical precision from the rolling window method (residual SD of the quietest 5% of 30-second windows in the continuous timeseries). Dark blue bars: median Allan deviation computed within individual flux measurement windows. Empirical precision generally exceeds (i.e., is worse than) the manufacturer specification, with substantial inter-instrument variability. LGR1 shows the largest discrepancy (Allan deviation 3.6 ppb vs. 0.9 ppb specification), while LGR2 CH4 precision (0.6--0.7 ppb) is actually better than the datasheet value.

---

### Figure 2. CH4 flux distributions by tissue component under quality filtering

![component_flux_by_mdf_filter](plots/precision/component_flux_by_mdf_filter.png)

**Caption.** Distribution of CH4 flux rates (nmol m-2 s-1) at Harvard Forest by tissue component (Stem, Branch, Leaf), shown for each of 13 quality filter criteria applied independently. Rows are ordered from most permissive (No filter, 136/136 measurements retained) to most restrictive (CH4 R-squared > 0.9, 51/136 retained). Density ridgelines show the flux distribution using all original flux values; the fill color indicates the percentage of measurements retained under each filter for that component (blue = high retention, red = low retention). Jittered points below each ridge show individual measurements: grey points are retained at their original flux value, red points have been set to zero by the filter. Vertical red lines mark the component mean flux (computed with zeroed values included). The x-axis uses an asinh transformation to accommodate both near-zero and large positive fluxes. Y-axis labels report the number and percentage of measurements passing each filter across all components. Stem fluxes are robust across filter definitions; branch and leaf means are more sensitive, particularly to Christiansen MDF thresholds and R-squared filters.

---

### Figure 3. Integrated stand-level CH4 budget under quality filtering

![budget_stacked_by_filter](plots/precision/budget_stacked_by_filter.png)

**Caption.** Sensitivity of the integrated stand-level CH4 budget to quality filter choice, with outliers retained (top row) and removed (bottom row, 3 x IQR criterion). Each stacked bar shows the integrated flux (nmol m-2 ground s-1) partitioned among four tissue compartments: stem < 2 m (dark brown), stem > 2 m (tan), branches (blue), and leaves (green). Integrated flux for each compartment equals the filtered mean per-tissue flux rate multiplied by the corresponding Whittaker & Woodwell (1967) surface area index. Bars are ordered from most permissive (No filter, 100%) to most restrictive (Christiansen 99%, 38%) with the percentage of retained measurements in parentheses. Annotated values above each bar give the total integrated flux. With outliers retained, the total budget ranges from 0.692 (no filter) to 0.582 nmol m-2 ground s-1 (Christiansen 99%), a maximum reduction of 16%. With outliers removed, the budget ranges from 0.539 to 0.425, a maximum reduction of 21%, reflecting the additional loss of high-flux stem measurements. The shared y-axis (0--0.8) facilitates direct comparison between rows.

---

### Figure 4. Relative sensitivity of component flux rates to quality filtering

![flux_rate_uncertainty_by_component](plots/precision/flux_rate_uncertainty_by_component.png)

**Caption.** Sensitivity of per-tissue CH4 flux rate estimates to quality filter definition, shown in absolute units (top) and as a percentage of the no-filter baseline (bottom). Each point represents the mean flux rate for one compartment under one of 13 quality filter criteria, with filled circles for outliers-retained and filled triangles for outliers-removed analyses. Large open symbols mark the median across all filter definitions. Vertical bars span the full range. In absolute terms, stem < 2 m fluxes dominate (approximately 1.0 nmol m-2 s-1) and are insensitive to filter choice (CV = 2.2%). Relative sensitivity increases with decreasing flux magnitude: stem > 2 m (CV = 4.8%, range 88--100%), branches (CV = 4.9%, range 82--100%), and leaves (CV = 9.0%, range 78--100% with outliers retained). With outliers removed, leaf sensitivity increases dramatically (0--100% range), reflecting dependence on a single high-flux measurement. Stem > 2 m sensitivity also roughly doubles (76--100% range) without outliers.

---

## 3. Results

### 3.1 Empirical precision versus manufacturer specification

Empirical CH4 precision varied substantially across instruments and estimation methods (Fig. 1). For LGR1, the median Allan deviation (3.6 ppb) was four times the manufacturer specification (0.9 ppb), and the rolling window estimate (2.5 ppb) was nearly three times the specification. LGR3 showed a more moderate discrepancy (Allan deviation 1.5 ppb, rolling window 1.2 ppb). LGR2 was the exception: both empirical estimates (Allan deviation 0.6 ppb, rolling window 0.7 ppb) were better than the datasheet value. The Yale-Myers instrument (LGR3 redeployed) yielded empirical precision of 1.2--1.3 ppb, consistent with its Harvard Forest performance. These results confirm that manufacturer specifications can substantially underestimate real-world noise for some instruments, and that relying on a single datasheet value across instruments introduces uncontrolled variability into detection limits.

CO2 precision showed a similar pattern: empirical estimates ranged from 0.5 to 1.8 ppm across instruments, compared to the manufacturer specification of 0.35 ppm.

### 3.2 Minimal Detectable Flux

The three MDF approaches produced detection thresholds spanning roughly an order of magnitude for CH4 (Table in quality criteria summary). The manufacturer-based MDF was the most permissive, flagging 20% of CH4 measurements as below detection. Wassmann thresholds (which substitute empirical global precision for the datasheet value) flagged 26--34%, depending on confidence level. The Christiansen approach — which combines per-measurement Allan deviation with both a factor of 3 and a t-distribution critical value — was the most conservative, flagging 52--62% as below detection. Nearly all CO2 fluxes (>95%) exceeded detection thresholds under every method.

### 3.3 Quality filter sensitivity: component flux rates

Applying 13 quality filters independently to the Harvard Forest dataset (n = 136) revealed that stem flux rates (both below and above 2 m) were robust to filter definition, while branch and leaf rates were more sensitive (Fig. 2, Fig. 4).

**Stem < 2 m** fluxes — the largest component (~1.0 nmol m-2 s-1) — varied by only 7% across all filter definitions (range 93--100% of the no-filter baseline, CV = 2.2%). These fluxes are large relative to any detection threshold and are therefore insensitive to the choice of quality criterion.

**Stem > 2 m** fluxes (range 88--100%, CV = 4.8%) and **branch** fluxes (range 82--100%, CV = 5.0%) showed moderate sensitivity. The most restrictive filters (Christiansen 99%, CH4 R2 > 0.9) reduced mean branch flux by up to 18%, primarily by zeroing measurements near the detection limit.

**Leaf** fluxes were the most sensitive component (range 78--100%, CV = 9.0% with outliers retained). Because leaf CH4 fluxes are small and many individual measurements fall near or below detection, the choice of filter criterion has a proportionally larger effect on the component mean. With outliers removed, leaf sensitivity increased dramatically (0--100% range), reflecting dependence on a single high-flux measurement that, when retained, anchors the leaf mean above zero across all filter definitions.

### 3.4 Quality filter sensitivity: integrated stand-level budget

Despite the variable sensitivity of individual components, the integrated stand-level CH4 budget was remarkably stable across filter definitions (Fig. 3). With outliers retained, the total budget ranged from 0.692 nmol m-2 ground s-1 (no filter) to 0.582 nmol m-2 ground s-1 (Christiansen 99%), a maximum reduction of 16%. This stability reflects the dominance of stem < 2 m fluxes — which are insensitive to filtering — in the integrated budget after scaling by surface area indices.

With outliers removed, the budget ranged from 0.539 to 0.425 nmol m-2 ground s-1 (maximum reduction 21%). The larger relative sensitivity arises for two reasons: first, most outliers are high-flux stem > 2 m measurements, so removing them makes that component substantially more sensitive to subsequent filtering (CV roughly doubles from 4.8% to 9.9%); second, outlier removal shrinks the total budget baseline, so the same absolute filtering effect represents a larger percentage reduction.

The empirical SNR filters (SNR > 2 and SNR > 3) produced intermediate results: 96% and 95% of the unfiltered budget, respectively, placing them between the Wassmann and Christiansen thresholds in stringency.

---

## 4. Discussion

### 4.1 Manufacturer precision is not sufficient

The four-fold discrepancy between LGR1's empirical noise and the manufacturer specification demonstrates that datasheet precision cannot be assumed to hold under field conditions. Vibration during canopy access, temperature fluctuations in the forest environment, and optical degradation over multi-season deployments all plausibly contribute. The inter-instrument variability (LGR2 outperforming the spec while LGR1 substantially underperformed) further argues against using a single manufacturer value. Any MDF calculation that relies on the datasheet value — including the default goFlux implementation — will systematically underestimate detection limits for noisier instruments.

### 4.2 The choice of quality filter matters most for small fluxes

The sensitivity analysis reveals a clear hierarchy: filter choice matters in proportion to how close a component's flux is to the detection limit. Stem < 2 m fluxes, which are 5--10x larger than any MDF threshold, are effectively immune to the choice of quality criterion. Leaf fluxes, by contrast, are comparable in magnitude to the detection limit and are therefore highly sensitive — the most sensitive of any component, with mean flux ranging from 0 to 100% of the unfiltered baseline when outliers are removed. This extreme sensitivity reflects both the small number of leaf measurements (n = 21) and the presence of a single outlier that anchors the component mean; its inclusion or exclusion alone spans the full range of possible outcomes.

Critically, even fluxes that fall below detection limits may scale to a substantive fraction of the total stand-level budget when multiplied by the large leaf area index (LAI = 4.5 m2 m-2 ground — the largest surface area index of any compartment). A leaf flux that is individually undetectable could, integrated over the full canopy, contribute meaningfully to whole-tree emissions. This means that the leaf component simultaneously contributes the highest methodological uncertainty and potentially the largest unresolved flux term in the budget.

Beyond detection limits, leaf CH4 fluxes — if present — are likely sensitive to environmental drivers that our measurements did not control, including light, vapor pressure deficit, and ambient and leaf-internal gas concentrations. Our measurements, like all components in this study, also cover only a narrow range of species, conditions, and times of year. We therefore do not attempt to generalize about the magnitude or role of foliar fluxes in whole-tree budgets; rather, the analysis highlights the compounding uncertainties in their measurement — instrumental (near-detection fluxes), environmental (uncontrolled drivers), and sampling (few measurements, limited species and conditions). This is not a limitation unique to our study; it reflects a broader gap across the field, where chamber-based approaches operating near their detection limits cannot definitively distinguish small biological fluxes from instrument noise in high-surface-area, low-flux compartments.

### 4.3 Budget stability despite component sensitivity

The integrated stand-level budget varies by at most 16% across 13 filter definitions (with outliers retained), despite individual component sensitivities reaching 22%. This stability arises because the budget is dominated by compartments with robust flux estimates. Branches (area index 1.70 m2 m-2) and stems (0.45 m2 m-2) carry moderate-to-large per-area fluxes that are well above detection limits and therefore insensitive to filtering. Leaves have the largest surface area index (LAI = 4.5 m2 m-2) but their per-area flux rates are so small — often near or below detection — that the integrated leaf contribution remains a minor fraction of the total budget despite the large canopy area.

This stability suggests that whole-tree CH4 emission estimates are not critically dependent on the specific quality criterion adopted — a reassuring result for cross-study comparisons where different groups may apply different detection thresholds. However, the 22% budget reduction when outliers are removed (compared to 16% from filtering alone) indicates that individual high-flux measurements exert a stronger influence on the budget than the choice of detection methodology.

### 4.4 Empirical SNR as an integrative quality metric

The empirical SNR — defined as the absolute flux divided by the per-measurement Allan deviation noise floor in flux units — offers several advantages over alternative quality metrics. Unlike the manufacturer-based MDF, it reflects the actual instrument noise during each specific measurement. Unlike R-squared thresholds, it is independent of the flux model (linear vs. Hutchinson-Mosier) and does not penalize measurements with genuinely small but real fluxes. And unlike the Christiansen MDF, which scales the Allan deviation by both a factor of 3 and a t-distribution critical value (yielding effective multipliers of ~5--8x), the SNR provides a direct, interpretable ratio of signal to noise.

In practice, SNR > 2 and SNR > 3 filtered 68% and 63% of measurements respectively, placing them in the middle of the filter stringency spectrum — more conservative than the manufacturer MDF (80% retained) but substantially more permissive than Christiansen 99% (38% retained). For datasets where preserving statistical power is important but manufacturer specifications are unreliable, empirical SNR thresholds offer a principled compromise.

---

## 5. Multi-instrument analysis (Harvard Forest stem flux dataset)

The following sections extend the analysis above using the Harvard Forest stem flux monitoring dataset (n = 1640 measurements: 1142 GLA131/UGGA measurements from 2023--2024, 498 LI-COR LI-7810 measurements from 2025), which provides sufficient sample size and instrument contrast to evaluate how below-MDF treatment and instrument precision interact with quality filtering.

### 5.1 Additional methods

**Instruments.** The LI-COR LI-7810 was deployed for stem flux measurements beginning in spring 2025. Its manufacturer-specified CH4 precision (0.6 ppb, 1-sigma at 1 s) is comparable to the GLA131 specification (0.9 ppb). Allan deviation was computed per-measurement from 1-Hz raw timeseries for both instruments: pre-segmented CSVs for the LI-7810 (60-second measurement windows), and continuous recordings segmented by field log timestamps for the GLA131 (variable duration, typically 180--300 seconds).

**Chamber geometry.** All stem measurements used permanent collars with radius 5.08 cm (surface area = 0.00811 m2). Chamber volume was computed from measured tree diameters and collar depth for each tree (tree_volumes.csv), with an additional 0.028 L for connecting tubing. The flux term (nmol / surface area, encoding P, V, and T via the ideal gas law) was computed per measurement.

**Instrument-specific global precision for Wassmann MDF.** Because the two instruments have very different empirical noise characteristics, the Wassmann global SD was computed separately per instrument (median Allan deviation: 3.96 ppb for GLA131, 0.095 ppb for LI-7810) rather than pooling across the full dataset.

**Below-MDF treatment comparison.** For each of the seven MDF thresholds (Manufacturer, Wassmann 90/95/99%, Christiansen 90/95/99%), we applied three treatments to measurements flagged as below detection:

1. **Remove:** exclude below-MDF measurements from the dataset entirely.
2. **Set to zero:** replace below-MDF flux values with zero, preserving sample size.
3. **Keep original (flag only):** retain all measured flux values; flag each measurement as above or below its per-measurement MDF.

We compared the resulting flux distributions, means, and medians against the unfiltered baseline, both pooled and split by instrument.

### 5.2 Additional figures

#### Figure 5. Per-measurement Allan deviation across instruments

![ch4_allan_deviation_by_instrument](plots/precision/ch4_allan_deviation_by_instrument.png)

**Caption.** Per-measurement CH4 Allan deviation (ppb, log scale) for GLA131/UGGA measurements (2023--2024 stem flux, n = 1206) and LI-7810 measurements (2025 stem flux, n = 378). Violins show the full distribution with quartile lines; individual measurements are jittered within each violin. Dashed horizontal lines mark manufacturer specifications (GLA131: 0.9 ppb, orange; LI-7810: 0.6 ppb, blue). The two instruments occupy non-overlapping ranges on a log scale. The GLA131 median Allan deviation (4.0 ppb) is 4.4x worse than its manufacturer specification, while the LI-7810 median (0.095 ppb) is 6.3x better than its specification. This approximately 40-fold difference in real-world CH4 precision dominates all downstream detection limit calculations and quality filtering outcomes.

---

#### Figure 6. Positive bias from removing below-MDF measurements

![ch4_mdf_treatment_bias](plots/precision/ch4_mdf_treatment_bias.png)

**Caption.** Mean (top) and median (bottom) CH4 flux under three below-MDF treatments across seven detection thresholds (n = 1640 stem flux measurements, 2023--2025). Grey: unfiltered baseline (all measurements retained at original values). Blue: below-MDF measurements set to zero. Red: below-MDF measurements removed entirely. The unfiltered mean (1.12 nmol m-2 s-1) and median (0.065 nmol m-2 s-1) serve as the reference. Setting below-MDF values to zero changes the mean by less than 5% across all thresholds, because symmetric noise around zero cancels in the mean. However, it pulls the median toward zero, reaching exactly zero for Christiansen thresholds (where >50% of measurements are flagged). Removing below-MDF measurements inflates the mean by 17% (Manufacturer MDF) to 135% (Christiansen 99%), and the median by 40% to 319%, because removal selectively eliminates near-zero measurements from both sides of the distribution, retaining only larger positive fluxes.

---

#### Figure 7. Bias quantification: % change from below-MDF treatment

![ch4_mdf_treatment_pct_change](plots/precision/ch4_mdf_treatment_pct_change.png)

**Caption.** Percent change in mean (top) and median (bottom) CH4 flux relative to the unfiltered baseline, for the two active treatments (set to zero, blue; remove, red) across seven detection thresholds. Removal consistently inflates both statistics, with the bias scaling monotonically with filter stringency. Setting to zero introduces minimal mean bias (< 5% at all thresholds) but suppresses the median by up to 100% at stringent thresholds where the majority of measurements fall below detection. Neither treatment is unbiased for both the mean and the median simultaneously; retaining original values (the unfiltered baseline) is the only approach that preserves both.

---

#### Figure 8. CH4 flux vs. per-measurement noise floor

![ch4_flux_vs_noise_floor](plots/precision/ch4_flux_vs_noise_floor.png)

**Caption.** CH4 flux (nmol m-2 s-1, asinh-scaled y-axis) plotted against the per-measurement noise floor (nmol m-2 s-1, log-scaled x-axis) for all measurements with Allan deviation estimates (n = 1416). The noise floor is defined as sigma_Allan / t x flux.term, converting the 1-sigma Allan deviation in concentration units to flux units for each measurement. Orange: GLA131/UGGA (2023--2024); blue: LI-7810 (2025). Dotted diagonal lines mark flux = +/- noise floor. Points within the dotted wedge are indistinguishable from zero at the 1-sigma level. The two instruments occupy distinct regions: GLA131 noise floors cluster around 0.02--1 nmol m-2 s-1, while LI-7810 noise floors are approximately 0.002--0.01 nmol m-2 s-1 --- an order of magnitude lower. Many GLA131 measurements, including the majority of negative fluxes, fall within the noise floor wedge. The LI-7810, with its lower noise floor, resolves small positive fluxes that would be indistinguishable from zero on the GLA131.

---

#### Figure 9. Below-MDF treatment effect by instrument

![ch4_mdf_treatment_by_instrument](plots/precision/ch4_mdf_treatment_by_instrument.png)

**Caption.** Effect of below-MDF treatment on CH4 flux distributions, shown separately for GLA131/UGGA (left column) and LI-7810 (right column) across seven detection thresholds (rows). Grey: unfiltered; blue: set below-MDF to zero; red: remove below-MDF. Removal causes a progressive rightward shift and flattening of the GLA131 distribution as successively more near-zero measurements are excluded --- by Christiansen 99%, the remaining distribution bears little resemblance to the original. The set-to-zero treatment instead develops a growing spike at zero while preserving the right tail shape. The LI-7810 distributions are nearly indistinguishable across all three treatments at every threshold, because the instrument's high precision means very few measurements fall below even the most stringent MDF. This demonstrates that the bias from removal is predominantly an artifact of instrument noise, not a property of the underlying flux distribution.

---

### 5.3 Additional results

#### 5.3.1 Multi-instrument empirical precision

Empirical CH4 precision differed by approximately 40-fold between the two instrument types (Fig. 5). The GLA131 median Allan deviation across 1206 stem flux measurements was 4.0 ppb (IQR: 3.2--4.9 ppb), consistent with the 3.6 ppb observed for LGR1 in the canopy dataset and 4.4x worse than the manufacturer specification of 0.9 ppb. The LI-7810 median was 0.095 ppb (IQR: 0.086--0.108 ppb), 6.3x better than its manufacturer specification of 0.6 ppb, with a remarkably tight distribution (IQR spanning less than a factor of 1.3). The two distributions were completely non-overlapping on a log scale.

This gap is far larger than the 1.5-fold difference implied by manufacturer specifications (0.9 vs. 0.6 ppb). The LI-7810 also outperformed the GLA131 for CO2 precision (0.32 vs. 1.37 ppm median Allan deviation), though both were within an order of magnitude and both outperformed the LI-7810 CO2 specification of 3.5 ppm.

#### 5.3.2 Instrument-dependent MDF outcomes

The instrument precision gap propagated directly into detection limit calculations. The instrument-specific Wassmann global SD (GLA131: 3.96 ppb, LI-7810: 0.095 ppb) produced thresholds roughly 40x higher for GLA131 measurements. Consequently, Wassmann 99% flagged 46% (531/1142) of GLA131 measurements but only 15% (74/498) of LI-7810 measurements as below detection.

#### 5.3.3 Negative fluxes: mostly instrument noise

The stem flux dataset contained 284 negative CH4 flux measurements (17.3% of 1640), distributed asymmetrically between instruments: 275/1142 (24.1%) for the GLA131 vs. 9/498 (1.8%) for the LI-7810. Three lines of evidence indicate that the majority of GLA131 negative fluxes reflect instrument noise rather than biological CH4 uptake:

1. **Apples-to-apples seasonal comparison.** The GLA131 covered year-round measurements (Jun 2023--Dec 2024) while the LI-7810 covered the growing season only (Apr--Oct 2025). Restricting the comparison to overlapping calendar months (Jun--Oct), the GLA131 still showed 20.3% negative fluxes (169/832) vs. 1.7% (6/353) for the LI-7810. Same trees, same season, same chambers --- the difference is predominantly instrument noise.

2. **Noise floor analysis.** Of 246 GLA131 negative fluxes with Allan deviation estimates, 52% had absolute values within 1-sigma of the per-measurement noise floor, 74% within 2-sigma, and 84% within 3-sigma. The median absolute negative flux (0.049 nmol m-2 s-1) was comparable to the median noise floor (0.051 nmol m-2 s-1). These values are statistically indistinguishable from zero (Fig. 8).

3. **Distribution shape.** On the GLA131, the near-zero flux distribution was roughly symmetric around a small positive value, consistent with Gaussian noise broadening a weakly positive signal. On the LI-7810, the same trees resolved as a tight, positively-skewed distribution with almost no negative tail.

The GLA131 negative flux rate did show a seasonal pattern --- higher in dormant months (35.7% Nov--Apr) than growing months (20.5% May--Oct) --- consistent with real fluxes approaching zero during dormancy. But even in the growing season, 20% negative fluxes from the GLA131 vs. <2% from the LI-7810 for the same months indicates that the dominant source of negative values is noise scatter, not biology. The seasonal pattern likely reflects real fluxes approaching zero (bringing more of the distribution within the noise wedge), rather than true biological CH4 uptake.

#### 5.3.4 Below-MDF treatment comparison

**Removal inflates both mean and median.** Removing below-MDF measurements shifted the mean flux upward by 17% (Manufacturer MDF) to 135% (Christiansen 99%) relative to the unfiltered baseline (Fig. 6, Fig. 7). The median was even more severely affected: +40% to +319%. This positive bias arises because removal selectively eliminates near-zero measurements from both sides of the distribution --- both small positive and small negative fluxes --- leaving a distribution enriched in larger positive values. The bias scales monotonically with filter stringency.

**Setting to zero preserves the mean but distorts the median.** Replacing below-MDF measurements with zero changed the mean by less than 5% at all thresholds, because symmetric noise scattering around zero cancels in the arithmetic mean. However, the median dropped toward zero and reached exactly zero for all Christiansen thresholds, where more than 50% of measurements were flagged as below detection. This reflects an artificial pile-up of values at exactly zero that was not present in the original distribution.

**Keeping original values preserves both.** Retaining all measurements at their original values --- flagging but not modifying those below detection --- is the only treatment that preserves both the mean and median of the flux distribution.

**The bias is instrument-dependent (Fig. 9).** The treatment effect was almost entirely confined to the GLA131. For the LI-7810, all three treatments produced nearly indistinguishable distributions at every MDF threshold. For the GLA131, removal at Christiansen 99% inflated the instrument-specific mean from 1.44 to 4.61 nmol m-2 s-1 (a 220% increase), while the LI-7810 mean changed by only 51% (0.38 to 0.57) under the same threshold. This confirms that the removal bias is driven by instrument noise, not by the biological flux distribution.

#### 5.3.5 R-squared filters preferentially retain high-flux measurements

CH4 R2 > 0.5 and Christiansen 99% each retained the same number of measurements (665/1640, 40.5%), but with very different compositions. R2 > 0.5 retained only 1.4% negative fluxes, while Christiansen 99% allowed 5.6% negatives. R-squared selects measurements where the signal is clear above the noise; MDF selects measurements where the precision is sufficient relative to the measurement geometry and duration. A tree with a true flux of 0.02 nmol m-2 s-1 measured on a GLA131 with 4 ppb noise will produce a concentration change of ~0.3 ppb over a 5-minute measurement --- far below the noise floor, and thus never yielding a high R-squared, even though it is a valid measurement of a small flux.

---

## 6. Synthesis: question-dependent filtering recommendations

The results above demonstrate that no single filtering approach is universally optimal. The appropriate treatment of below-detection measurements depends fundamentally on the scientific question being asked. We identify two broad categories of questions that call for different filtering strategies.

### 6.1 Scaling and budgeting questions: retain all fluxes

**Question type:** "What is the mean flux or flux distribution for this ecosystem component?" Examples include computing stand-level CH4 budgets, estimating seasonal emission totals, comparing mean flux rates across species or sites, and scaling chamber measurements to landscape-level estimates.

**Recommended approach: keep original values, flag detection status.**

For questions about means, totals, and distributions, removing below-MDF measurements introduces a systematic positive bias that can exceed 100% at stringent thresholds (Section 5.3.4). Setting below-MDF values to zero avoids this bias in the mean but distorts the median and creates an artificial mode at zero. Retaining original values is the only treatment that preserves both the mean and median of the flux distribution.

The key insight is that instrument noise scatters symmetrically around the true flux value. For a tree emitting 0.01 nmol m-2 s-1, individual measurements on a noisy analyzer will scatter across a range that includes negative values --- but the mean of those measurements will converge on the true flux. Removing the below-detection scatter (which includes both the negative and the small-positive values) selectively retains the high-scatter tail, inflating the mean. This is the same logic underlying the standard treatment of non-detects in analytical chemistry (e.g., Helsel 2005), where substituting values for censored observations (including zero substitution) is recognized as biased, and maximum-likelihood or Kaplan-Meier approaches that use the full dataset are preferred.

The practical recommendation is:

- **Tier 1 (hard QC):** Remove measurements with clear evidence of equipment failure or procedural error --- chamber seal failures (negative CO2 flux from live tissue), equipment malfunctions noted in field logs, extreme Allan deviation outliers indicating instrument instability. These are quality flags based on what went wrong, not statistical thresholds on the flux value.
- **Tier 2 (detection flagging):** Compute the per-measurement MDF using the Christiansen method with empirical Allan deviation. Flag each measurement as above or below its individual detection limit. Report the fraction below detection. Retain all flux values.
- **Tier 3 (sensitivity analysis):** For key reported results, demonstrate robustness by showing results under all three treatments (unfiltered, set-to-zero, remove). If a conclusion changes depending on treatment, that is important information --- it indicates dependence on measurements that cannot be confidently distinguished from noise.

This approach is directly supported by the Harvard Forest canopy analysis (Section 3.4), which showed that the integrated stand-level budget was stable to within 16% across all filter definitions even when the set-to-zero approach was used --- and would be even more stable if original values were retained, because the set-to-zero treatment slightly underestimates the contribution of near-detection compartments like leaves.

### 6.2 Detection and biological signal questions: filter to the confident subset

**Question type:** "Is the flux from this surface significantly nonzero, and in what direction?" Examples include testing whether tree stems are net CH4 sources or sinks, whether leaves emit or consume CH4, whether a given species produces nonzero fluxes, and whether fluxes differ significantly from zero across seasons or treatments.

**Recommended approach: filter to above-MDF measurements only.**

For questions about whether a flux exists at all --- and especially about its sign --- the relevant subset is the measurements where the instrument can actually resolve the signal. Including below-detection measurements in a test of "is the flux nonzero?" adds noise without adding information, reducing statistical power. And including a population of noise-scattered near-zero values when asking "is this surface a source or a sink?" can obscure a real directional signal.

In this context, the positive bias from removal is not a bug --- it is the expected behavior. If you ask "among measurements where we confidently detected a flux, what is the distribution?", the answer should exclude the unresolvable measurements. The bias relative to the unfiltered distribution is the point: the unfiltered distribution includes noise that cannot inform the biological question.

The practical recommendation is:

- Apply the Christiansen MDF at the desired confidence level (99% for conservative detection claims, 90% for exploratory analyses).
- Report results for the above-MDF subset only, with the sample size and fraction above detection clearly stated.
- Interpret the resulting distribution as "the distribution of confidently detected fluxes," not as the population distribution. The two are different, and the distinction matters.

This is the appropriate approach for claims like: "We detected significant CH4 emission from X% of leaf measurements (Christiansen 99% MDF), with a median flux of Y among detected measurements." It is not appropriate for: "The mean leaf CH4 flux at this site is Y" --- that estimate should use the full dataset (Section 6.1).

### 6.3 R-squared as a detection metric, not a quality filter

R-squared thresholds deserve separate consideration because they interact with the two question types differently. R-squared measures how well a model fits the concentration timeseries, which is a function of both the flux magnitude and the instrument noise. A measurement with R2 > 0.7 is one where the concentration trend was clear above the noise --- a "detection" in the same spirit as the MDF threshold.

However, R-squared is a biased detection metric because it is not independent of flux magnitude. Two measurements with the same instrument noise but different true fluxes will have different R-squared values; the larger flux always produces a higher R2. This means R-squared filters preferentially retain high-flux measurements and systematically exclude low emitters --- exactly the opposite of what is needed for unbiased population estimates (Section 6.1), and a subtly different selection than MDF (which accounts for chamber geometry and duration, not just signal magnitude).

R-squared is most useful as a confidence flag for individual measurement interpretation: "this particular measurement showed a clear monotonic trend." It is least useful as a population filter for computing means or budgets.

### 6.4 Instrument precision as the dominant methodological variable

The approximately 40-fold difference in real-world CH4 precision between the GLA131 and LI-7810 (Fig. 5) overwhelmed every other source of variability in our quality filtering analysis. The LI-7810 flux distribution was nearly invariant across all MDF thresholds, all treatments, and all filter criteria. Every methodological question addressed here --- the choice of MDF approach, the treatment of below-detection values, the impact of R-squared thresholds --- was driven almost entirely by the GLA131 data.

This has two implications. First, for studies primarily concerned with small CH4 fluxes (e.g., dormant-season stem emissions, low-emitting species, leaf fluxes), the choice of analyzer may be more consequential than the choice of quality filter. A precise instrument renders the entire filtering question nearly moot, because almost all measurements exceed detection. Second, when combining data from instruments of differing precision, the noisier instrument contributes disproportionately to the fraction of below-detection measurements and thus to the sensitivity of results to filtering methodology. Multi-instrument datasets should present instrument-specific empirical precision and detection limit statistics alongside any pooled quality-filtered results.

### 6.5 Summary of recommendations

| Question | Goal | Below-MDF treatment | Appropriate filters |
|----------|------|---------------------|---------------------|
| Mean flux / budgets / scaling | Unbiased population estimate | Keep original values + flag | Hard QC only (seal failures, equipment errors) |
| Seasonal / spatial patterns in mean flux | Unbiased comparison of distributions | Keep original values + flag | Hard QC only; report fraction below MDF per group |
| Is this flux nonzero? | Detection confidence | Remove below MDF | Christiansen MDF at 90/95/99% |
| Is this a source or sink? | Sign determination | Remove below MDF | Christiansen MDF + consider sign test on above-MDF subset |
| Which trees are big emitters? | Identify high-flux individuals | Remove below MDF (or R2 filter) | R2 > 0.5 or Allan SNR > 3 appropriate here |
| Model fitting / driver analysis | Confident flux estimates for regression | Remove below MDF | Christiansen or Allan SNR; report selection and fraction retained |
| Cross-study comparison | Comparable methodology | Report all three treatments | Show sensitivity; flag instrument and empirical precision |

---

## 7. Data transformations for statistical analysis

CH4 flux data present several challenges for standard statistical methods: distributions are typically right-skewed with long tails, variances differ by orders of magnitude across ecosystem components, and --- critically --- fluxes can be negative (representing CH4 uptake/oxidation). The choice of data transformation interacts directly with the quality filtering considerations above, because transformations that cannot handle negative values effectively impose an additional, often unrecognized, data filter.

### 7.1 The problem with log transformations

The log transformation (log10 or ln) is the most common variance-stabilizing transformation in ecology and biogeochemistry. For strictly positive data, it compresses long right tails and stabilizes variance across groups with different means. However, the log of a non-positive number is undefined, and applying a log transform to CH4 flux data silently drops all zero and negative values.

For datasets where negative fluxes are rare or arise solely from instrument noise (e.g., canopy CH4 emissions measured with a precise analyzer), this data loss may be minor. But for datasets where negative fluxes represent real biological processes --- CH4 oxidation in upland soils, net uptake by methanotrophs, or consumption at aerobic surfaces --- log transformation discards a scientifically meaningful portion of the distribution.

We demonstrate this with three CH4 flux datasets spanning a range of negative-flux prevalence (Fig. 10):

- **Wetland soil fluxes** (Yale-Myers Forest, semirigid chambers, n = 288): 83% of measurements are negative, reflecting methanotrophic CH4 consumption dominating over methanogenesis at most plots. A log10 transform retains only 49 of 288 measurements (17%) --- the rare, large emission events --- and produces a distribution that bears no resemblance to the underlying flux population.

- **Tree stem fluxes** (Yale-Myers Forest, semirigid chambers, n = 369): 31% of measurements are negative. Log transformation drops 115 values, systematically removing the low-emission and uptake portion of the distribution and inflating the apparent mean.

- **Tree canopy fluxes** (Harvard Forest, n = 136): 10% negative. Log transformation drops 14 measurements. The effect on distribution shape is modest but still biases the population estimate upward.

The data loss from log transformation is not merely a sample size issue --- it is a systematic bias. The dropped values are not randomly distributed; they are concentrated at the low end of the flux distribution, where both near-zero emissions and true uptake occur. Any statistical analysis performed on the log-transformed subset (regression, ANOVA, mixed models) describes only the behavior of above-zero fluxes and cannot be generalized to the population of all flux measurements. A regression of log(CH4 flux) against temperature, for example, would describe "how large positive fluxes vary with temperature" rather than "how CH4 flux varies with temperature" --- a fundamentally different question.

### 7.2 The arcsinh transformation

The inverse hyperbolic sine, arcsinh(x) = ln(x + sqrt(x^2 + 1)), provides a variance-stabilizing transformation that handles the full real line. Its key properties for flux data are:

1. **Defined for all real numbers.** Unlike log, arcsinh accepts negative values, zero, and positive values without modification or data loss.

2. **Approximately linear near zero.** For |x| << 1, arcsinh(x) ~ x, preserving the scale and sign of small fluxes, including the distinction between small positive emissions and small negative uptake.

3. **Approximately logarithmic for large |x|.** For |x| >> 1, arcsinh(x) ~ sign(x) * ln(2|x|), compressing long tails in both the positive and negative directions. This provides the same variance-stabilizing benefit as the log transform for large fluxes, without its limitations for small or negative values.

4. **Smooth transition.** The shift from linear to logarithmic behavior is continuous, with no threshold or breakpoint to specify.

5. **Symmetric.** arcsinh(-x) = -arcsinh(x), so the transformation treats positive and negative fluxes symmetrically. A distribution of fluxes centered near zero but including both sources and sinks remains centered near zero after transformation.

The arcsinh transformation is closely related to the "generalized logarithmic" or "started logarithmic" transforms sometimes used in ecology (log(x + c) for some constant c), but avoids the arbitrary choice of the shift constant c and handles negative values naturally.

Fig. 11 (arcsinh behavior) shows the transformation function: it is indistinguishable from the identity function near the origin and converges to sign(x) * log(2|x|) for large arguments.

### 7.3 Empirical comparison

**Distribution shape (Fig. 10).** For the wetland soil dataset, the raw distribution is strongly left-skewed with a long right tail (Shapiro-Wilk W = 0.15, p < 10^-33). The log10 transform produces a roughly normal distribution (W = 0.97, p = 0.33) but retains only 17% of the data. The arcsinh transform retains all 288 values and substantially improves symmetry (W = 0.83, p < 10^-16) --- not perfectly normal, but far more tractable for parametric methods than the raw data, and without discarding the majority of the observations.

**Q-Q plots (Fig. 12).** For the wetland soil data, the raw Q-Q plot shows extreme departure from normality with heavy tails in both directions. The log10 Q-Q plot appears well-behaved but represents only the 49 positive values. The arcsinh Q-Q plot shows moderate departure in the tails (driven by a few large emission events) but tracks the normal reference line through the central 80% of the distribution.

**Variance stabilization.** Across the three datasets, raw variance spans a 1334-fold range (wetland soil: 139, tree canopy: 0.57, tree stem: 0.10). After arcsinh transformation, the variance ratio drops to 29-fold (wetland soil: 1.06, tree canopy: 0.23, tree stem: 0.04). This improvement is critical for mixed models or ANOVAs that compare flux rates across ecosystem components, where heteroscedasticity in raw data would violate model assumptions.

### 7.4 Practical recommendations for statistical analysis

The arcsinh transformation substantially improves variance homogeneity and distributional symmetry relative to raw data, but it does not fully normalize CH4 flux distributions --- particularly for datasets with genuine bimodality (e.g., a mixture of uptake-dominated and emission-dominated plots) or extreme outlier events. It is best understood as a variance-stabilizing and visualization tool, not as a route to normality. For formal statistical inference, non-parametric and model-based approaches that do not assume normality are often more appropriate.

**Arcsinh for visualization and variance stabilization.** Use arcsinh-scaled axes for plotting (as in the ridgeline and scatter plots throughout this analysis) and as a transformation in exploratory analysis. Its primary advantages are preserving negative values, stabilizing variance across groups (1334-fold to 29-fold in our datasets), and compressing long tails without arbitrary constants. Where a transformation is needed for model fitting (e.g., to improve residual behavior in a linear mixed model), arcsinh is preferable to log because it does not exclude data.

**Non-parametric methods for group comparisons.** For testing whether fluxes differ between groups (species, treatments, seasons, landscape positions), non-parametric approaches avoid distributional assumptions entirely: Wilcoxon rank-sum tests for two-group comparisons, Kruskal-Wallis tests for multiple groups, and permutation tests for more complex designs. These methods use all observations regardless of sign or magnitude and are robust to the heavy-tailed, mixed-sign distributions typical of CH4 flux data.

**Bootstrapped confidence intervals for means and medians.** For estimating population parameters and their uncertainty, bootstrap resampling (e.g., BCa intervals) provides valid confidence intervals without requiring normality. This is particularly useful for mean flux estimates used in scaling (Section 6.1), where the arithmetic mean of untransformed data is the quantity of interest but the distribution is too skewed for standard parametric CIs.

**Generalized linear models for regression.** When modeling CH4 flux as a function of continuous predictors (temperature, moisture, stem diameter), generalized linear models (GLMs) or generalized linear mixed models (GLMMs) can accommodate non-normal error distributions directly. Gaussian GLMs with an identity link on arcsinh-transformed data are one option; alternatively, Tweedie or gamma families handle right-skewed positive data, though these require positive responses and thus face the same exclusion problem as log transforms for mixed-sign data. For mixed-sign flux data, a Gaussian family on arcsinh-transformed responses, or a two-part (hurdle) model separating uptake from emission, may be most appropriate.

**If log-transforming, report the data loss explicitly.** State the number and percentage of values excluded, their distribution (are they noise-scattered near-zero values, or true biological uptake?), and whether the remaining subset is representative of the population. A statement like "we log-transformed CH4 fluxes, excluding 31% of measurements with non-positive values" is informative; silently dropping negatives is not.

**For back-transformation and reporting:** arcsinh-transformed means can be back-transformed via sinh() for interpretation in original units. However, as with log-transformed data, the back-transformed mean of the transformed data is not equal to the arithmetic mean of the original data (it approximates the geometric mean). For unbiased population means, use untransformed data (Section 6.1). Use transformations for variance stabilization in models, not for computing summary statistics.

**Interaction with quality filtering:** The arcsinh transformation and the "retain all values" recommendation (Section 6.1) are complementary. Retaining below-MDF measurements preserves the population distribution; arcsinh transformation makes that distribution more tractable for visualization and modeling. Conversely, if data are first filtered to remove below-MDF values and then log-transformed, both steps independently remove near-zero and negative values, compounding the positive bias documented in Section 5.3.4.

---

### Figure 10. Transformation comparison across CH4 flux datasets

![ch4_transform_comparison](plots/precision/ch4_transform_comparison.png)

**Caption.** Histograms of CH4 flux under three transformations (columns: raw, log10, arcsinh) for three datasets (rows: wetland soil, tree stem, tree canopy). Sample sizes in each panel report the number of valid (non-NA) values; for the log10 column, dropped values and their percentage are noted. The wetland soil dataset loses 83% of its measurements under log10 because the majority of fluxes are negative (CH4 uptake). Tree stems lose 31% and tree canopy 10%. The arcsinh transformation retains all values across all datasets while compressing the long tails that make raw distributions difficult to analyze. Dashed vertical lines mark zero.

---

### Figure 11. arcsinh transformation behavior

![arcsinh_behavior](plots/precision/arcsinh_behavior.png)

**Caption.** Properties of the arcsinh function. Left: arcsinh(x) (blue solid) compared to the logarithmic approximation sign(x) * log(2|x|) (red dashed) over the range [-200, 200]. The two converge for |x| > ~10. Right: near the origin, arcsinh(x) is approximately equal to x (green dashed line), meaning small fluxes --- including the critical distinction between small positive (emission) and small negative (uptake) values --- are preserved on their original scale without compression.

---

### Figure 12. Q-Q plots under three transformations (wetland soil CH4)

![ch4_transform_qq](plots/precision/ch4_transform_qq.png)

**Caption.** Normal Q-Q plots for wetland soil CH4 flux under raw (left), log10 (center), and arcsinh (right) transformations. Red lines show the best-fit linear reference. The raw distribution departs severely from normality with heavy tails. The log10 Q-Q plot appears approximately normal but represents only the 49 positive values (17% of 288 total measurements); the remaining 83% of the distribution --- all negative fluxes representing CH4 uptake --- has been silently excluded. The arcsinh Q-Q plot retains all 288 values and tracks the normal reference through the central portion of the distribution, with moderate departure only in the extreme tails (driven by a few large emission or uptake events).
