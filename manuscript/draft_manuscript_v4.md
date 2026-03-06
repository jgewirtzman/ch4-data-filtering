# Consequential choices and hidden bias in chamber CH₄ flux data processing: considerations for detection limits, filtering, and transformation

**[[ Authors ]]**

**[[ Affiliations / journal submission line ]]**

# Abstract

Chamber-based CH₄ flux measurements increasingly capture ecosystem components whose fluxes are small, near-detection, or intermittently negative, making data processing choices consequential for reported results. Researchers processing such data face a sequence of largely unconsidered decisions: how to estimate instrument detection limits, which measurements to retain, what to do with below-detection values, and how to transform data for statistical analysis. Each choice interacts with the others, and errors compound downstream. **[[ 2--3 sentences on empirical approach and datasets. ]]** Empirical CH₄ precision estimated from the continuous analyzer timeseries differed substantially from manufacturer specifications, with discrepancies ranging from four-fold worse to better-than-spec across instruments; this variation propagated directly into detection limit calculations, with three published minimal detectable flux (MDF) methods flagging 20--62% of measurements as below detection depending on the precision estimate and confidence level used. Sensitivity analysis across 13 quality filter criteria applied to a multi-tissue canopy dataset showed that stand-level CH₄ budgets were stable within 16% across all filter definitions, but individual component flux rates --- particularly leaves --- varied by up to 100% depending on filter choice and outlier treatment. Critically, removing below-detection measurements from the dataset inflated the mean flux by 17--135% and the median by 40--319% relative to the unfiltered baseline, with bias scaling monotonically with filter stringency; retaining original values and flagging detection status is the only treatment that preserves both statistics. Finally, the log transformation commonly applied to flux data silently drops 10--83% of measurements across ecosystem types because it cannot accommodate negative values, introducing a systematic upward bias; the arcsinh transformation is a direct replacement that handles negative values, stabilizes variance across components, and introduces no data loss. We provide question-dependent recommendations for each processing step and a reporting checklist to support cross-study comparability.

# 1. Introduction

Chamber-based measurements of CH₄ flux are used across a broad range of ecosystem components, from wetland soils with large, temporally variable emissions to tree stems, branches, and leaves whose fluxes are often small and close to the detection limit of field-deployable analyzers. The methodological literature has devoted considerable attention to the design of chamber systems, the selection of flux-fitting models, and the appropriate temporal resolution of measurements **[[ citations ]]**. Comparatively little attention has been paid to the sequence of data processing choices that follow raw flux computation: how to characterize instrument precision and set detection thresholds, which measurements to retain or discard, how to treat measurements that fall below detection, and how to transform flux distributions for statistical analysis.

These choices are not merely technical details. Each step in the processing pipeline introduces potential bias, and the choices interact: an overly conservative detection limit combined with removal of below-detection values and log transformation of the remainder compounds three independent sources of positive bias into the reported flux estimate. The magnitude of these biases depends on the fraction of measurements near detection --- which in turn depends on the true flux magnitude, the instrument precision, and the chamber geometry --- but can reach or exceed 100% of the true flux for small-flux components such as leaf or dormant-season stem emissions.

Three interconnected issues motivate this paper. First, minimal detectable flux (MDF) calculations require a precision estimate, and manufacturer specifications are routinely used as a convenient default. We show that empirical precision estimated from the analyzer timeseries can differ substantially from the datasheet value in either direction, and that this discrepancy propagates directly into detection thresholds and the fraction of measurements flagged as below detection. Second, quality filtering involves two distinct choices that are often conflated: the choice of detection criterion and the treatment of measurements that fail it. We show that the treatment --- whether to remove below-detection measurements, set them to zero, or retain them with a flag --- is more consequential for reported means and medians than the choice of detection threshold itself. Third, the log transformation applied to flux data for variance stabilization and normalization silently discards negative values, which are common in mixed-sign distributions such as those from methanotrophic soils or near-detection stem fluxes; the arcsinh transformation provides a drop-in replacement that handles the full real line without data loss.

We address each issue empirically using chamber CH₄ flux datasets from Harvard Forest and Yale-Myers Forest, spanning ecosystem components with very different flux magnitudes and negative-value prevalence. **[[ Brief description of datasets and instruments introduced here --- expand once site/instrument methods are finalized. ]]** Our goal is not to prescribe a single processing protocol, but to document the consequences of each choice and provide question-dependent recommendations that allow researchers to select the approach appropriate for their scientific objective.

# 2. Methods

## 2.1 Study sites and datasets

**[[ Full site descriptions for Harvard Forest and Yale-Myers Forest to be added here. ]]** Three datasets are used throughout this analysis. The *Harvard Forest canopy dataset* (n = 136 measurements) spans four tissue compartments --- stem below 2 m, stem above 2 m, branches, and leaves --- sampled using ABB GLA131-series Microportable Greenhouse Gas Analyzers (UGGA; hereafter GLA131). The *Harvard Forest stem monitoring dataset* (n = 1,640 measurements) comprises stem flux measurements from 2023--2024 using GLA131 instruments (n = 1,142) and 2025 using the LI-COR LI-7810 (n = 498). The *Yale-Myers Forest datasets* include stem flux measurements (n = 369) and wetland soil flux measurements (n = 288) collected with GLA131 instruments; these are used primarily for the data transformation comparison (Section 5). **[[ Add measurement dates, tree species, collar installation details, and any other site-level metadata needed for methods reproducibility. ]]**

## 2.2 Instruments and chamber geometry

GLA131-series analyzers use off-axis integrated cavity output spectroscopy (OA-ICOS) and report dry-mole-fraction CH₄ and CO₂ at 1 Hz. The manufacturer specifies CH₄ precision of 0.9 ppb and CO₂ precision of 0.35 ppm (1-sigma at 1 s averaging). The LI-COR LI-7810 is a cavity ring-down spectrometer with manufacturer-specified CH₄ precision of 0.6 ppb (1-sigma at 1 s).

All stem measurements used permanent collars with radius 5.08 cm (surface area = 0.00811 m²). Chamber volume was computed from measured tree diameters and collar depth for each tree, with an additional 0.028 L for connecting tubing. The flux term --- encoding ambient pressure, chamber volume, surface area, and temperature via the ideal gas law --- was computed per measurement using the *goFlux* R package **[[ (citation) ]]**. For the Harvard Forest canopy dataset, chamber dimensions varied by tissue type; **[[ add canopy chamber dimensions and any tissue-specific geometric corrections here. ]]**

## 2.3 Empirical precision estimation

### 2.3.1 Allan deviation (per-measurement)

For each flux measurement, we estimated Allan deviation at τ = 1 s from the dry-mole-fraction timeseries within the measurement window. Because concentrations change monotonically during a chamber closure, the flux-driven trend was removed by first-differencing before estimating noise, yielding a measurement-specific precision estimate that is independent of the flux model fit.

### 2.3.2 Rolling window scan (per-instrument)

We scanned the full continuous timeseries for each instrument using a 30-second rolling window. For each window we computed the absolute slope (via linear regression) and the residual standard deviation around the fit. Windows were ranked by a combined normalized score of \|slope\| + residual SD; the bottom 5% (flattest and quietest) were selected as representative of the true instrument noise floor. This yields a single global precision estimate per instrument that is used in the Wassmann MDF calculation (Section 2.4).

## 2.4 Minimal detectable flux (MDF) calculation approaches

We computed MDF using three published approaches, all sharing the general form:

MDF = precision_term / t × flux_term

where *t* is measurement duration in seconds and flux_term encodes chamber volume-to-area ratio and temperature via the ideal gas law, computed per measurement by *goFlux*. The three approaches differ in their precision term and confidence scaling:

-   **Manufacturer MDF (goFlux).** Uses the GLA131 datasheet precision (0.9 ppb CH₄): MDF = precisionₘₐₙᵤᶠᵃᶜᵗᵘʳᵉʳ / t × flux_term.

-   **Wassmann et al. (2018).** Uses the global empirical precision from the rolling window method, scaled for the desired confidence level: MDF = (z × SDᵍˡᵒᵇᵃˡ) / t × flux_term, where z = 2.576, 1.960, or 1.645 for 99%, 95%, and 90% confidence respectively.

-   **Christiansen et al. (2015) MQL variant.** Uses the per-measurement Allan deviation, a factor of 3, and the t-distribution critical value: MDF = (SDₐˡˡᵃₙ × 3 × tᶜʳᵉᵗ) / t × flux_term, where tᶜʳᵉᵗ is the two-tailed t critical value at the specified confidence level with df = t − 2. Computed at 99%, 95%, and 90% confidence.

**[[ Add brief note on measurement duration effects: MDF scales with 1/t; figure showing MDF vs. t at empirical precision levels for representative chamber geometry would go here or in supplementary material. ]]**

## 2.5 Additional quality filter criteria

Beyond MDF, we evaluated several additional quality filters commonly applied to chamber flux data:

-   **CH₄ R² thresholds ( 0.7,  0.9):** Coefficient of determination from the best-fit model (linear or Hutchinson-Mosier, as selected by goFlux).

-   **CO₂ R² threshold ( 0.7):** Applied to the concurrent CO₂ flux from the same measurement as an independent check on measurement quality.

-   **Signal-to-noise ratio (SNR) thresholds ( 2,  3):** Defined as \|CH₄ flux\| / noise_floor, where noise_floor = σₐˡˡᵃₙ / t × flux_term converts the per-measurement Allan deviation to flux units. This empirical, per-measurement definition captures instrument-specific and condition-specific noise rather than relying on manufacturer specifications.

## 2.6 Quality filter sensitivity analysis

To assess how the choice of quality filter propagates through component and whole-tree flux scaling, we applied each of 13 criteria independently to the Harvard Forest canopy dataset (n = 136 measurements). For each filter, measurements failing the criterion had their CH₄ flux set to zero; passing measurements retained their original flux values. We then computed mean per-tissue flux rates for four compartments (stem \< 2 m, stem ≥ 2 m, branches, leaves) and scaled these by Whittaker and Woodwell (1967) surface area indices (stem bark = 0.45, branch bark = 1.70, LAI = 4.5 m² m⁻² ground) to obtain integrated stand-level budgets. A cone taper model (26 m tree, 40 cm DBH) was used to partition stem bark area into below- and above-2 m fractions. The analysis was repeated with outliers removed (measurements exceeding Q3 + 3 × IQR per compartment) to assess the combined effect of quality filtering and outlier sensitivity.

## 2.7 Below-MDF treatment comparison

For each of seven MDF thresholds (Manufacturer, Wassmann 90/95/99%, Christiansen 90/95/99%), we applied three treatments to measurements flagged as below detection: (1) *Remove*: exclude below-MDF measurements from the dataset entirely; (2) *Set to zero*: replace below-MDF flux values with zero, preserving sample size; and (3) *Keep original (flag only)*: retain all measured flux values and flag each measurement as above or below its per-measurement MDF. We compared the resulting means and medians against the unfiltered baseline, pooled and split by instrument, using the Harvard Forest stem monitoring dataset (n = 1,640).

## 2.8 Data transformations

We compared three transformations --- raw (untransformed), log₁₀, and arcsinh --- across three datasets spanning a range of negative-flux prevalence (wetland soil, tree stem, tree canopy). For each transformation we report data retention (fraction of measurements with defined transformed values), the Shapiro-Wilk W statistic as a measure of distributional symmetry, and the variance ratio across components as a measure of heteroscedasticity. **[[ Add any additional statistical details on the transformation comparison here, e.g. whether formal tests of variance homogeneity were applied. ]]**

**[[ ADDITIONAL METHODS NEEDED: (1) Temporal autocorrelation / degrees of freedom treatment --- brief justification for df = t − 2 in Christiansen MDF, or note on whether this was examined. (2) Model selection (linear vs. Hutchinson-Mosier) --- note on how goFlux selects models and whether model choice interacts with filter outcomes. (3) Any seasonal or environmental covariates on instrument precision --- was Allan deviation stable across seasons / temperatures? (4) Cross-instrument Wassmann SD pooling decision --- expand the justification for computing instrument-specific vs. pooled global SD. ]]**

# 3. Detection limits

## 3.1 Empirical precision vs. manufacturer specification

Empirical CH₄ precision varied substantially across instruments and estimation methods (Fig. 1). For LGR1, the median Allan deviation (3.6 ppb) was four times the manufacturer specification (0.9 ppb), and the rolling window estimate (2.5 ppb) was nearly three times the specification. LGR3 showed a more moderate discrepancy (Allan deviation 1.5 ppb, rolling window 1.2 ppb). LGR2 was the exception: both empirical estimates (Allan deviation 0.6 ppb, rolling window 0.7 ppb) were better than the datasheet value. The Yale-Myers instrument (LGR3 redeployed) yielded empirical precision of 1.2--1.3 ppb, consistent with its Harvard Forest performance.

These results confirm that manufacturer specifications can substantially underestimate real-world noise for some instruments, while overestimating it for others. The direction of discrepancy is not predictable from instrument type or age alone. Relying on a single datasheet value across instruments therefore introduces uncontrolled, instrument-specific variability into detection limits. CO₂ precision showed a similar pattern: empirical estimates ranged from 0.5 to 1.8 ppm across instruments, compared to the manufacturer specification of 0.35 ppm.

In the multi-instrument stem monitoring dataset, the GLA131 median Allan deviation across 1,206 measurements was 4.0 ppb (IQR: 3.2--4.9 ppb), 4.4× worse than the manufacturer specification. The LI-7810 median was 0.095 ppb (IQR: 0.086--0.108 ppb), 6.3× better than its specification, with a remarkably tight distribution spanning less than a factor of 1.3 across the interquartile range (Fig. 5). The two distributions were completely non-overlapping on a log scale, representing an approximately 40-fold difference in real-world precision --- far larger than the 1.5-fold difference implied by the manufacturer specifications (0.9 vs. 0.6 ppb).

**[[ Figure 1 reference: CH4 precision manufacturer vs. empirical bar chart (already produced). Figure 5 reference: per-measurement Allan deviation violin plot by instrument (already produced). ]]**

## 3.2 Consequences for MDF thresholds and detection rates

The three MDF approaches produced detection thresholds spanning roughly an order of magnitude for CH₄ (Fig. **[[ X ]]**). The manufacturer-based MDF was the most permissive, flagging 20% of CH₄ measurements as below detection. Wassmann thresholds (which substitute empirical global precision for the datasheet value) flagged 26--34%, depending on confidence level. The Christiansen approach --- which combines per-measurement Allan deviation with a factor of 3 and a t-distribution critical value, yielding effective multipliers of approximately 5--8× the Allan deviation --- was the most conservative, flagging 52--62% as below detection. Nearly all CO₂ fluxes (95%) exceeded detection thresholds under every method.

In the multi-instrument stem dataset, instrument-specific Wassmann global SDs (GLA131: 3.96 ppb, LI-7810: 0.095 ppb) produced thresholds roughly 40× higher for GLA131 measurements. Consequently, Wassmann 99% flagged 46% (531/1,142) of GLA131 measurements but only 15% (74/498) of LI-7810 measurements as below detection. This demonstrates that when empirical precision is used --- as it should be --- the fraction of below-detection measurements is itself instrument-dependent and cannot be treated as a property of the flux distribution alone.

**[[ ADDITIONAL ANALYSIS NEEDED: Figure showing MDF as a function of measurement duration t at empirical precision levels for representative chamber geometry. This would concretely illustrate the tradeoff between duration and detectability and give readers a practical tool for study design. ]]**

## 3.3 Negative fluxes: noise or biology?

The stem flux dataset contained 284 negative CH₄ flux measurements (17.3% of 1,640), distributed asymmetrically between instruments: 275/1,142 (24.1%) for the GLA131 vs. 9/498 (1.8%) for the LI-7810. Three lines of evidence indicate that the majority of GLA131 negative fluxes reflect instrument noise rather than biological CH₄ uptake:

-   **Apples-to-apples seasonal comparison.** Restricting the comparison to overlapping calendar months (June--October), the GLA131 still showed 20.3% negative fluxes (169/832) vs. 1.7% (6/353) for the LI-7810 --- same trees, same season, same chambers.

-   **Noise floor analysis.** Of 246 GLA131 negative fluxes with Allan deviation estimates, 52% had absolute values within 1σ of the per-measurement noise floor, 74% within 2σ, and 84% within 3σ. The median absolute negative flux (0.049 nmol m⁻² s⁻¹) was comparable to the median noise floor (0.051 nmol m⁻² s⁻¹) (Fig. 8).

-   **Distribution shape.** On the GLA131, the near-zero flux distribution was roughly symmetric around a small positive value, consistent with Gaussian noise broadening a weakly positive signal. On the LI-7810, the same trees resolved as a tight, positively-skewed distribution with almost no negative tail.

The GLA131 negative flux rate showed a seasonal pattern --- higher in dormant months (35.7% November--April) than growing months (20.5% May--October) --- consistent with real fluxes approaching zero during dormancy. But even in the growing season, 20% negative fluxes from the GLA131 vs. \<2% from the LI-7810 for the same months indicates that the dominant source of negative values is noise scatter, not biology. The seasonal pattern likely reflects real fluxes approaching zero (bringing more of the distribution within the noise wedge) rather than true biological CH₄ uptake.

This finding has direct implications for quality filtering: sign-based exclusion of negative fluxes is not appropriate without instrument-specific noise characterization, because the majority of negative values from a noisy instrument may be statistically indistinguishable from zero and should not be interpreted as uptake.

# 4. Quality filtering

## 4.1 Two distinct quality questions

Quality filtering of chamber flux data conflates two genuinely distinct questions that require different tools and different logic. The first is a question of *measurement validity*: was this closure conducted correctly, with no leaks, equipment malfunctions, or procedural errors that would make the measured concentration change unrepresentative of the true surface flux? The second is a question of *signal detectability*: was the flux large enough relative to instrument noise to be distinguishable from zero? These are independent. A perfectly valid closure can yield a below-detection flux (common for leaves, dormant stems, and low-emitting surfaces). A corrupted closure can produce a concentration change that passes every statistical detection criterion. Conflating the two --- using R² or MDF thresholds as a proxy for both simultaneously --- is where most quality-filtering confusion originates.

The distinction is particularly important for CH₄ relative to CO₂. For CO₂ from living tissue, fluxes are almost always large relative to instrument noise, so detectability is rarely in question; the interesting quality problem is entirely one of validity. For CH₄, both questions are live simultaneously --- the flux may be near-detection even from a perfectly valid closure --- and the appropriate tools for each question are different. The same considerations apply to N₂O and other trace gases measured at low mixing ratios. Throughout this paper we focus primarily on the detectability question (Sections 3--5), but the validity question must be addressed first, before detection-limit analyses are meaningful.

## 4.2 Assessing measurement validity

The most useful diagnostic for detecting invalid closures is **CO₂ flux as a process tracer for live tissue**. Live respiring tissue must produce a positive CO₂ flux; a zero or negative CO₂ rate is strong evidence of a seal failure, collapsed chamber, or dead tissue, regardless of what the CH₄ signal shows. This is the basis for the CO₂ R² filter included in our analysis, which is better understood as a validity criterion than a detection criterion: it tests whether a coherent respiratory signal was present, not whether the CH₄ flux was large enough to detect. Note, however, that CO₂ as a validity tracer fails for surfaces without net respiration --- standing water, bare mineral soil in some conditions, dead wood --- where a valid closure may produce near-zero CO₂ flux.

A critical limitation of R²-based validity checks is that *a constant leak produces a perfectly linear concentration change*. If ambient CH₄ concentration differs from the chamber headspace, a steady leak will drive a monotonic concentration change indistinguishable from a real flux by R² alone --- and will yield a high R² for both CH₄ and CO₂ if the leak is large. Detecting leaks therefore requires direct physical evidence, not statistical criteria applied to the concentration timeseries. For **manual chamber systems**, this means: physical leak testing at collar installation (pressurization, soap bubble tests), consistent chamber seating protocols, field notes recording any observed seal problems or disturbances, and operator familiarity with expected concentration dynamics for the surface type being measured. For **automated chamber systems**, additional automated checks become feasible: lid position sensors to confirm full closure, comparison of pre- and post-closure ambient concentrations for consistency, tracer gas injection with recovery monitoring (e.g., SF₆ or a CH₄ spike) to quantify leak rates per closure, and outlier detection across the repeated-closure time series for each collar position. These automated checks substantially reduce the irreducible uncertainty that manual systems face.

There is an irreducible challenge for low-flux CH₄ measurements on manual systems: a small but steady leak on a surface with a near-detection flux may produce a concentration change that (a) passes R² thresholds, (b) exceeds MDF, and (c) is not flagged by CO₂ if the surface has low or variable respiration (e.g., dormant woody tissue, recently disturbed soil). In this regime, no post-hoc statistical filter can reliably distinguish a valid below-detection flux from a corrupted above-detection artifact. This is a genuine limitation that should be acknowledged rather than obscured by the application of filters that were not designed for this purpose. The appropriate response is investment in physical system integrity --- collar installation quality, seating protocols, and repeated leak checks --- rather than reliance on concentration-timeseries statistics as a proxy for validity.

A further and distinct source of near-zero and below-detection fluxes is genuine biological suppression: fluxes that are small not because of instrument noise or measurement error, but because the true emission rate is low. Temperate soil CH₄ and CO₂ fluxes approach zero or become negligible during winter, drought, or other periods of biological dormancy. Woody tissue CH₄ emissions may drop near zero in the dormant season. Methanotrophic soils may shift from net uptake to near-zero net exchange under certain conditions. In all of these cases, a low R², a small SNR, and a below-MDF flag are expected and correct outcomes --- the measurement is valid, the instrument is functioning, and the near-zero result is the scientific finding. Applying a quality filter that removes or zeros such measurements does not correct for noise; it removes a biologically real state from the dataset. The consequence is not random attenuation but structured, ecologically correlated bias: excluded measurements are systematically associated with cold temperatures, dry conditions, or dormant seasons, so the filtered dataset over-represents high-flux conditions. For annual or multi-season budgets, this is a direct and predictable overestimate. The distinction between "near-zero because of noise" and "near-zero because of biology" cannot be resolved by any statistical criterion applied to the concentration timeseries alone; it requires knowledge of the ecosystem, the season, and the measurement context.

## 4.3 Sensitivity of component flux rates to filter choice

Applying 13 quality filters independently to the Harvard Forest canopy dataset (n = 136) revealed that stem flux rates were robust to filter definition, while branch and leaf rates were more sensitive (Fig. 2, Fig. 4).

Stem \< 2 m fluxes --- the largest component (\~1.0 nmol m⁻² s⁻¹) --- varied by only 7% across all filter definitions (range 93--100% of the no-filter baseline, CV = 2.2%). These fluxes are large relative to any detection threshold and are therefore insensitive to the choice of quality criterion. Stem  2 m fluxes (range 88--100%, CV = 4.8%) and branch fluxes (range 82--100%, CV = 5.0%) showed moderate sensitivity. The most restrictive filters (Christiansen 99%, CH₄ R²  0.9) reduced mean branch flux by up to 18%, primarily by zeroing measurements near the detection limit.

Leaf fluxes were the most sensitive component (range 78--100%, CV = 9.0% with outliers retained). Because leaf CH₄ fluxes are small and many individual measurements fall near or below detection, the choice of filter criterion has a proportionally larger effect on the component mean. With outliers removed, leaf sensitivity increased dramatically (0--100% range), reflecting dependence on a single high-flux measurement that, when retained, anchors the leaf mean above zero across all filter definitions.

A clear hierarchy emerges: filter choice matters in proportion to how close a component's flux is to the detection limit. Stem \< 2 m fluxes, which are 5--10× larger than any MDF threshold, are effectively immune to the choice of quality criterion. Leaf fluxes, by contrast, are comparable in magnitude to the detection limit and are therefore highly sensitive.

## 4.4 Sensitivity of the stand-level budget

Despite the variable sensitivity of individual components, the integrated stand-level CH₄ budget was remarkably stable across filter definitions (Fig. 3). With outliers retained, the total budget ranged from 0.692 nmol m⁻² ground s⁻¹ (no filter) to 0.582 nmol m⁻² ground s⁻¹ (Christiansen 99%), a maximum reduction of 16%. This stability reflects the dominance of stem \< 2 m fluxes --- which are insensitive to filtering --- in the integrated budget after scaling by surface area indices.

With outliers removed, the budget ranged from 0.539 to 0.425 nmol m⁻² ground s⁻¹ (maximum reduction 21%). The larger relative sensitivity when outliers are removed arises for two reasons: first, most outliers are high-flux stem  2 m measurements, so removing them makes that component substantially more sensitive to subsequent filtering (CV roughly doubles from 4.8% to 9.9%); second, outlier removal shrinks the total budget baseline, so the same absolute filtering effect represents a larger percentage reduction. Outlier removal was therefore more consequential for the integrated budget than the choice of detection methodology.

The empirical SNR filters (SNR  2 and SNR  3) produced intermediate results: 96% and 95% of the unfiltered budget, respectively, placing them between the Wassmann and Christiansen thresholds in stringency.

## 4.5 R² as a detection metric, not a quality filter

R² thresholds deserve separate consideration because they interact with the two question types differently. R² measures how well a model fits the concentration timeseries, which is a function of both the flux magnitude and the instrument noise. A measurement with R²  0.7 is one where the concentration trend was clear above the noise --- a detection criterion in the same spirit as the MDF threshold, but with a key difference: R² is not independent of flux magnitude.

Two measurements with the same instrument noise but different true fluxes will have different R² values; the larger flux always produces a higher R². This means R² filters preferentially retain high-flux measurements and systematically exclude low emitters. CH₄ R²  0.5 and Christiansen 99% each retained the same number of measurements (665/1,640, 40.5%), but with very different compositions: R²  0.5 retained only 1.4% negative fluxes, while Christiansen 99% allowed 5.6% negatives. A tree with a true flux of 0.02 nmol m⁻² s⁻¹ measured on a GLA131 with 4 ppb noise will produce a concentration change of \~0.3 ppb over a 5-minute measurement --- far below the noise floor, and thus never yielding a high R², even though it is a valid measurement of a small flux.

R² is most useful as a confidence flag for individual measurement interpretation: "this particular measurement showed a clear monotonic trend." It is least useful as a population filter for computing means or budgets, where it introduces systematic selection bias toward high-flux measurements.

Critically, R² filters do not merely bias toward high-flux *measurements* --- they bias toward high-flux *conditions*. Because true flux magnitude determines R², an R² threshold will preferentially retain measurements made during warm, wet, or biologically active periods and preferentially exclude measurements from cold, dry, or dormant periods when fluxes are genuinely small. This is a form of seasonal or environmental confounding that operates independently of instrument noise: even with a perfect instrument, a R² threshold would exclude winter soil measurements, drought-period emissions, and dormant-season stem fluxes simply because the biological rates are low. For any analysis involving annual budgets, seasonal comparisons, or environmental driver attribution, this structured exclusion of low-flux conditions introduces a systematic overestimate that no amount of careful gap-filling or scaling can fully correct, because the excluded data are not missing at random --- they represent a coherent ecological state.

## 4.6 Question-dependent filtering recommendations

The sensitivity analysis reveals that no single filtering approach is universally optimal. The appropriate treatment of below-detection measurements depends fundamentally on the scientific question being asked (Table 1). We identify two broad categories:

**Scaling and budgeting questions** (e.g., what is the mean flux or flux distribution for this ecosystem component; stand-level CH₄ budgets; seasonal emission totals; cross-site comparisons): retain all flux values and flag detection status. Removing below-MDF measurements introduces systematic positive bias that can exceed 100% at stringent thresholds. Setting below-MDF values to zero avoids this bias in the mean but distorts the median.

**Detection and biological signal questions** (e.g., is the flux from this surface significantly nonzero; is a given surface a net source or sink; does flux differ significantly between treatments): filter to above-MDF measurements only. Including below-detection measurements in a test of "is the flux nonzero?" adds noise without adding information. The positive bias from removal is, in this context, expected behavior: the relevant subset is measurements where the instrument could actually resolve the signal.

Table 1. Question-dependent filtering recommendations.

| Question / goal | Objective | Below-MDF treatment | Appropriate filters |
|---|---|---|---|
| Mean flux / budgets / scaling | Unbiased population estimate | Keep original values + flag | Hard QC only (seal failures, equipment errors) |
| Seasonal / spatial patterns in mean flux | Unbiased comparison of distributions | Keep original values + flag | Hard QC only; report fraction below MDF per group |
| Is this flux nonzero? | Detection confidence | Remove below MDF | Christiansen MDF at 90/95/99% |
| Is this a source or sink? | Sign determination | Remove below MDF | Christiansen MDF + consider sign test on above-MDF subset |
| Which trees are big emitters? | Identify high-flux individuals | Remove below MDF (or R² filter) | R² > 0.5 or Allan SNR > 3 appropriate here |
| Model fitting / driver analysis | Confident flux estimates for regression | Remove below MDF | Christiansen or Allan SNR; report selection and fraction retained. If below-MDF fraction is high, consider censored regression with per-measurement MDF as censoring threshold. Do not set to zero. |
| Cross-study comparison | Comparable methodology | Report all three treatments | Show sensitivity; flag instrument and empirical precision |

**[[ ADDITIONAL RECOMMENDATION NEEDED: A reporting checklist --- what to include in any chamber CH₄ methods section to enable cross-study comparison: instrument model, empirical precision estimate method, n measurements, fraction below MDF, filter applied, below-MDF treatment, transformation used. Could be a second small table or a numbered list here or in Section 6. ]]**

The table above does not list "set below-MDF values to zero" as a primary recommendation for any question type. This is intentional. Setting below-MDF values to zero preserves the arithmetic mean (because symmetric noise cancels), but it creates a heterogeneous error structure: measurements well above detection retain their full symmetric noise distribution, while near-detection measurements have their lower tail truncated and replaced with a point mass at zero. This selective distortion of the residual distribution is inconsequential for simple mean estimation but is problematic for regression, where OLS assumes a homogeneous error structure across all observations. The set-to-zero treatment is retained as a sensitivity analysis option --- to be shown alongside the keep-original and remove treatments for key reported results --- but should not be the primary treatment for any inferential analysis.

# 5. Treatment of below-detection measurements

## 5.1 Statistical consequences of the three treatments

Across seven MDF thresholds applied to the Harvard Forest stem monitoring dataset (n = 1,640), the three below-MDF treatments produced markedly different outcomes for both the mean and median CH₄ flux (Figs. 6, 7).

Removing below-MDF measurements shifted the mean flux upward by 17% (Manufacturer MDF) to 135% (Christiansen 99%) relative to the unfiltered baseline. The median was even more severely affected: +40% to +319%. This positive bias arises because removal selectively eliminates near-zero measurements from both sides of the distribution --- both small positive and small negative fluxes --- leaving a distribution enriched in larger positive values. The bias scales monotonically with filter stringency because more stringent filters flag a larger fraction of near-zero values.

Setting below-MDF measurements to zero changed the mean by less than 5% at all thresholds, because symmetric noise scattering around zero cancels in the arithmetic mean. However, the median dropped toward zero and reached exactly zero for all Christiansen thresholds, where more than 50% of measurements were flagged as below detection. This reflects an artificial accumulation of values at exactly zero that was not present in the original distribution and that would distort any statistical analysis of the median or distribution shape.

Retaining all measurements at their original values --- flagging but not modifying those below detection --- is the only treatment that preserves both the mean and median of the flux distribution.

## 5.2 Why removal looks reasonable but is not, for population estimates

The intuition behind removal is straightforward: "if I cannot confidently measure a flux, I should not include it." This logic is correct for detection questions (Section 4.6) but is wrong for questions about population means or flux budgets.

The key insight is that instrument noise scatters symmetrically around the true flux value. For a tree emitting 0.01 nmol m⁻² s⁻¹, individual measurements on a noisy analyzer will scatter across a range that includes negative values --- but the mean of those measurements will converge on the true flux. Removing the below-detection scatter (which includes both the negative and the small-positive values) selectively retains the high-scatter tail, inflating the mean. This is analogous to the well-documented problem of non-detect treatment in analytical chemistry (Helsel 2005), where substituting values for censored observations --- including zero substitution --- is recognized as biased, and methods that use the full dataset are preferred.

The practical recommendation is a tiered approach: (1) *Hard QC* --- remove measurements with clear evidence of equipment failure or procedural error (chamber seal failures evidenced by negative CO₂ flux from live tissue, equipment malfunctions noted in field logs, extreme Allan deviation outliers indicating instrument instability). These are quality flags based on what went wrong, not statistical thresholds on the flux value. (2) *Detection flagging* --- compute the per-measurement MDF using the Christiansen method with empirical Allan deviation; flag each measurement as above or below its individual detection limit; report the fraction below detection; retain all flux values. (3) *Sensitivity analysis* --- for key reported results, demonstrate robustness by showing results under all three treatments. If a conclusion changes depending on treatment, that is important information.

## 5.3 Instrument-dependence of the removal bias

The treatment effect was almost entirely confined to the GLA131. For the LI-7810, all three treatments produced nearly indistinguishable distributions at every MDF threshold (Fig. 9). For the GLA131, removal at Christiansen 99% inflated the instrument-specific mean from 1.44 to 4.61 nmol m⁻² s⁻¹ (a 220% increase), while the LI-7810 mean changed by only 51% (0.38 to 0.57 nmol m⁻² s⁻¹) under the same threshold.

This demonstrates that the removal bias is driven by instrument noise, not by the biological flux distribution. The LI-7810's high precision means very few measurements fall below even the most stringent MDF, so the three treatments converge. This has an important implication for multi-instrument datasets: pooling data across instruments of different precision and then applying a uniform removal criterion will preferentially distort the noisier instrument's data, producing a biased composite.

# 6. Data transformation

## 6.1 The problem with log transformation

The log transformation is the most common variance-stabilizing transformation in ecology and biogeochemistry. For strictly positive data, it compresses long right tails and stabilizes variance across groups with different means. However, the log of a non-positive number is undefined, and applying a log transform to CH₄ flux data silently drops all zero and negative values.

The data loss from log transformation is not merely a sample size issue --- it is a systematic bias. The dropped values are concentrated at the low end of the flux distribution, where both near-zero emissions and true uptake occur. Any statistical analysis performed on the log-transformed subset describes only the behavior of above-zero fluxes and cannot be generalized to the population of all flux measurements.

We demonstrate this with three CH₄ flux datasets spanning a range of negative-flux prevalence (Fig. 10). Wetland soil fluxes (Yale-Myers Forest, n = 288) are 83% negative, reflecting methanotrophic CH₄ consumption dominating at most plots; log₁₀ transformation retains only 49 of 288 measurements (17%). Tree stem fluxes (Yale-Myers Forest, n = 369) are 31% negative; log transformation drops 115 values. Tree canopy fluxes (Harvard Forest, n = 136) are 10% negative; log transformation drops 14 measurements. In each case, the dropped values are not random --- they are systematically concentrated at the low end of the flux distribution, biasing the transformed population upward relative to the true flux distribution.

This problem compounds with the removal bias documented in Section 5: if data are first filtered to remove below-MDF values and then log-transformed, both steps independently remove near-zero and negative values, compounding the positive bias.

Several alternative transformations have been proposed for data that are long-tailed and cross zero. The *signed log* \[sign(x) × log(\|x\| + 1)\] is common in older ecological literature and is simple to implement; the +1 softens behavior near zero and the sign is preserved, but the function is not smooth at the origin and its asymptotic properties are less clean than the arcsinh. The *Yeo-Johnson transform* is a parametric family (estimated via maximum likelihood) that generalizes Box-Cox to handle zero and negative values; it is useful when a data-driven choice of functional form is preferred over a fixed one, and is available in standard statistical software. Quantile or rank transforms map data to a target distribution nonparametrically and are completely robust to zeros, skew, and tails, but sacrifice interpretability of scale. For CH₄ flux data specifically --- where interpretability matters for reporting flux magnitudes and where the near-zero behavior is scientifically meaningful --- the arcsinh transformation is the current best-practice recommendation and has the largest supporting literature for this data type.

## 6.2 The arcsinh transformation

The inverse hyperbolic sine, arcsinh(x) = ln(x + √(x² + 1)), provides a variance-stabilizing transformation that handles the full real line. Its key properties for flux data are: (1) defined for all real numbers, accepting negative values, zero, and positive values without modification or data loss; (2) approximately linear near zero, preserving the scale and sign of small fluxes including the distinction between small positive emissions and small negative uptake; (3) approximately logarithmic for large \|x\|, compressing long tails in both directions; (4) a smooth transition between these regimes with no threshold or breakpoint to specify; and (5) symmetric --- arcsinh(−x) = −arcsinh(x) --- so distributions centered near zero remain centered near zero after transformation (Fig. 11).

The arcsinh transformation is related to the "started logarithmic" transforms sometimes used in ecology \[log(x + c) for some constant c\], but avoids the arbitrary choice of the shift constant c and handles negative values naturally by construction rather than by adding an offset.

A scaled variant, arcsinh(x/theta), allows the sensitivity near zero to be tuned by the scaling parameter theta. For CH₄ flux data the near-zero behavior depends on measurement units (nmol vs. umol vs. mg CH₄-C m-2 s-1), and we recommend checking sensitivity to the choice of theta or using arcsinh(x/SD(x)) --- normalizing by the standard deviation before transforming --- to make the scaling data-driven and unit-independent. For datasets with a very high proportion of near-zero or zero values, a two-part (hurdle) model that explicitly separates the probability of a nonzero flux from the distribution of nonzero flux magnitudes may be more appropriate than any single transformation; this approach naturally accommodates mixed-sign distributions but requires sufficient sample size in both components.

## 6.3 Empirical comparison

For the wetland soil dataset, the raw distribution is strongly left-skewed (Shapiro-Wilk W = 0.15, p \< 10⁻³³). The log₁₀ transform produces a roughly normal distribution (W = 0.97, p = 0.33) but retains only 17% of the data. The arcsinh transform retains all 288 values and substantially improves symmetry (W = 0.83, p \< 10⁻¹⁶) --- not perfectly normal, but far more tractable for parametric methods than the raw data, and without discarding the majority of observations.

Across the three datasets, raw variance spans a 1,334-fold range (wetland soil: 139, tree canopy: 0.57, tree stem: 0.10 \[nmol m⁻² s⁻¹\]²). After arcsinh transformation, the variance ratio drops to 29-fold (wetland soil: 1.06, tree canopy: 0.23, tree stem: 0.04). This improvement in variance homogeneity is critical for mixed models or ANOVAs that compare flux rates across ecosystem components.

**[[ Figure 12 reference: Q-Q plots under three transformations for wetland soil CH4 (already produced). Figure 11 reference: arcsinh transformation behavior plot (already produced). ]]**

## 6.4 Practical recommendations for statistical analysis

We recommend the following approach for statistical analysis of chamber CH₄ flux data:

-   **Arcsinh for visualization and variance stabilization.** Use arcsinh-scaled axes for plotting and as a transformation for exploratory analysis and model fitting (e.g., to improve residual behavior in linear mixed models). It is preferable to log because it does not exclude data.

-   **Non-parametric methods for group comparisons.** Wilcoxon rank-sum tests for two-group comparisons, Kruskal-Wallis for multiple groups, and permutation tests for complex designs. These use all observations regardless of sign or magnitude and are robust to heavy-tailed, mixed-sign distributions.

-   **Bootstrapped confidence intervals for means used in scaling.** Bootstrap resampling (e.g., BCa intervals) provides valid CIs without requiring normality, and is particularly appropriate for mean flux estimates used in budget calculations.

-   **GLMs/GLMMs for regression.** Gaussian GLMs with an identity link on arcsinh-transformed data work well when the below-MDF fraction is modest, because instrument noise scatters symmetrically and residuals are approximately homogeneous. When a large fraction of measurements fall below detection, the appropriate model is *censored regression* (Tobit-style), treating each below-MDF measurement as left-censored at its per-measurement detection limit rather than as an observed value or a zero. This is the statistically correct treatment: it uses all the data and models the censoring explicitly without substituting a value. Note that setting below-MDF values to zero is *not* recommended for regression: it selectively truncates the lower tail of the error distribution for near-detection observations only, creating heterogeneous residuals across the dataset that violate OLS assumptions in a way that is difficult to diagnose or correct. For mixed-sign flux data with moderate below-MDF fractions, two-part (hurdle) models separating the probability of a detectable flux from the distribution of detectable flux magnitudes are an alternative worth considering.

-   **If log-transforming, report data loss explicitly.** State the number and percentage of values excluded, their distribution (noise-scattered near-zero or true biological uptake), and whether the remaining subset is representative of the population.

# 7. Synthesis and recommendations

The analyses above document a processing pipeline in which each choice feeds into the next and errors compound downstream. Relying on manufacturer precision specifications inflates detection rates by underestimating the noise floor of noisier instruments. Conservative MDF methods, when combined with removal of below-detection values, can inflate reported mean fluxes by more than 100% relative to the true population mean. Log transformation then discards the near-zero and negative values that survived both previous steps, compounding the upward bias further. For small-flux ecosystem components --- leaves, dormant-season stems, low-emitting species --- the cumulative effect of these three choices can render reported flux estimates nearly unrecognizable relative to the underlying data.

The key unifying principle is that the appropriate processing choice depends on the scientific question. For population estimates and budget calculations, preserving all measurements is essential --- the noise in below-detection values is informative because it scatters symmetrically around the true flux. For detection and sign questions, filtering to the confidently-detected subset is appropriate because the scientific question is specifically about the resolved signal. These two objectives are genuinely different and require genuinely different treatments; conflating them by applying detection-oriented filtering to population-estimation problems is the root cause of most of the biases documented here.

A further dimension of this problem --- not fully captured by noise-based framing --- is that near-zero fluxes arise from two distinct sources: instrument noise scattering around a small true flux, and genuine biological suppression of flux under cold, dry, or dormant conditions. Quality filters based on R², SNR, or MDF cannot distinguish between these sources, yet their consequences for filtered datasets differ. Noise-driven near-zero values scatter symmetrically and their mean converges on the true flux when retained; retaining them is therefore unbiased. Biologically real near-zero values represent a specific ecological state --- winter soil, drought-stressed plant, dormant stem --- that is as scientifically meaningful as any high-flux measurement. Excluding them is not a noise correction; it is a structured removal of an entire condition from the dataset, with direct and predictable consequences for any analysis that spans seasons, years, or environmental gradients. The practical implication is that quality filtering criteria developed and evaluated on warm-season or high-flux datasets --- where near-zero values are rare and mostly noise-driven --- may perform very differently when applied to year-round or climatically diverse datasets where genuine near-zero states are common.

**[[ REPORTING CHECKLIST to add here (or as a numbered list / small table): minimum set of statistics to report in any chamber CH4 methods section --- instrument model; empirical precision estimate (method + value); n total measurements; fraction below MDF (method used); filter applied; below-MDF treatment; transformation used; outlier criteria if any. This enables cross-study comparison even when methodologies differ. ]]**

# 8. Conclusions

Chamber CH₄ flux data processing involves a sequence of choices --- precision estimation, detection limit calculation, quality filtering, below-detection treatment, and data transformation --- each of which introduces predictable, quantifiable bias if handled carelessly. Manufacturer precision specifications are unreliable substitutes for empirical estimates and can underestimate real instrument noise by up to fourfold or more. The choice of detection threshold matters less than what is done with the measurements that fall below it: removing below-detection values inflates reported means by 17--135% and medians by 40--319%, while retaining original values with a detection flag preserves both statistics. The log transformation, near-universally applied for variance stabilization, silently discards 10--83% of measurements depending on the prevalence of negative fluxes; the arcsinh transformation is a direct replacement that handles negative values, stabilizes variance, and introduces no data loss.

These processing choices are not independent: they interact, and the appropriate choice at each step depends on the scientific question being asked. We recommend a question-dependent framework in which population estimates and budget calculations use all retained measurements with detection status flagged, while detection and sign questions filter to the confidently-detected subset. Regardless of the approach chosen, we recommend that authors report empirical instrument precision, the fraction of measurements below detection, the below-MDF treatment applied, and the transformation used --- sufficient information to allow readers to assess comparability across studies and, where necessary, reprocess data under alternative assumptions.

# Data and code availability

**[[ Data availability statement --- datasets from Harvard Forest and Yale-Myers Forest; to be deposited in \[repository\] upon acceptance. Code for MDF calculations and arcsinh comparison available at \[repository / goFlux package citation\]. ]]**

# Acknowledgments

**[[ Acknowledgments to be added. ]]**

# References

**[[ Full reference list to be added. Key citations needed: goFlux package; Wassmann et al. 2018; Christiansen et al. 2015; Whittaker & Woodwell 1967; Helsel 2005 (non-detects in environmental data); instrument manufacturer specs for GLA131 and LI-7810. ]]**
