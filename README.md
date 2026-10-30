# Oxygen Reduction — Joint Model


## Recharge terminology

Throughout this analysis, **infiltration** refers specifically to the **land-surface recharge (LSR)** pathway (soil and unsaturated-zone percolation). The distinct **river-recharge** pathway is water derived from the river. Do not describe oxygen loss associated with river recharge as "infiltration loss."

- **River DO (monitoring reference):** measured river-water DO, represented centrally by 7.5 mg/L.
- **River-recharge DO at aquifer entry:** effective DO of the river-derived component entering the regional groundwater system; not directly measured.
- **River-recharge DO deficit:** monitored-river DO minus river-recharge DO at aquifer entry; a scenario parameter **not a measured riverbed reaction loss**. It could reflect riverbed passage and other source/boundary-condition differences.
- **Effective LSR DO:** land-surface-derived recharge DO at aquifer entry after any soil/unsaturated-zone consumption.
- **Apparent age-dependent reduction:** DO depletion associated with groundwater MRT; it does not uniquely locate oxygen consumption in the saturated aquifer.

Script 07 retains its original `07_infiltration_loss_sensitivity.R` filename and `infiltration_*.csv/png` output filenames **only for backwards compatibility**. Plot titles, model labels, tabular variables and scientific interpretation use the corrected terminology.


A standalone two-script workflow for the 30 October 2026 Heretaunga groundwater analysis. It is independent of the earlier **fixed-rainfall-DO** rate calculations.

## Files

- `scripts/01_three_tracer_mixing.R` — source fractions from δ¹⁸O, chloride and excess air; three-tracer agreement screen and rainfall-endmember sensitivity.
- `scripts/02_joint_oxygen_partition.R` — simultaneous fitting of source-associated oxygen deficit and aquifer reduction, comparing first- and zero-order models using leave-one-well-out prediction and a conditional bootstrap.
- `data/groundwater.csv` — **not committed**. The source research repository is private and this repository is public; use the groundwater.csv provided in the separately prepared workflow ZIP.

## Run in Positron

Install `tidyverse` if needed (`install.packages('tidyverse')`). Open this repository as your working directory, place the prepared dataset at `data/groundwater.csv`, then run:

```r
source("scripts/01_three_tracer_mixing.R")
source("scripts/02_joint_oxygen_partition.R")
```

Results are written to `outputs/` and figures to `plots/`. These directories are Git-ignored.

## Assumptions

- Mixing uses site 2801 as the river **tracer endpoint**, and 1032, 4294 and the 108 m HB Prison sample as the land-surface **tracer endpoint**.
- For each well, the median of the clipped rainfall fractions from δ¹⁸O, chloride and excess air is used. Reliability requires all three tracers and a range of **unclipped** fractions ≤ 0.35. No minimum river recharge fraction.
- 29 of 42 records pass the production screen (five rainfall-dominated).
- Monitored river DO reference = 7.5 mg/L; **assumed** pre-infiltration *land-surface recharge* DO reference = 9.5 mg/L.
- For the first-order model:
  
  DO_pred = [fR * 7.5 + fP * (9.5 - L)] * exp(-k * MRT)

  L and k are fitted together. For zero-order reduction, the post-mixing concentration is reduced by k * MRT, bounded below by zero.

Only **effective** rainfall-source DO (9.5 - L) is constrained by the data. An inferred land-surface source-associated deficit cannot, by itself, be assigned specifically to soil or unsaturated-zone oxygen consumption. The first-order rate is an apparent residence-time-dependent rate, not independent proof of saturated-zone reaction.

Excess air is inferred from dissolved gases; its treatment as a conservative mixing tracer is an **assumption needing investigation**. This analysis does not model excess N₂ or denitrification kinetics.

## Verification

The prepared standalone model and selection were numerically cross-checked against the earlier analysis. R itself was not executed in the preparation environment; please run both scripts in Positron to validate the outputs.

## Bayesian uncertainty test (optional; run after script 01)

The main scientific question is whether groundwater residence time explains oxygen depletion when the DO entering the aquifer from each recharge source is uncertain. The Bayesian script compares **recharge-only (k = 0)** against **recharge plus a first-order age-dependent depletion coefficient**.

Run in Positron after creating the local groundwater.csv:

```r
source("scripts/01_three_tracer_mixing.R")
source("scripts/04_bayesian_recharge_DO.R")
```

This Bayesian script is a **separate sensitivity experiment**: it does not overwrite the joint-model estimates in script 02. Script 03 is an optional plotting script shared separately.

**Prior specification and interpretation**

- Population mean river DO: Normal(7.5, 0.7 mg/L), truncated to (0, 12). The 0.7 mg/L uncertainty of the mean is an **assumption**, not the observed river temporal SD.
- Between-sample river DO SD: half-normal(scale 1.9 mg/L), informed by river DO variability.
- Population mean effective rainfall DO: uniform(0, 9.5 mg/L).
- Between-sample rainfall DO SD: half-normal(scale 2.5 mg/L).
- Additional unexplained site variation: half-normal(scale 1.5 mg/L). The measurement SD is fixed at 0.2 mg/L.
- First-order depletion k: half-normal(scale 0.10 year^-1), or exactly zero in the competing model.

The between-sample source DO values are Gaussian random effects that are *integrated out analytically*. This is a modelling approximation: individual latent source DO is not explicitly truncated at physical bounds.

In addition to the baseline priors, a wider prior scenario tests k sensitivity to the assumptions. The script produces posterior summaries, a simple four-chain R-hat check, **matched five-fold cross-validated** predictive RMSE and log score, and diagnostic PNG plots. The cross-validation refits both models using only the training folds.

**Scientific limitation:** These models assume initial recharge DO variability is statistically independent of groundwater mean residence time, conditional on mixing fraction. A positive k establishes a need for *age-associated* depletion under those assumptions. It does **not** establish the site of oxygen consumption (recharge pathways versus saturated-aquifer transport). The three-tracer recharge fractions are held fixed; excess-air tracer suitability and source fraction uncertainty are not propagated here.

The initial posterior and predictive comparison was cross-checked using a separate Python implementation. **The R script has not been executed in this environment** and should be run in Positron.


## Confounding test: age versus recharge-pathway differences

Run after script 01:

```r
source("scripts/01_three_tracer_mixing.R")
source("scripts/05_recharge_pathway_confounding.R")
```

Script 05 fits four models, using the same hierarchical DO priors as script 04:

1. Source DO uncertainty alone, no age term.
2. Separate mean effective source DO for confined and unconfined wells, no age term.
3. Uncertain source DO plus first-order age term.
4. Both confinement-related source DO differences and first-order age term.

These *confinement-related differences* are a proxy sensitivity analysis, **not** proof that the recharge pathway is different or that depletion took place before water reached the aquifer.

Confinement and age are correlated; 29 wells include related Eastbourne, Gilbertson and Flaxmere observations. Therefore script 05 uses **grouped five-fold CV**, keeping these obvious groups together. Grouping by name is a heuristic, not a spatial or temporal independence guarantee.

Outputs:
- outputs/confounding_posterior_summary.csv
- outputs/confounding_group_CV_predictions.csv
- outputs/confounding_group_CV_summary.csv
- outputs/confounding_within_group_correlations.csv
- plots/confounding_DO_MRT_by_confinement.png
- plots/confounding_group_CV_DO.png

**Structural non-identifiability:** An alternative model in which river- and rainfall-recharge DO at entry vary as `R_i = R exp(-b MRT_i)` and `P_i = P exp(-b MRT_i)` produces *exactly the same predictions* as a first-order transport attenuation model with k=b. Fitted DO–MRT observations alone cannot locate oxygen depletion along recharge pathways versus subsequent saturated-aquifer transport.

The confinement-group analysis was calculated independently in Python for checking. The R script must still be run and diagnosed in Positron. The Bayesian random effects are Gaussian and can imply unphysical negative inlet DO at extreme tails; physical constraints and independent recharge DO observations would improve the next model.


## Age-adjusted river DO prior (script 06)

Run from Positron after the mixing script:

```r
source("scripts/01_three_tracer_mixing.R")
source("scripts/06_age_adjusted_river_DO_prior.R")
```

The model holds the river DO **mean** at an uncertain reference of 7.5 mg/L and permits age-related attenuation of **variance**. The reported 1.9 mg/L SD describes the monitored **daily means**, not the diurnal cycle. The daily-mean variance is partitioned between a persistent component and a seasonal sinusoid; the analysis tests 100%, 50% and 25% persistent variance. An additional 1.0 mg/L **diurnal input SD is assumed, not measured**, and is attenuated separately. The attenuation filter for either oscillation is `H(T,MRT) = 1/sqrt(1+(2*pi*MRT/T)^2)`, derived from an **assumed exponential travel-time distribution**.

The same six model variants are fitted and compared using matched five-fold cross-validation: each of the three river-prior variance scenarios with either k=0 or a fitted first-order age-related coefficient. Outputs are unique to this script:
- `outputs/age_adjusted_river_prior_posteriors.csv`
- `outputs/age_adjusted_river_prior_CV_summary.csv`
- `outputs/age_adjusted_river_prior_CV_predictions.csv`
- `outputs/age_adjusted_river_prior_convergence.csv`
- `plots/age_adjusted_river_prior_k.png`
- `plots/age_adjusted_river_prior_CV.png`

**Interpret carefully:** Only site 2801 has MRT=0 in the current 29-well selection; the remaining MRTs are ≥1 year. The year-specific attenuation is an assumed travel-time response; MRT alone does not determine dispersion or attenuation. Piston flow would simply delay the sinusoid. The daily-mean variance split cannot be estimated from the 29 well DO measurements, and should be constrained by time-frequency analysis of the actual river record where possible. The age term still cannot distinguish oxygen depletion along recharge pathways from subsequent saturated-aquifer reaction. The code is supplied for execution in Positron; the equivalent numerical sensitivity was independently evaluated in Python, but this R script has not been run here.


## River-derived recharge DO deficit: sensitivity and lower bounds (script 07)

Run after the three-tracer mixing model:

```r
source("scripts/01_three_tracer_mixing.R")
source("scripts/07_infiltration_loss_sensitivity.R")
```

Script 07 tests the remaining major uncertainty in the apparent oxygen reduction rate: the river's monitoring DO (reference 7.5 mg/L) is **not** necessarily the DO entering the saturated aquifer after passage through riverbed and connected sediments.

- Assume a river-derived recharge DO deficit of 0–7.5 mg/L in 0.25 mg/L increments relative to monitored river DO; for each assumed deficit refit effective land-surface recharge DO (0–9.5 mg/L) and first-order age coefficient k. Repeat all fits **with and without** young, high-DO site 2801.
- Use grouped five-fold CV (related Eastbourne, Gilbertson and Flaxmere wells together), reporting within-sample and held-out RMSE separately.
- Calculate the **minimum river-derived recharge DO deficit** required for each well if **no age-dependent loss occurs after recharge entry**, allowing the land-surface recharge component to be entirely anoxic. This favours the no-aquifer-reduction explanation; oxygen remaining in land-surface recharge would require a larger deficit in the river component.
- Repeat required-loss calculation for river references 5.6, 7.5 and 9.4 mg/L. Flag samples, notably site 2801 at 7.5 mg/L, whose DO exceeds the maximum implied by these fixed source bounds. A river sample can temporarily exceed its long-run mean; this is not evidence of oxygen production in the aquifer.
- Compare physically bounded k=0 models with shared effective source DO, with different river DO for confined versus unconfined wells, or with **both river and rainfall effective DO varying by confinement**; contrast these against an age-dependent model with river mean fixed at 7.5 mg/L.

Output tables:
- `outputs/infiltration_loss_profile.csv`
- `outputs/infiltration_minimum_loss_per_well.csv`
- `outputs/infiltration_required_loss_summary.csv`
- `outputs/infiltration_model_CV_predictions.csv`
- `outputs/infiltration_model_CV_summary.csv`

Output figures:
- `plots/infiltration_profile_CV.png`
- `plots/infiltration_profile_k.png`
- `plots/infiltration_minimum_loss_by_age.png`

**Inference limit:** The minimum-deficit calculation identifies what a no-age-dependent-reduction explanation *requires*, not where oxygen consumption occurred. The confinement model is only a crude proxy for differences in recharge entry conditions; no independent recharge-entry DO measurements are available. If initial recharge DO is permitted to vary systematically with MRT, a sufficiently flexible recharge-entry-only model can be exactly equivalent to the first-order age model. The independent Python calculations provide numerical results; run the R script in Positron to verify its outputs. Raw groundwater data have not been committed.


## Main comparison: aquifer-age effect vs LSR oxygen loss (script 08)

**This is the preferred analysis of the central scientific question.** Script 07 explored river-recharge entry deficits, which are not the present focus. "Infiltration" refers only to the land-surface recharge pathway; **no riverbed loss is estimated in script 08**.

Run in Positron after script 01:

```r
source("scripts/01_three_tracer_mixing.R")
source("scripts/08_LSR_vs_aquifer_oxygen_loss.R")
```

The same 29 selected wells and three-tracer source fractions are used, with a repeat analysis omitting site 2801. Four models are compared:

| Model | First-order groundwater-age effect | LSR oxygen loss |
|---|---|---|
| A | Estimated | Fixed at zero (effective LSR DO = 9.5 mg/L) |
| B | Estimated | One fitted common LSR loss |
| C | Estimated | Between-well variable LSR loss, hierarchically distributed |
| D | Fixed k=0 | Between-well variable LSR loss as in C |

LSR effective DO is restricted to **0–9.5 mg/L** in B/C/D. In C/D, site-specific LSR DO values come from a bounded logistic-normal source distribution, **integrated out**, rather than fitting one free loss parameter per observation. Model A/B give a fixed LSR DO for all wells. This tests whether the variable LSR hypothesis improves predictions without making the model automatically fit every well.

The river monitoring reference is 7.5 mg/L. Its population mean has an **assumed** N(7.5, 0.7) prior. Source DO temporal variance is based on the ~1.9 mg/L SD of **daily means**, with age-dependent attenuation under an **assumed exponential groundwater travel-time distribution**, 25% assumed persistent variance and an additional assumed 1 mg/L diurnal SD. The true distribution and variance split are not independently known. Shared residual model-discrepancy SD has a half-normal scale of 1.5 mg/L; DO measurement uncertainty is assumed 0.2 mg/L. The fitted apparent rate k has a half-normal scale of 0.10 year⁻¹. Model C/D LSR spatial variability is governed by a half-normal prior on its logistic SD (scale 1.5). These are **analytical assumptions**, not measured variances or reaction rates.

Script 08 fits four independently initialized Bayesian MCMC chains for each model, provides basic R-hat diagnostics, and compares posterior predictive distributions with grouped five-fold cross-validation (including obvious Eastbourne, Gilbertson and Flaxmere groups). The script generates:

- `outputs/LSR_loss_posterior_summary.csv`
- `outputs/LSR_loss_convergence.csv`
- `outputs/LSR_loss_CV_predictions.csv`
- `outputs/LSR_loss_CV_summary.csv`
- `plots/LSR_loss_posterior_k.png`
- `plots/LSR_loss_CV_predictions.png`
- `plots/LSR_loss_CV_residuals.png`

The model was also evaluated independently in Python; approximate five-fold Bayesian posterior predictive RMSEs using all 29 wells were 1.08 (A), 1.11 (B), 1.12 (C bounded), and 3.09 (D bounded) mg/L. These are **independent numerical checks** and should not be misrepresented as results of the R code; repeat the script in Positron to validate R performance. Small differences due to MCMC randomness are expected.

**Key inference:** Adding mean or variable LSR oxygen loss does not improve held-out prediction under these assumptions, and inferred k remains around 0.07 year⁻¹. The LSR loss parameters are poorly constrained: only five wells are LSR-dominated. The fit does **not** uniquely identify where the age-related depletion occurred; age-correlated source/boundary variation and uncertainty in mixing fractions or MRT remain relevant alternatives. The model uses fixed source fractions and MRT point estimates. The logistic-normal quadrature is approximate. No private groundwater data are stored in this public repository.


## DO residuals versus LSR contribution and mixing uncertainty (script 09)

This is a diagnostic of the EXISTING first-order model, not a new oxygen-loss process.

```r
source("scripts/01_three_tracer_mixing.R")
source("scripts/09_LSR_mixing_residual_diagnostics.R")
```

- Refit the current first-order model (monitored river DO fixed at 7.5 mg/L; common **effective LSR DO** bounded 0–9.5 mg/L); predict each held-out well by **leave-one-well-out** refitting.
- Test whether absolute/signed DO residuals correlate with LSR fraction or with **raw (unclipped)** disagreement among the three tracer-derived LSR fractions. Repeat excluding site 2801. Spearman p-values are exploratory.
- Hold river DO and the fitted apparent k fixed, then permit **each well's LSR entry DO to vary independently across the entire physical 0–9.5 mg/L range**. Compute whether observed DO falls within that well's permitted prediction interval. This is an extreme sensitivity bound, not a fitted hierarchical model.
- Refitting the same 29 selected wells with δ18O alone, chloride alone, and excess air alone checks whether alternative tracer weights explain DO residuals. Report both apparent k and inferred effective LSR DO: a stable k **does not imply** inferred LSR oxygen loss is stable.
- Compare DO prediction ranges obtained from the three tracer-specific fractions, without changing other model parameters.
- Optional spatial Moran's I diagnostic with five symmetric nearest neighbors and permutations. The **private** sample coordinates must be present locally, as `E/N` columns in the mixing file or in a local `data/Uwe_HP_data_prepared.csv`. This is exploratory with clustered samples. No coordinates are added to GitHub.

Outputs are `outputs/LSR_mixing_residuals_by_well.csv`, `outputs/LSR_mixing_residual_correlations.csv`, `outputs/LSR_mixing_alternative_tracers.csv` and optionally `outputs/LSR_mixing_residual_spatial_test.csv`. Plots are `plots/LSR_mixing_abs_residual_vs_LSR.png`, `plots/LSR_mixing_abs_residual_vs_tracer_disagreement.png` and `plots/LSR_mixing_possible_DO_by_well.png`.

**Independent Python diagnostic of the 29-well dataset:** Spearman correlations between absolute LOOWO DO error and LSR fraction -0.19 (p=0.32), and with intertracer disagreement -0.13 (p=0.51). Under the fitted fixed river DO and k, independently variable LSR entry DO within 0–9.5 mg/L can match only 6 of 29 observations. Refitting using individual tracers gives k about 0.0656–0.0716 yr^-1, but effective LSR DO 6.74–9.17 mg/L. This demonstrates uncertainty in the LSR DO-loss partition even though changing mixing fractions does not improve predictive DO RMSE. None of these results proves that the first-order rate belongs solely to saturated-aquifer reaction. The R script needs its first execution in Positron.

## Groundwater age-distribution sensitivity (script 10)

This tests the limitation of assigning one water age to every well based on its reported mean residence time (MRT). **We must first establish how MRT was estimated from the original age tracers before any particular age-distribution shape is interpreted physically.**

Run in Positron, in sequence:

    source("scripts/01_three_tracer_mixing.R")
    source("scripts/10_groundwater_age_distribution_sensitivity.R")

For the piston model the DO attenuation is exp(-k*MRT). For a gamma groundwater age distribution with **fixed** shape alpha and mean MRT, the integrated attenuation is (1+k*MRT/alpha)^(-alpha). Shape alpha=1 is an exponential age distribution; alpha approaching infinity is the piston limit. The coefficient of variation of groundwater ages is 1/sqrt(alpha).

Script 10 compares alpha = Inf, 100, 50, 20, 10, 5, 3, 2, 1. For each shape it refits effective LSR entry DO within 0–9.5 mg/L and a common apparent first-order coefficient k. The river reference is held at 7.5 mg/L; no river-recharge deficit is estimated. The test uses **bounded least squares without new Bayesian priors**. Source fractions, reported MRT and river DO are treated as fixed. The same grouped five-fold CV splits are applied across shapes, and all fits repeat without site 2801.

Files generated:
- outputs/groundwater_age_distribution_sensitivity.csv
- outputs/groundwater_age_distribution_CV_predictions.csv
- plots/groundwater_age_distribution_CV.png
- plots/groundwater_age_distribution_k.png
- plots/groundwater_age_distribution_residuals.png

**Independent Python check (all 29 wells, grouped CV):** Piston k=0.0701 yr^-1, RMSE 1.001 mg/L; gamma shape 10 k=0.0761, RMSE 1.005; shape 5 k=0.0824, RMSE 1.014; shape 2 k=0.1038, RMSE 1.065; exponential shape 1 k=0.1491, RMSE 1.192. **The Python and R randomized group folds are not identical**, so do not present the Python numbers as output from executing the R script.

**Interpretation:** Wider age distributions increase the fitted apparent reduction coefficient but do not improve predictive DO accuracy. This does not establish that the age distribution is piston flow or that oxygen consumption occurs exclusively in the saturated aquifer. The river-derived and LSR components may have different age distributions, and the reported MRT could already be conditional on a specific lumped-parameter model. Age-tracer dating details and uncertainty are needed before introducing Bayesian priors on age-distribution shape.

## Bayesian inference for age-distribution and oxygen-loss parameters (script 11)

**Preferred next step after fixed gamma-shape sensitivity in script 10:** Put priors on the uncertain **gamma shape** alpha and apparent first-order coefficient k, and fit them jointly with effective land-surface recharge (LSR) DO, mean river DO and residual variation. This does **not** investigate how the supplied groundwater MRT values were originally derived, nor does it fit a river-recharge loss term.

Run after the three-tracer mixing script:

    source("scripts/01_three_tracer_mixing.R")
    source("scripts/11_bayesian_age_parameter_priors.R")

The groundwater-age attenuation is **(1 + k*MRT/alpha)^(-alpha)** for gamma shape alpha, with piston-flow exp(-k*MRT) fitted as a separate Bayesian reference model. The R code uses 29 three-tracer-screened wells. Predictions under the gamma model are conditional on the given MRTs and the assumption that river and LSR fractions share the same age distribution.

**Specified priors (some assumptions, not independently measured):**
- Mean river DO: Normal(7.5, 0.7 mg/L), bounded to 0–12 mg/L. The prior SD of the *mean* is assumed; it is **not** the measured 1.9 mg/L SD of daily means.
- Effective LSR DO at aquifer entry: Uniform(0, 9.5 mg/L), representing an unknown common LSR oxygen deficit but **not** an independent LSR-loss parameter per well.
- Apparent first-order k: HalfNormal(SD = 0.2 yr^-1).
- Residual model-discrepancy SD: HalfNormal(SD = 1.5 mg/L), with an additional assumed measurement SD = 0.2 mg/L.
- Gamma shape: log(alpha) ~ Normal(log(median), 1.3), alpha restricted to 0.25–500. **Three median scenarios 2, 10 and 50** measure sensitivity to this consequential prior.

The script runs four independent MCMC chains per full-data scenario and a basic split-Rhat check. It repeats the reference piston and median-10 gamma posteriors without site 2801, and runs grouped five-fold held-out comparisons for piston and all three gamma priors (using two shorter chains per fold). Output tables are:
- outputs/bayesian_age_parameter_posterior_summary.csv
- outputs/bayesian_age_parameter_convergence.csv
- outputs/bayesian_age_parameter_correlations.csv
- outputs/bayesian_age_parameter_CV_predictions.csv
- outputs/bayesian_age_parameter_CV_summary.csv

Figures are:
- plots/bayesian_age_parameter_k_intervals.png
- plots/bayesian_age_parameter_alpha_prior_sensitivity.png
- plots/bayesian_age_parameter_CV.png

**Independent Python numerical test, not yet executed as R:** Gamma alpha posterior medians are approximately 5.0, 12.8, and 47.3 for alpha prior medians 2, 10, and 50. Corresponding k posterior medians are approximately 0.083, 0.077, and 0.072 yr^-1. The piston-flow Bayesian k posterior median is about 0.070 yr^-1. Under the median-10 prior, k 95% CrI is about 0.059–0.106 and alpha about 2.3–128. Basic split-Rhat values for the independent Python chains were <1.01. The independent Python grouped CV had RMSE 1.173 mg/L (piston) versus 1.189 mg/L (gamma, prior median 10), a **small difference**; the R implementation uses its own randomized fold arrangement, so compare its generated results rather than treating the Python scores as R outputs.

**Conclusion:** alpha and k are negatively correlated. DO and MRT data poorly identify groundwater age-distribution width without independent constraints: alpha posteriors shift considerably with alpha priors. k is more stable but its uncertainty broadens when alpha is allowed to vary. This analysis does **not** establish the physical location of oxygen depletion or account for uncertainty in mixing fractions and MRT. No groundwater samples are uploaded to the public repository.
