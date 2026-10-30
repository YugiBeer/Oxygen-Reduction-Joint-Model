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


## Joint Bayesian source-mixing, LSR DO and MRT model (script 12)

**This addresses the main joint-inference gap.** Earlier Bayesian oxygen-reduction scripts held the three-tracer source fractions and groundwater MRT fixed. Script 12 estimates their uncertainties **simultaneously** with the oxygen model, fitting the original measured tracer values rather than the previously calculated mixing fractions.

**Conservative mixing tracers:** chloride, δ18O and **excess air (xs_Air)**. **Excess N2 (xs_N2) is not treated as a conservative mixing tracer** because it can include denitrification products. It may be included as a fourth, *reactive* observation via an optional empirical age-related association; this is **not** a stoichiometric denitrification model and must not be interpreted causally without an independent nitrogen budget.

The Stan implementation is in `stan/12_joint_bayesian_mixing.stan`; the Positron runner is in `scripts/12_joint_bayesian_mixing.R`. Requires R packages `tidyverse`, `here`, `cmdstanr`, `posterior`, and a separate **CmdStan installation**. For example, if `cmdstanr` is not installed:

```r
install.packages("cmdstanr",
                 repos = c("https://stan-dev.r-universe.dev",
                           getOption("repos")))
cmdstanr::install_cmdstan()
```

Ensure the **local, private** file `data/groundwater.csv` exists (created by the earlier local workflow; never commit it). Then run:

```r
source("scripts/12_joint_bayesian_mixing.R")
```

The default fit uses **41 wells with nonmissing MRT and DO**, including groundwater samples previously excluded for disagreement among tracer-derived fractions. A missing δ18O or chloride measurement is skipped in the corresponding observation likelihood; it does not force the well out of the analysis. Each sample has a latent bounded LSR fraction with a Uniform(0,1) prior. Estimated tracer end members each have priors and observational errors, and a Student-t(4) tracer likelihood with fitted extra dispersion limits the influence of samples inconsistent with simple binary mixing.

**Default priors and sensitivities:**
- River δ18O: Normal(−7.79, 0.24‰); LSR δ18O: Normal(−6.58, 0.40‰).
- River Cl: Normal(4.9, 3 mg/L); LSR Cl: Normal(60, 15 mg/L), with physically bounded nonnegative values.
- River excess air: Normal(0.122, 0.35 mL/kg); LSR excess air: Normal(4.04, 1.4 mL/kg), physically bounded nonnegative values.
- Monitored river DO mean: Normal(7.5, 0.8 mg/L), bounded 0–12 mg/L.
- Effective LSR DO at aquifer entry: Uniform(0,9.5 mg/L), so its deficit relative to the assumed 9.5 mg/L pre-infiltration reference is estimated.
- Apparent first-order reduction k: HalfNormal(0.15 yr^-1); extra DO residual SD: HalfNormal(1.5 mg/L).
- Individual true MRT: log-normal with **expected value equal to the reported MRT** and assumed log SD 0.25 or 0.50. The log-location parameter is corrected by subtracting half the squared log SD. A reported zero MRT receives a near-zero half-normal prior.
- Optional spatial variation in LSR DO via a hierarchical logit-normal distribution (`cfg$vary_lsr_do = TRUE`); default is a common effective LSR DO across wells.
- Optional gamma age-distribution shape with log-normal prior median 10, log SD 1.3 (`cfg$default_gamma = TRUE`); default is piston attenuation using the latent MRT.
- Optional xs_N2: prior on river/LSR background excess N2, an empirical nonnegative coefficient relating N2 to apparent age-associated DO loss, and extra residual variance (`cfg$use_n2` scenario in R). **Do not interpret this coefficient as nitrogen production stoichiometry.**

**Important empirical-Bayes caveat:** The endmember prior centres derive from the same study samples used in the fit, so the prior and likelihood are not fully independent. This can effectively double-use evidence and narrow source-fraction posteriors. Test a larger `cfg$tracer_prior_scale`, anchor with external measurements, and do not treat the endmember posteriors as independent evidence of source composition.

With `run_prior_sensitivity = FALSE`, the runner fits three scenarios: MRT log SD 0.25 without N2, MRT log SD 0.25 with N2 as reactive response, and MRT log SD 0.50 without N2. With **`run_prior_sensitivity = TRUE` (currently the script default)**, it adds four more Stan fits: doubled end-member prior SDs, Beta(2,5) oxygen-poor LSR prior, Beta(5,2) oxygen-rich LSR prior, and a hierarchical varying-LSR-DO model. Seven four-chain Stan fits can take substantial time. It saves posterior summaries, individual source-fraction intervals, latent MRT intervals, fitted DO, MCMC diagnostics, plot and serialized fits under gitignored `outputs/` and `plots/`. The detailed well-level outputs remain local.

**Independent numerical prototype** (Python custom-MH integrating over latent fLSR and true MRT, NOT execution of Stan/R code): under MRT log SD 0.25, apparent k approximately 0.072 (95% CrI 0.049–0.104) yr^-1 and common LSR loss approximately 2.44 (0.15–7.04) mg/L. Adding reactive excess N2 gives k approximately 0.071 (0.047–0.101). MRT log SD 0.50 without N2 gives k approximately 0.085 (0.056–0.121). These runs use a different numerical implementation (quadrature plus Metropolis rather than explicit latent-variable HMC); their exact posteriors and credible intervals need not match the Stan fit. The Python fit also indicates that inferred site LSR fractions can differ appreciably from three-tracer medians, including wells 63, 1247 and 3698 Flaxmere. The R/Stan model **has not yet been compiled or run in the user's environment**.

**Identifiability limit:** A broad posterior on LSR source DO and source fractions means the partition between oxygen depleted during LSR and oxygen depleted in subsequent groundwater transport is not fixed by DO observations alone. The inferred k remains a model-conditional **age-associated** coefficient, not direct proof that all oxygen reduction happened within the saturated aquifer. Need posterior predictive checks, Rhat/ESS assessment, endmember-prior sensitivity and matched group CV before treating any of these outputs as a definitive result.


## Posterior parameter trade-off, prior sensitivity, and modelled oxygen loss (script 13)

After the full joint Bayesian Stan analysis in script 12 has completed, run:

```r
source("scripts/13_joint_posterior_partition_diagnostics.R")
```

Script 13 reads **local** fitted Stan RDS files under `outputs/` and calculates:
- Joint posterior Pearson and Spearman correlations between LSR DO deficit and apparent k, river DO and k, and overall tracer-inferred LSR fraction and k.
- Per-well posterior distributions of model-implied **LSR-pathway oxygen loss** `fLSR*(9.5 - LSR_entry_DO)` and **apparent age-associated oxygen depletion** `entry_DO - mean_DO`. Neither term uniquely identifies the physical location of consumption: the latter is not necessarily saturated-aquifer loss.
- Arithmetic summaries across the sampled wells, **not volumetric aquifer budgets**. The source-reference LSR DO 9.5 mg/L is an assumption.
- R-hat and effective sample size for global parameters, with warnings for R-hat above 1.05 or bulk ESS under 400.
- Joint posterior plots of LSR DO loss vs k, river DO vs k, and interval plots of modelled oxygen loss components under different priors.

Output tables: `outputs/joint_bayes_partition_correlations.csv`, `outputs/joint_bayes_partition_summary.csv`, `outputs/joint_bayes_partition_by_well.csv`, `outputs/joint_bayes_partition_convergence.csv`. These and all posterior draws are gitignored. Output figures use the `plots/joint_bayes_partition_*.png` prefix.

**Independent Python diagnostics using 41 wells, no reactive N2, and mean-centred MRT log-normal prior (log SD 0.25):** A long four-chain fit with the Uniform(0,9.5) LSR entry DO prior gives apparent k median **0.075** (95% CrI **0.050–0.105**) yr^-1, LSR loss median **2.45** (CrI **0.12–7.14**) mg/L; max basic split-Rhat **1.011**. A separate long four-chain fit with an oxygen-poor **Beta(2,5)** prior on LSR entry DO/9.5 gives k median **0.067** (CrI **0.045–0.098**) and LSR loss median **5.40** (CrI **2.45–8.30**); max split-Rhat **1.021**. In the first model the LSR loss–k posterior correlation is **−0.391** and river DO–k correlation is **+0.363**. These numbers were obtained by the Python quadrature/MH prototype, **NOT by running the R/Stan implementation**.

A separate importance-weighting screen of the original posterior showed that doubling empirical prior SDs for chloride, δ18O and excess-air endmembers hardly shifted median k or LSR oxygen loss; a strongly oxygen-poor LSR DO prior had a much greater effect. Importance-weighting ESS is not MCMC effective sample size; direct refitting provides a stronger check.

**Mixing fraction insight:** On the 41-well prototype, switching from tracer-only to tracer+DO posterior mixing fractions shifted median fLSR by less than 0.001 at essentially every well (median absolute change about 0.00005). Thus the DO observations contribute almost no information about the mixing fractions in this model. The larger differences from the old 29-well three-tracer median mixing estimates reflect the changed tracer likelihood, estimated endmembers, and inclusion of excluded wells—not oxygen-data feedback. This is model-conditional.

**Unresolved:** Posterior predictive error and cross-validation of the final joint model have **not** yet been established. The original tracer endmember centres come from the same study, and two-recharge-component mixing may be too simple. Conclusions about LSR-versus-aquifer DO reduction remain conditional on the assumed river DO, MRT prior, mixing structure and common age-associated k. Do not interpret a positive k as independent proof of saturated-aquifer oxygen loss.

## Held-out DO validation: common versus variable LSR DO (script 14)

**Decisive predictive test for the added LSR-DO spatial-variability parameter.**

The joint Stan model now has a **has_DO** observation mask. This allows strictly holding out groundwater DO while retaining the independent chloride, δ18O, excess-air observations and reported MRT for the same wells. Script 12 passes has_DO=1 for all wells; script 14 sets has_DO=0 for the withheld wells AND overwrites their DO values with zeros before passing data to Stan. Excess N2 is excluded (use_n2=0) so no potentially reactive nitrogen information leaks into DO prediction.

Run in Positron after CmdStan/Stan dependencies are installed:

    source("scripts/14_joint_Bayesian_LSR_DO_grouped_CV.R")

The two alternatives use identical priors and the exact same grouped five-fold splits (related Eastbourne, Gilbertson, Flaxmere and exact repeat Site_ID samples are kept together):

- **Common effective LSR DO:** the global effective LSR DO parameter applies to all wells; prior Uniform(0,9.5 mg/L).
- **Variable effective LSR DO:** each well has a logistic-normal source DO between 0 and 9.5 mg/L, controlled by a common central DO and hierarchical SD tau; prior on tau is HalfNormal(SD=1.5) on the logit scale. No site-specific DO values are fitted from held-out groundwater DO.

Both models estimate source fractions and tracer end members jointly; reported MRT is mean-centered under an assumed log-normal error with log SD=0.25. The fit uses all **41 wells with usable DO and MRT**, not only the 29 that pass the old tracer-agreement screen. Posterior predicted mean DO is used for RMSE and MAE, with posterior predictive intervals for uncertainty calibration. Four-chain HMC diagnostics (Rhat, bulk ESS and divergences) are checked after every fit.

Outputs (local only, gitignored):
- outputs/joint_LSR_DO_CV_folds.csv
- outputs/joint_LSR_DO_CV_predictions.csv
- outputs/joint_LSR_DO_CV_convergence.csv
- outputs/joint_LSR_DO_CV_summary.csv
- outputs/joint_LSR_DO_CV_paired_errors.csv
- plots/joint_LSR_DO_grouped_CV_predictions.png

**Evidence standard:** an improvement must be greater than fold-to-fold and MCMC uncertainty, and not depend entirely on the assumed spatial-variability prior. On 41 wells, a small RMSE difference is not evidence of identified spatial LSR DO variation. Source-specific DO loss remains conditional on an assumed 9.5 mg/L LSR reference, an assumed 7.5 mg/L river DO reference, the source mixing structure, and the MRT prior.

**Implementation status:** script 14 and the Stan observation-mask modification are committed, but R/Stan has not been executed here. An independent Python joint posterior integration/CV check is being used to assess expected predictive behavior; its results must not be represented as R/Stan output. Since the source end-member priors were derived partly from the same study, the validation conditions on those empirical prior centres and is not a fully independent external-tracer validation.

### Independent five-fold numerical check of common versus variable LSR DO

An independent Python approximation to the same 41-well joint Bayesian model was fitted using **all three conservative tracers plus MRT uncertainty**. Reported MRT was mean-centered under a log-normal error model with log SD=0.25. Posterior mixing fraction and MRT uncertainty were numerically marginalized using fixed quadrature; site-specific LSR DO variation was integrated under the same bounded logistic-normal distribution used in Stan. Held-out DO was not included in each training posterior; held-out tracers and MRT remained known. The models used five matched groups of 9/8/8/8/8 wells (keeping related Eastbourne, Flaxmere, Gilbertson and repeat Site_ID wells together).

**Matched grouped-CV results (mg/L):**

| Subset | n | Common LSR DO: RMSE / MAE | Variable LSR DO: RMSE / MAE |
|---|---:|---:|---:|
| All usable wells | 41 | **1.8148 / 1.3084** | **1.8123 / 1.3094** |
| LSR-dominated (tracer-derived fraction >= 0.5) | 5 | 1.9251 / 1.5223 | 1.9245 / 1.5272 |
| Other wells | 36 | 1.7989 / 1.2787 | 1.7962 / 1.2791 |

The full-data overall RMSE difference is only **−0.0024 mg/L** (variable minus common). Variable LSR DO has lower absolute error in 18 wells and higher error in 23 wells. A descriptive *five-fold block bootstrap* for the RMSE difference spans approximately **−0.0095 to +0.0121 mg/L**; with only five groups this is not a reliable formal uncertainty interval. The improvement is too small to support adding LSR DO heterogeneity for predictive purposes. The apparent gain is likewise negligible among the **five LSR-dominated wells**.

**MCMC qualification:** Three independent Python chains per fold/model were run with 850 warmup and 1000 retained draws. Two problematic fits were rerun with four chains, 1650 warmup and 1750 retained draws. Maximum basic split-Rhat after these reruns was 1.056 (variable-LSR model, one fold); the remaining 9 fits were at most 1.035. Split-Rhat alone is not a full effective-sample-size or convergence assessment. All numbers above come from **Python, not a compiled/executed R/Stan model**. The R/Stan code in script 14 is provided for a separately reproducible validation and may yield slightly different results.

**Interpretation:** This does not prove spatial variation in LSR oxygen loss is absent. Rather, with only five strongly LSR-dominated wells, the larger model provides no observed improvement in out-of-sample DO prediction under the current tracer mixing structure, river and LSR DO priors, and MRT uncertainty assumptions. Do not use its extra parameters to claim an identified LSR-versus-aquifer oxygen-loss split. Endmember prior centres remain partly empirical from the same study.


## Tracer-only Bayesian mixing-model validation (script 15; 31 Oct 2026)

**Scientific rationale:** Before interpreting LSR versus groundwater-age DO reduction, establish how well a two-source model can predict the three measured, nominally conservative tracers. DO, groundwater MRT, nitrate and **excess N2** are deliberately excluded from this comparison. The third mixing tracer is **excess air (xs_Air)**, not excess N2. Excess air itself may reflect gas exchange and recharge physics, and its use as a conservative mixing tracer is a hypothesis being tested.

**Three candidate models fit the same source-fraction and endmember priors:**
- **A: common end members:** Two sources with one global river and one global LSR value for each tracer, plus study-level tracer residual SDs. This provides a normal-likelihood reference for the earlier two-endmember model; it is not mathematically identical to the earlier robust Student-t mixing likelihood.
- **B: variable LSR isotope signature:** Add a half-normal-prior between-well SD for the LSR **δ18O** source value. The LSR isotope random effect is analytically integrated into the δ18O likelihood, contributing `(fLSR*tau_O18)^2` to variance; no site-specific isotope measurement is fabricated.
- **C: variable LSR δ18O plus nonnegative chloride enrichment:** Adds an independent half-normal chloride concentration term with unknown population scale. It is analytically integrated using a skew-normal likelihood. This is an empirical mixing-model discrepancy term, **not proof of anthropogenic contamination**.

Source-endmember priors are centred on study-derived endmember samples and have finite SDs. These empirical prior centres are not independent field measurements. The script exposes `cfg$endmember_prior_scale` to widen all such priors; compare scale 1 versus 2 before interpretation.

**Stan and R files:**
- `stan/15_tracer_only_mixing_models.stan`
- `scripts/15_tracer_only_mixing_CV.R`

Run in Positron with CmdStanR, CmdStan, tidyverse, posterior and here already installed:

```r
source("scripts/15_tracer_only_mixing_CV.R")
```

**Predictive design:** Every sample retains at least two observed tracers. Across THREE reproducible, balanced and group-aware validation rotations, EACH observed chloride, δ18O and excess-air concentration is withheld exactly once. Each posterior fit has only that site's other measured tracers available to infer the source fraction. The held-out value is masked in the likelihood **and overwritten with zero in the Stan data**. All three mixing models use the SAME masks; this evaluates actual tracer prediction instead of in-sample flexibility. Group-aware masks keep Eastbourne/Gilbertson/Flaxmere repeats together. No groundwater DO is supplied to the model.

The script performs 9 four-chain Stan validation fits, optionally followed by 3 all-data fits for per-well fLSR posterior intervals, saving local gitignored files:
- `outputs/tracer_only_CV_masks.csv`
- `outputs/tracer_only_CV_predictions.csv`
- `outputs/tracer_only_CV_scores_by_tracer.csv`
- `outputs/tracer_only_CV_scores_overall.csv`
- `outputs/tracer_only_CV_convergence.csv`
- `outputs/tracer_only_model_source_fractions.csv`
- `outputs/tracer_only_model_endmembers.csv`
- `plots/tracer_only_heldout_predictions.png`
- `plots/tracer_only_normalised_CV_error.png`

**Independent preliminary screening, NOT a converged Bayesian result:** A separate SciPy *maximum-a-posteriori* calculation with source fractions numerically marginalized under a uniform prior was run on 42 private records, with 124 withheld tracer measurements. It used the same three structural alternatives and matched 3-rotation observation masks; its normal error model is deliberately simple. The held-out point-prediction RMSEs were:

| Tracer | A common | B variable δ18O | C plus Cl enrichment |
|---|---:|---:|---:|
| Chloride (mg/L, 41 observations) | 12.27 | 10.66 | 10.66 |
| δ18O (‰, 41 observations) | 0.389 | 0.385 | 0.385 |
| Excess air (original units, 42 observations) | **1.177** | 1.209 | 1.209 |

The aggregate held-out log predictive score **conditional on MAP parameters** (higher is better) was A **−279.52**, B **−281.58**, C **−281.58**; differences around 2 log-score units over 124 observations are not decisive. For this particular screen, fitted positive chloride enrichment approached zero, and B/C gave almost indistinguishable predictions. Variable LSR δ18O modestly improved Cl while worsening excess-air prediction. Endmember estimates shifted appreciably between masking rotations, highlighting weak identification. **Do not report these MAP values as Stan posterior or as validated Bayesian CV.**

A separate independent custom Metropolis posterior attempt showed unsatisfactory chain mixing (early split-Rhat maxima ~1.13–1.40), so it has **NOT** been used to draw conclusions. R/Stan code has been structurally checked but **not compiled or executed** in this environment. Credible intervals, latent mixing fractions and posterior predictive comparisons must come from the compiled script 15 with acceptable Rhat, ESS, divergences and 95% predictive coverage.

**Interpretation:** Present evidence does not justify adding an unconstrained chloride enrichment term to the current mixing model. However the strong δ18O/Cl/excess-air disagreement in several wells remains real. If model B/C fails to improve held-out prediction, investigate whether excess air behaves as a source tracer and/or whether a two-recharge-component model is inadequate rather than adding further uncertain source fractions. Do NOT update the oxygen reduction model with new fractions until this comparison has passed convergence and predictive checks.

## Original CSV: other candidate tracers and pairwise mixing (script 16)

**Source:** private original/prepared groundwater CSV (currently 43 records, 42 with reported MRT). Script 16 uses the private local `data/groundwater.csv` and needs no Stan installation. It compares the three possible two-tracer mixing estimates and the prediction of the omitted tracer, audits additional analytes, and assesses effects on the same set of screened wells in a simple first-order DO fit.

Run:

    source("scripts/16_pairwise_mixing_and_extra_tracer_screen.R")

The fixed end members follow the original definition: river sample 2801; LSR endmember is the mean of sites 1032 and 4294 and the 108-m sample from HB Prison. Each individual tracer gives a raw fLSR, clipped to [0,1] when forming the pairwise mean. The omitted tracer is predicted without using that well's observed value in the fraction estimate. Scores are reported with and without the four endmember samples; **because endmembers are derived from the study dataset, this is not an independent or fully Bayesian CV**. For refitted Bayesian endmember models with masked tracer observations see script 15.

**Independent fixed-endmember check of current 43-row prepared data:**

| Tracers used | Predicted tracer | n | Held-out tracer RMSE | Excluding four endpoint samples |
|---|---|---:|---:|---:|
| δ18O + Cl | excess air | 41 | 0.894 (xs_Air units) | 0.909 |
| δ18O + excess air | chloride | 41 | 11.82 mg/L | 12.29 mg/L |
| Cl + excess air | δ18O | 41 | 0.415‰ | 0.432‰ |

These RMSEs have different units and should not be compared as though commensurable. An alternative check used the **same 30 wells** that currently pass the three-tracer agreement QC (threshold 0.35). The fixed-endmember, conditional first-order k estimates are 0.0701 yr^-1 (three-tracer median), 0.0714 (Cl+air), 0.0674 (δ18O+air), and 0.0673 (δ18O+Cl). Corresponding effective LSR DO estimates are 8.46, 9.03, 7.46 and 7.42 mg/L. This reiterates that the apparent age association is comparatively robust while the LSR DO deficit is not.

**Other analytes:** the raw/prepared file contains 20 negative-valued H2 results, very likely stable hydrogen isotope δ²H rather than dissolved H2, but metadata/units must be verified and the river endpoint is absent. δ²H and δ18O have Pearson r=0.945 in the 20 paired records, so both may help diagnose recharge provenance and evaporation but are not automatically independent mixing constraints. Sodium (n=42) correlates very strongly with chloride (r=0.975), and both electrical conductivity columns are largely redundant with other ions. Mg, Ca, sulfate, K and silica may help evaluate water-rock reaction or source-dependent solute composition, but cannot simply be treated as conservative. Nitrate, Fe, Mn, NH3, methane and **excess N2** should remain distinct reactive/redox variables at this stage. The existing categorical `source` field labels 28 samples as rain, 14 as river and one as missing; the classification is **not independent ground truth** until its provenance is established.

**Sample-count warning:** older scripts assume 42 input records and 29 selected wells. The currently inspected prepared file contains 43 records, 41 with all three mixing tracers, and 30 passing the previous QC. Its missing-MRT sample is Frimley (also missing chloride). Keep track of data version before directly comparing fitted results; do not silently interpret a changed sample selection as a model improvement.

Local gitignored outputs:
- `outputs/mixing_pairwise_endmembers.csv`
- `outputs/mixing_pairwise_all_fractions.csv`
- `outputs/mixing_pairwise_predictions.csv`
- `outputs/mixing_pairwise_scores.csv`
- `outputs/mixing_pairwise_effect_on_DO_model.csv`
- `outputs/mixing_additional_tracer_screen.csv`
- `outputs/mixing_source_label_comparison.csv`
- `plots/mixing_pairwise_predict_third_tracer.png`
- `plots/mixing_pairwise_source_fraction_disagreement.png`

**Next priority:** establish whether xs_Air behaves as a conservative mixing tracer and whether variable δ18O endmembers improve predictions after full Bayesian model-15 convergence; don't add sodium, δ²H, conductivity, sulfate or redox indicators to the mixing likelihood without independent physical justification and suitable endmember data.

## Four-tracer source mixing with sodium (script 18; 31 October 2026)

**Research question:** Does adding sodium to the existing δ18O, chloride and excess-air tracers improve source-mixing estimates and predict other withheld tracers, without confusing the tightly correlated major ions with independent endmember information?

This is still **TWO recharge sources** (river recharge and land-surface recharge, LSR), **not four sources**. It does **not** use DO, groundwater MRT, nitrate, excess N2 or the preassigned source categories.

New files:
- `stan/18_four_tracer_joint_mixing.stan`
- `scripts/18_four_tracer_mixing_CV.R`

Run with CmdStanR and CmdStan installed, using the locally prepared **private** `data/groundwater.csv`:

```r
source("scripts/18_four_tracer_mixing_CV.R")
```

The Stan model provides a **matched three-tracer control** (δ18O, Cl, excess air) and a **four-tracer** model (δ18O, Cl, Na, excess air) with identical priors on shared endmembers and tracer noise. Both include a bounded latent LSR fraction per well. **Temporal/spatial LSR δ18O variability** is represented as a zero-mean random endmember deviation integrated out of the likelihood: its additional δ18O variance is `fLSR^2 * tau_LSR_O18^2`, with `tau_LSR_O18 ~ HalfNormal(0.8‰)`. This does not estimate rainfall isotope seasonality or hydrograph-based recharge history.

The joint Na–Cl likelihood uses a bivariate normal conditional formulation when both ions are present, accounting for a fitted **residual correlation after conditioning on source mixing**, with an informative zero-centred prior on rho. When one ion is withheld, the remaining observed ion enters the other ion's conditional posterior predictive distribution. This is essential for fair validation of correlated ions. An assumed analytical error model assigns minimum analytical SDs of **1 mg/L Cl**, **0.5 mg/L Na**, and 3% of concentration if larger. These analytical uncertainties are **assumptions**, not supplied by the source CSV. The model fits additional study-level discrepancy SDs and uses the available xs_Air uncertainty. Na can be altered by cation exchange; no explicit exchange stoichiometry is inferred here.

The priors for river Na and LSR Na (centres 6.6 and 44.33 mg/L respectively), together with chloride and δ18O prior centres, are derived from the **same groundwater study**. These priors are not an external endmember survey. Repeat the experiment using `cfg$endmember_prior_scale = 2` in the R script to examine prior sensitivity, and consider external recharge-water time series before interpreting the mixture as ground truth.

**Validation:** Four deterministic rotations mask exactly one measured tracer per well in each fit; every available tracer measurement is withheld exactly once over the four rotations. Related Eastbourne, Gilbertson and Flaxmere well samples share the same held-out tracer within each rotation. A withheld tracer is removed from the likelihood **and its input value is overwritten with zero** before the fit. The held-out concentration is predicted from the other measurements; four-tracer Na and Cl forecasts use their appropriate conditional distributions. The fair three-versus-four model ranking uses ONLY shared targets: δ18O, chloride, excess air. **Predictive Na results are a separate assessment of the four-tracer model**, rather than scored as an advantage over a baseline that never uses Na.

Compare **per-tracer RMSE, MAE, bias, normalised RMSE, 95% predictive coverage, paired absolute errors, MCMC Rhat/ESS and divergences**. The optional full-data fits compare posterior LSR fraction medians and intervals. Do not conclude sodium improves mixing simply because the four-tracer model fits Na well *in sample*.

Outputs (local, all gitignored):
- `outputs/four_tracer_CV_masks.csv`
- `outputs/four_tracer_CV_predictions.csv`
- `outputs/four_tracer_CV_convergence.csv`
- `outputs/four_tracer_CV_scores_by_tracer.csv`
- `outputs/four_tracer_CV_scores_shared_outcomes.csv`
- `outputs/four_tracer_CV_paired_predictions.csv`
- `outputs/four_tracer_full_data_LSR_fractions.csv`
- `outputs/four_tracer_full_data_endmembers.csv`
- `outputs/four_tracer_LSR_fraction_change.csv`
- `plots/four_tracer_CV_predictions.png`
- `plots/four_tracer_CV_normalised_error.png`

**Known original-data diagnostic:** The currently inspected prepared CSV contains 43 records; 42 contain Na and Cl. Na–Cl Pearson correlation is **0.975**, but this primarily describes the common concentration gradient, not necessarily conditional residual correlation. The endpoint mixing line is Na(mg/L) = **3.244 + 0.6848×Cl(mg/L)**. Seven of 42 sodium concentrations differ from that fixed-endpoint line by over 3 mg/L, supporting a test of mixing and chemical departures. This is descriptive evidence only, not a completed Bayesian predictive result. The data include isotope variation and source-endmember uncertainties that may remain underidentified.

**Implementation status:** The Stan/R scripts were committed and statically checked, but **R, CmdStan and CmdStanR are not available in this execution environment**. No posterior predictive scores, credible intervals or chain-convergence results have been fabricated. Please run the script in Positron; don't use or propagate revised source fractions to the DO model until there is evidence of held-out improvement and adequate convergence.


### Independent four-tracer numerical results (31 October 2026)

An **independent Python Bayesian importance-sampling implementation** of the same structural two-source mixing model was run on **43 records** from the original private `Uwe_HP_data_prepared.csv`. Source fractions fLSR were integrated out numerically under a Uniform(0,1) prior using Gaussian quadrature; endmember and error parameters were integrated by importance sampling with a multivariate-t proposal. The original dataset was never added to this public repository.

**Important implementation differences from Stan:** this independent implementation uses 21 quadrature nodes, approximates excess-air analytical error as a constant 0.15 in original units rather than per-sample values, and approximates posterior means by resampling weighted global-parameter draws. It does not produce Stan NUTS R-hat, ESS, divergences or full predictive coverage. Keep results separate from the not-yet-executed R/Stan model.

Four group-aware matched masking rotations withheld each of the available δ18O, Cl, Na, and xs_Air observations exactly once. During a held-out δ18O/Cl/air prediction both models were scored against the SAME outcome and observations; Na could contribute to the four-tracer model's conditioning measurements. Na forecast scores were assessed separately. All endmember prior centres were empirical from the same groundwater study. The common three-tracer comparison used the same allowance for variable LSR isotope signature and fitted extra tracer scatter as the four-tracer model.

**Posterior-averaged held-out RMSE, default endpoint prior widths:**

| Held-out tracer | n | Three tracers | Four tracers | Relative change |
|---|---:|---:|---:|---:|
| δ18O (‰) | 42 | 0.388 | **0.366** | −5.8% |
| Cl (mg/L) | 42 | 10.713 | **3.684** | −65.6% |
| excess air (source units) | 43 | 1.101 | **0.918** | −16.6% |
| Na (mg/L; four-tracer model only) | 42 | not modelled | **2.541** | n/a |

The **enormous chloride gain is expected**: sodium measured at the same well predicts chloride well, given their strong correlation. The meaningful nonion evidence is improved δ18O and excess-air prediction, and its magnitude is smaller. With endpoint prior SDs doubled, the corresponding three- and four-tracer RMSEs were δ18O 0.389→0.366‰, Cl 10.889→3.674 mg/L, and excess air 1.136→0.933. The conclusion thus persists under that one prior-width perturbation.

The importance-sampling effective sample size (ESS) across the eight default CV fits was **minimum 209**, with 12,000 candidate draws per proposal round. Doubled-prior-width fits used 9,000 candidate draws per round and had minimum ESS **96**, so the wider-prior results require stronger confirmation. Importance ESS does not establish convergence of independent Stan chains. A matched **group-level nonparametric bootstrap** for the RMSE difference (four minus three) has approximate descriptive 2.5–97.5% ranges: chloride −9.70 to −4.71 mg/L, δ18O −0.046 to +0.005‰, excess air −0.391 to −0.003 (source units). These bootstrap intervals ignore global endmember prior dependence and numerical posterior approximations and are not formal model-selection probabilities.

Excluding the four wells used to define the original endmembers (2801, 1032, 4294, and HB Prison depth 108m) still yields four-tracer advantages: Cl 9.59→3.49 mg/L, δ18O 0.400→0.378‰ and excess air 1.048→0.888. This is NOT fully external validation because the endpoint prior centres still came from these same wells.

**Full-data posterior source fractions** (43 wells): mean fraction changed by median absolute **0.0138**, mean absolute **0.0239**, maximum absolute **0.0963** between the three- and four-tracer models. Median width of an approximate 95% LSR-fraction posterior interval fell from **0.172 to 0.108**. Six samples were LSR-dominated in both versions. The four-tracer joint Na–Cl residual correlation posterior was broad: median about **−0.09**, approximate 95% range **−0.62 to +0.39**. The LSR δ18O variability parameter was also weakly identified, with four-tracer median about **0.19‰** and approximate 95% interval **0.007–0.80‰**. The modest shift and narrower posterior interval should not be confused with independent proof of conservative sodium transport.

**Conditional oxygen-model sensitivity (not a new joint Bayesian oxygen model):** using the same 30 three-tracer-QC-screened wells and fitting a common effective LSR DO and first-order k to fixed posterior-mean source fractions gave k **0.0716 yr−1** (3-tracer fractions) versus **0.0714 yr−1** (4-tracer fractions). Leave-one-well-out DO RMSE **1.007 versus 1.030 mg/L**, respectively. Thus **the improved mixing-tracer predictions have not translated into improved groundwater DO prediction** in this screening calculation. The fitted effective LSR DO was close to its upper bound (9.5 mg/L) under both posterior-mean fractions, demonstrating that LSR-source loss remains weakly identifiable. This screening holds fractions fixed and must NOT be interpreted as a joint Bayesian oxygen-loss partition.

**Crucial local-input warning:** The older oxygen-project \`data/groundwater.csv\` is a **42-row reduced dataset without Na**. Script 18 now prefers local \`data/Uwe_HP_data_prepared.csv\` and gives an actionable error if sodium is missing. This 43-record dataset is stored in the user's **PRIVATE** original \`Lincoln-Environmental/Oxygen-Reduction\` repository under \`data/Uwe_HP_data_prepared.csv\`. Copy it into the public modelling repo's local \`data/\` folder without committing it. Keep all groundwater data files private and gitignored.

**Research interpretation:** Sodium adds predictive hydrochemical information to this formulation and sharpens mixing fractions, but Na can be modified by cation exchange; the uncertain conditional Na–Cl residual correlation and study-derived endmembers limit causal interpretation. The true physical LSR fractions cannot be independently validated from internal held-out tracer scores. Confirm with a compiled Stan run, a wider-endmember sensitivity, posterior predictive coverage and ideally external endmember time series (especially rainfall isotope composition) before transferring these fractions into the formal Bayesian oxygen budget.


## Fully joint Bayesian DO refit: three versus four tracers (script 19)

**Question:** Does the sodium-related reduction in posterior source-mixing uncertainty constrain LSR-associated oxygen depletion separately from the apparent age-dependent DO coefficient, *and* does it improve held-out groundwater DO prediction?

### Design

- `stan/19_joint_three_vs_four_tracer_DO.stan` — identical two-recharge-source joint models distinguished by `use_na=0` (δ18O, chloride, excess air) versus `use_na=1` (add Na). Both include uncertain tracer endmembers, per-well fLSR, fLSR-dependent isotope endmember variability, mean-centred uncertain true MRT, river entry DO, common effective LSR entry DO, and k. Na and Cl use a **joint correlated likelihood**, not an independently counted Na/Cl ratio.
- `scripts/19_joint_three_vs_four_tracer_DO.R` — fits two 4-chain full models, exports global and per-well posteriors with MCMC diagnostics, and runs **five-fold matched grouped held-out DO validation** with the same wells, DO priors and MRT assumptions. Held-out DO is removed from the Stan likelihood and **zeroed in its data input**; source tracers and reported MRT remain available for held-out wells. Posterior predictive DO intervals include Gaussian observation discrepancy.
- The original **private** `data/Uwe_HP_data_prepared.csv` is required, not the reduced oxygen-workflow CSV lacking sodium. The private input is gitignored. Current prepared data: **43 samples, 42 with DO+MRT**, compared with 41 in an earlier reduced dataset. The river-endmember site 2801 has DO 10.2 mg/L and reported MRT 0; this challenges the river entry DO prior centred at 7.5 mg/L. Script 19 has an optional flag to exclude only that DO observation in a sensitivity fit; baseline retains it.

Run in Positron from the repository root after installing CmdStanR, CmdStan and posterior:

```r
source("scripts/19_joint_three_vs_four_tracer_DO.R")
```

The default performs **12 Stan fits** (two full and ten grouped CV), which is computationally expensive. Set `cfg$run_CV=FALSE` for the two full fits first or `cfg$run_only_fold=1L` for a CV smoke test. Setting `cfg$run_prior_sensitivity=TRUE` also fits the four-tracer model with doubled endmember prior SDs. The default common LSR entry DO prior is Uniform(0,9.5) mg/L; the pre-infiltration LSR oxygen reference 9.5 mg/L is **assumed**, not measured.

Local (gitignored) outputs include:
- `outputs/joint_DO_19_posterior_parameters.csv` — median and central 95% credible intervals.
- `outputs/joint_DO_19_per_well_posterior.csv` — site-specific source fractions, DO and partition component summaries.
- `outputs/joint_DO_19_source_fraction_change.csv`.
- `outputs/joint_DO_19_CV_scores.csv`, `outputs/joint_DO_19_CV_predictions.csv`, `outputs/joint_DO_19_CV_paired_errors.csv`.
- `outputs/joint_DO_19_full_diagnostics.csv`, `outputs/joint_DO_19_CV_diagnostics.csv`.
- `outputs/joint_DO_19_full_Na0.rds`, `outputs/joint_DO_19_full_Na1.rds` (posterior draws retained locally).
- `plots/joint_DO_19_posterior_comparison.png`, `plots/joint_DO_19_grouped_CV_predictions.png`.

### Independent numerical execution (31 October 2026): provisional

The same full joint source-tracer/DO/age model was independently implemented in Python with **21-node Gauss–Legendre integration over per-well fLSR**, **five-node age integration**, and **four-chain adaptive multivariate Metropolis sampling of shared model parameters**. Posterior latent f and age were also explicitly sampled conditionally when reporting per-well mixing fractions and modelled DO losses. This is a genuine Bayesian posterior approximation, **not a Stan/HMC run**; it uses a constant excess-air analytical SD of 0.15 instead of per-well values. The independent run uses mean-centred lognormal MRT log SD 0.25, the same empirical endmember prior centres, and the common beta(1,1) LSR entry DO prior. All private per-well draws/observations are kept off public GitHub.

**Full-data posterior comparison (approximate 95% credible intervals):**

| Quantity | 3 tracer | 4 tracer |
|---|---:|---:|
| Apparent k (yr^-1) | **0.074** (0.051–0.104) | **0.074** (0.050–0.105) |
| Effective LSR DO deficit (mg/L) | **2.47** (0.13–7.02) | **2.17** (0.15–6.15) |
| Arithmetic mean modelled LSR-pathway DO loss per well (mg/L) | 0.70 (0.04–2.01) | 0.59 (0.04–1.68) |
| Arithmetic mean apparent age-associated DO loss per well (mg/L) | 4.54 (3.23–5.81) | 4.53 (3.17–5.79) |
| Approximate LSR share of modelled losses averaged across sampled wells | 0.131 (0.008–0.364) | 0.117 (0.007–0.329) |

These modelled loss fractions are **not volume-weighted catchment or aquifer oxygen budgets**; they are averages over a nonrandom set of wells. An apparent age-associated loss is not proof of oxygen consumption exclusively in the saturated aquifer. The LSR DO deficit is strongly prior-sensitive and depends on the assumed 9.5 mg/L reference.

Adding Na reduced the median site-level fLSR 95% posterior interval width from ~**0.172** to ~**0.108**, with a median absolute change in posterior mean fraction of ~**0.014**. This improvement in precision depends on the assumed source mixing model and Na transport behavior.

**Five-fold matched held-out DO predictions (42 wells; fold sizes 9/9/8/8/8):**

| Model | RMSE DO (mg/L) | MAE DO (mg/L) |
|---|---:|---:|
| 3 tracer | **1.8458** | **1.3196** |
| 4 tracer | 1.8461 | 1.3214 |

The 4–3 difference in RMSE is **+0.00035 mg/L**, essentially zero. Among 42 held-out wells, the 4-tracer model has lower absolute error in 20 and higher error in 22. This differs from older DO-CV figures for differently screened data and models; compare ONLY the matched conditions above.

**Convergence caveat:** the independent Python full-model maximum split-Rhat was approximately 1.058 for the 3-tracer global parameters and 1.285 for the 4-tracer parameters. The largest issues concerned isotope heterogeneity, ion-correlation and other tracer nuisance variables; full-model k split-Rhat was 1.005 (three) / 1.018 (four), and LSR DO 1.014 / 1.005. Some grouped CV tracer nuisance parameters had split-Rhat **as large as 1.53**, although the key k/LSR DO parameters showed substantially better agreement. This means the full inference **has not passed convergence checks**, and the numeric confidence intervals and tiny CV differences are **provisional**. The independent numerical approximation does NOT validate the still-uncompiled R/Stan fit, which may give different numerical results.

**Scientific conclusion, conditional on assumptions:** Sodium sharpens estimated source fractions, but the inferred LSR DO deficit is still broad and dependent on its prior; k and held-out DO predictive error are essentially unchanged. Better source mixing alone is insufficient to identify the physical oxygen-consumption partition. Validate the Stan fits and then focus on source-entry DO and variation in reaction/groundwater-age processes, rather than adding another weakly constrained tracer term.


## Paired Bayesian prediction of DO and excess N₂ (script 20)

**31 October 2026 — new joint experiment.** This tests whether the reactive gas signature adds information to the apparent oxygen-depletion history and improves predictive performance, rather than simply adding another conservative mixing tracer. The mixing model still has **two recharge sources and four mixing tracers** (δ18O, Cl, Na, excess air). Excess N₂ is a **separate reactive outcome**; its response is NOT a stoichiometric measurement of oxygen depletion.

New files:
- `stan/20_joint_paired_DO_excessN2.stan` — DO + N₂ joint likelihood with four-tracer mixing, Na–Cl residual correlation, source isotope variability and mean-centred uncertain true MRT.
- `scripts/20_joint_paired_DO_excessN2.R` — full-data posterior comparison and paired five-fold holdout of **BOTH DO and excess N₂** at exactly the same wells, retaining measured conservative tracers and reported MRT. Produces pair-specific predictive scores, uncertainty coverage, well-specific modelled depletion and convergence diagnostics.

Three competing models for **excess N₂**, plus a DO-only reference:
1. **DO_only**: unchanged four-tracer joint oxygen model; excess-N₂ measurements not used.
2. **N2_age**: excess N₂ as a study-conditional source/background plus β_age log(1 + true MRT). This is an empirical **uncoupled** reference.
3. **N2_DOloss**: baseline + positive β_link times the **modelled apparent cumulative DO depletion**. This extra coupling is a statistical association, **not an O₂/N₂ reaction stoichiometry**.
4. **N2_lowDO**: baseline + positive β_link times the modelled residence duration below a **fixed** DO threshold, default 2 mg/L (expressed in multiples of 20 years). Evaluate thresholds 1, 2, and 3 mg/L. This assumes a highly simplified piston-flow DO trajectory and does **not** capture microscale oxic/anoxic conditions.

Excess N₂ has a **Student-t(df=4)** observation likelihood with the provided `xs_N2_sig` and separately estimated extra discrepancy SD, so a single unusual N₂ measurement is not forced to determine the rate. This is a robust statistical approximation; negative tails of replicated N₂ should not be interpreted as physically meaningful negative denitrification. The n2_R/n2_LSR terms model source-associated background differences *empirically*; they do not assert excess N₂ is conservative. No measured nitrate, derived `initial_NO3N` or derived `reaction_progress` is used as a predictor. The latter two already contain the excess-N₂ observation, and would be **data leakage** during held-out N₂ prediction. Nitrate availability at recharge, redox microenvironments, and other electron-acceptor reactions remain unobserved limitations.

**Run in Positron** with the original private `data/Uwe_HP_data_prepared.csv` and CmdStanR/CmdStan installed:

```r
source("scripts/20_joint_paired_DO_excessN2.R")
```

Default fits four four-chain full-data models and 20 four-chain grouped-CV fits. To reduce computation, set `cfg$run_CV=FALSE` for the full posterior comparison only, or `cfg$run_only_fold=1L` for a single validation fold. The script's `cfg$exclude_16758=TRUE` option reruns everything without site 16758, which reports **excess N₂ 8.8, DO 0.3 mg/L, MRT 70 years**. The full original dataset has 43 records, of which 42 have usable DO, excess N₂ and MRT (Frimley lacks MRT). Excluding 16758 leaves 41 usable samples. `cfg$low_DO_threshold` allows 1, 2 and 3 mg/L threshold tests. Results have descriptive filenames incorporating the n, outlier exclusion and threshold.

Outputs are **local, gitignored**:
- `outputs/joint_DO_N2_20_full_posteriors_*`, `joint_DO_N2_20_per_well_*`, `joint_DO_N2_20_convergence_*`.
- `outputs/joint_DO_N2_20_CV_predictions_*`, `joint_DO_N2_20_CV_summary_*`, `joint_DO_N2_20_CV_diagnostics_*`.
- `plots/joint_DO_N2_20_CV_*`, `plots/joint_DO_N2_20_posterior_*`.

### Independent preliminary paired-CV screening (not Bayesian)

A separate Python SciPy **regularized nonlinear least-squares** screen used private 42-well DO, MRT, excess N₂ and fixed source fractions inferred from the independent four-tracer **tracer-only** model. It compared the age-only N₂ reference with the cumulative DO-loss and post-2-mg/L-exposure additions. **Both outcomes were held out together** in five matched grouped folds. It does not propagate fLSR, age or endmember posterior uncertainty, and it uses an approximately Gaussian N₂ discrepancy rather than the new Student-t likelihood. Consequently it must not be conflated with the Stan posterior or the previous full joint Bayesian results.

| Model | DO RMSE (mg/L) | excess N₂ RMSE (dataset units) |
|---|---:|---:|
| Uncoupled source + age | **1.713** | 1.369 |
| + apparent cumulative DO loss | 1.716 | 1.380 |
| + duration below 2 mg/L DO | 1.718 | **1.350** |

The slight ~0.019 N₂ RMSE improvement in the low-DO model is not convincing, especially because it worsens DO prediction very slightly. **Excluding site 16758 (n=41)** produced excess-N₂ RMSE **0.734, 0.737 and 0.736** for the three models respectively, and DO RMSE **1.733, 1.733, 1.732**, so no systematic advantage remains. These 41-well results are a second controlled validation, not evidence that 16758 is erroneous.

This outcome is scientifically plausible. Excess N₂ reports accumulated nitrate reduction, requiring nitrate availability, and can persist after oxygen concentrations change or water bodies mix. Present well DO does not measure how much reduction occurred over the groundwater history. An empirical N₂ response term is not a physical oxygen-budget constraint. **N₂ may inform source mixing/age or redox propensity without necessarily constraining effective LSR entry DO, k, or their partition uniquely**.

**Validation status:** The two Stan/R files are committed and have been structurally reviewed, but not executed/compiled in this environment. No Stan posterior credible interval or 5-fold Bayesian prediction score is yet available. Review chain Rhat, bulk ESS, divergences, threshold sensitivity, and results with/without well 16758 before treating joint predictions as evidence of identifiable redox pathways.
