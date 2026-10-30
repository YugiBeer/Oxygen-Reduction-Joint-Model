# Oxygen Reduction — Joint Model

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
- River DO = 7.5 mg/L; **assumed** pre-infiltration rainfall DO reference = 9.5 mg/L.
- For the first-order model:
  
  DO_pred = [fR * 7.5 + fP * (9.5 - L)] * exp(-k * MRT)

  L and k are fitted together. For zero-order reduction, the post-mixing concentration is reduced by k * MRT, bounded below by zero.

Only **effective** rainfall-source DO (9.5 - L) is constrained by the data. An inferred source-associated deficit cannot, by itself, be assigned uniquely to infiltration or the unsaturated zone. The first-order rate is an apparent residence-time-dependent rate, not independent proof of saturated-zone reaction.

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

**Scientific limitation:** These models assume initial recharge DO variability is statistically independent of groundwater mean residence time, conditional on mixing fraction. A positive k establishes a need for *age-associated* depletion under those assumptions. It does **not** establish the site of oxygen consumption (infiltration pathway versus saturated aquifer). The three-tracer recharge fractions are held fixed; excess-air tracer suitability and source fraction uncertainty are not propagated here.

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

**Structural non-identifiability:** An alternative model in which river- and rainfall-recharge DO at entry vary as `R_i = R exp(-b MRT_i)` and `P_i = P exp(-b MRT_i)` produces *exactly the same predictions* as a first-order transport attenuation model with k=b. Fitted DO–MRT observations alone cannot locate oxygen reduction in the infiltration zone versus the saturated aquifer.

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

**Interpret carefully:** Only site 2801 has MRT=0 in the current 29-well selection; the remaining MRTs are ≥1 year. The year-specific attenuation is an assumed travel-time response; MRT alone does not determine dispersion or attenuation. Piston flow would simply delay the sinusoid. The daily-mean variance split cannot be estimated from the 29 well DO measurements, and should be constrained by time-frequency analysis of the actual river record where possible. The age term still cannot distinguish infiltration loss from saturated-aquifer reaction. The code is supplied for execution in Positron; the equivalent numerical sensitivity was independently evaluated in Python, but this R script has not been run here.


## Riverbed oxygen-loss sensitivity and infiltration-only bounds (script 07)

Run after the three-tracer mixing model:

```r
source("scripts/01_three_tracer_mixing.R")
source("scripts/07_infiltration_loss_sensitivity.R")
```

Script 07 tests the remaining major uncertainty in the apparent oxygen reduction rate: the river's monitoring DO (reference 7.5 mg/L) is **not** necessarily the DO entering the saturated aquifer after riverbed infiltration.

- Fix riverbed oxygen loss at 0–7.5 mg/L in 0.25 mg/L increments; for each fixed loss refit effective rainfall-recharge DO (0–9.5 mg/L) and first-order age coefficient k. Repeat all fits **with and without** young, high-DO site 2801.
- Use grouped five-fold CV (related Eastbourne, Gilbertson and Flaxmere wells together), reporting within-sample and held-out RMSE separately.
- Independently calculate the *minimum* riverbed oxygen loss each well requires if **no loss occurs after infiltration**, allowing its rainfall-derived water to be entirely anoxic. This is an extremely favourable assumption for the infiltration-only explanation: oxygen in rainfall recharge would require yet more loss from the river component.
- Repeat required-loss calculation for river references 5.6, 7.5 and 9.4 mg/L. Flag samples, notably site 2801 at 7.5 mg/L, whose DO exceeds the maximum implied by these fixed source bounds. A river sample can temporarily exceed its long-run mean; this is not evidence of oxygen production in the aquifer.
- Compare physically bounded k=0 models with shared effective source DO or different effective river DO for confined versus unconfined wells, against an age-dependent model with river mean fixed at 7.5 mg/L.

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

**Inference limit:** The minimum-loss calculation identifies what a no-aquifer-reduction explanation *requires*, not where consumption actually happened. The confinement model only tests a coarse proxy for variable infiltration conditions; no coordinates or independent infiltration DO observations are available. If initial DO is permitted to vary systematically with MRT, a sufficiently flexible infiltration-only model can be exactly equivalent to the first-order age model. The independent Python calculations provide numerical results; run the R script in Positron to verify its outputs. Raw groundwater data have not been committed.
