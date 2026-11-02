# Heretaunga groundwater oxygen reduction — two-endmember DO model

## Primary analysis: `01_DO_dual_pathway_kinetics`

We model two recharge endmembers—**river recharge** and **land-surface recharge (LSR)**—using only **δ¹⁸O and chloride** to constrain their mixing. Dissolved oxygen (DO) is reactive: each recharge pathway loses oxygen during aquifer residence, and the resulting concentrations are then mixed. The two competing aquifer models are **zero-order** and **first-order** kinetics, fitted independently with uncertainty.

**Sodium and excess air are excluded**, following diagnostic experiments. The model contains no sodium or excess-air parameters. The former model and ablation experiment are preserved for reproducibility, not used in the primary analysis.

### Active files

| Path | Purpose |
|---|---|
| [`scripts/00_setup.R`](scripts/00_setup.R) | **Setup only**: parameters, private data validation, filtered wells, O18/Cl observations and `stan_data`; no MCMC |
| [`scripts/01_DO_dual_pathway_kinetics.R`](scripts/01_DO_dual_pathway_kinetics.R) | **Primary R runner**: zero and first order, diagnostics, outputs, optional grouped CV; sources setup automatically |
| [`scripts/02_DO_plotting.R`](scripts/02_DO_plotting.R) | **Plot-only**: publication-style figures from existing saved posterior CSVs; no refitting |
| [`scripts/03_high_DO_audit.R`](scripts/03_high_DO_audit.R) | **High-DO audit**: observed saturation, physical bounds, posterior comparison for 2130; no MCMC |
| [`scripts/04_river_DO_prior_sensitivity.R`](scripts/04_river_DO_prior_sensitivity.R) | **Controlled experiment**: 3 river-entry DO priors × 2 kinetic orders; unchanged 2801 observations |
| [`stan/04_river_DO_prior_sensitivity.stan`](stan/04_river_DO_prior_sensitivity.stan) | River DO as fixed or one shared uncertain boundary (no extra parameter in fixed scenario) |
| [`scripts/05_river_DO_grouped_CV.R`](scripts/05_river_DO_grouped_CV.R) | **Grouped held-out DO CV**: fixed versus broad river DO prior, both kinetic orders; 2801 masked by default |
| [`stan/01_DO_dual_pathway_kinetics.stan`](stan/01_DO_dual_pathway_kinetics.stan) | **Primary Stan model**: δ¹⁸O/Cl mixing and dual recharge-pathway DO kinetics |
| [`docs/01_DO_dual_pathway_kinetics.md`](docs/01_DO_dual_pathway_kinetics.md) | Scientific protocol, assumptions, outputs, limitations |
| [`data/README.md`](data/README.md) | Private data requirements |

Earlier models `23` and `23a`, their corresponding Stan code and the former `00_check_setup.R` preflight are preserved unchanged in [`archive/`](archive/README.md), alongside older experiments. **Scripts 00 (setup), 01 (baseline fitting), 02 (plotting), 03 (high-DO audit), and 04 (river-entry prior sensitivity) are active.** Setup is separate from sampling. The earlier `archive/scripts/00_check_setup.R` is a different, obsolete preflight.

### Run from Positron

1. Open the **Oxygen-Reduction-Joint-Model** repository folder as the working directory.
2. Copy the unchanged private source `Uwe_HP_data_prepared.csv` into `data/`; it must not be committed.
3. Ensure R packages `tidyverse`, `posterior`, `cmdstanr`, `viridisLite` and CmdStan are installed.
4. To inspect or prepare the model **without sampling**, run:

```r
source("scripts/00_setup.R")
print(cfg)
str(stan_data, max.level=1)
```

5. To **fit** both aquifer kinetic orders, run:

```r
source("scripts/01_DO_dual_pathway_kinetics.R")
```

Script 01 automatically sources `00_setup.R` again to guarantee fresh inputs. All configuration—including `run_CV`, the DO boundaries, MCMC settings and file paths—is defined **only** in `scripts/00_setup.R`. Running `00_setup.R` alone never compiles or samples Stan.

6. **Create or re-create all figures without rerunning MCMC:**

```r
source("scripts/02_DO_plotting.R")
```

Script 02 sources the setup **only to read model inputs** and reads existing CSV summaries from `outputs/01_DO_dual_pathway_kinetics/`. It does **not** compile Stan, load posterior fit objects or call `$sample()`. It prints figures to Positron's **Plots** pane when run in the **R console** and independently saves all PNG files to `plots/01_DO_dual_pathway_kinetics/`.

The initial run fits **two full models only** (four chains per model). The model-runner saves each fitted posterior immediately before producing reports. Results appear in:

```text
outputs/01_DO_dual_pathway_kinetics/
plots/01_DO_dual_pathway_kinetics/
```

First inspect `fit_status.csv` and `sampler_chains.csv`. Then examine `global_posteriors.csv`, `well_posteriors.csv` and `identifiability.csv`. Run **02** to produce the plots. Do **not** interpret model parameters if full-fit diagnostics fail.

**Grouped held-out DO comparison:** After both full fits converge, set `run_CV = TRUE` in the top-level `cfg` list of `scripts/00_setup.R`. It runs five matched grouped folds for both kinetic orders (10 extra fits), saves CV diagnostics and predictions, and performs a matched group bootstrap only if all folds converge. Source-tracer observations are retained for the held-out wells; their DO observations are masked.

### Fixed scientific assumptions

- Rainfall starts at approximate atmospheric DO saturation at **15 °C** and an assumed atmospheric pressure ratio of 1.
- River-to-aquifer entry DO is **fixed at 7.5 mg/L**.
- River-proximal well **2801** retains its **measured DO of 10.2 mg/L** in setup and all data objects. Baseline 01 excludes its measured DO from the likelihood; experiment 04 includes it as one groundwater observation in all three prior scenarios. The **7.5 mg/L value is a river-entry assumption**, never a replacement of the 2801 measurement.
- Far-south anomalous well **16758** is excluded.
- LSR enters the water table with an inferred DO concentration between 0 and rainfall saturation. The result is an **LSR oxygen deficit**, not an identifiable LSR kinetic rate without independent vadose travel time.
- Positive measured MRTs receive non-centred, mean-preserving lognormal uncertainty; proxy 2801 has MRT=0. The model assumes MRT effectively represents aquifer residence time for both recharge sources.
- Zero-order rates are reported in **mg/L/year**; first-order rates in **year⁻¹**. These are alternative aquifer models, not simultaneous reactions.
- Arithmetic averages over wells are **not** flow-weighted catchment oxygen budgets.

A converged sampler demonstrates computational reliability, **not** structural identifiability or validation of the endmembers. Check grouped predictive performance and sensitivity to endmembers, MRT and boundaries before physical interpretation or publication.

**Code verification:** The new files were statically checked and aligned with the previously converged no-sodium/no-air model. Compilation and MCMC execution must be verified in the local R/CmdStan environment.

### Obsolete local `0_setup.R` (one zero)

The previous local `scripts/0_setup.R` is **not** the new `scripts/00_setup.R`. It was an incorrectly named copy of the MCMC runner (contained `cmdstan_model()` and `$sample()`). It is not in the current GitHub tree and cannot be removed remotely from a local working copy. Preserve it outside the repository if needed, but **do not run it**. The correct setup file has **two zeros**: `scripts/00_setup.R`.

### High-DO audit (no MCMC)

The original prepared data have **exactly two wells with measured DO above the fixed 7.5 mg/L river boundary**: **2801 = 10.2 mg/L**, temperature 10.2 °C, MRT 0; and **2130 = 8.8 mg/L**, temperature 16.5 °C, MRT 8 years. Both are approximately 91% of sample-temperature oxygen saturation using the same freshwater polynomial. Proxy 2801's measured DO is masked by the model and must not be assessed as a fitted DO response. Site 2130 is included in the DO likelihood and underpredicted.

Run `source("scripts/03_high_DO_audit.R")` after the two model fits. Script 03 does **not** compile a model or sample Stan. It reads saved `fit_order0.rds`, `fit_order1.rds` and posterior summaries to check:

- source observations and oxygen saturation relative to observed sample temperatures;
- an optimistic DO bound at site 2130 with **no pre-water-table loss** and 100% of the highest-oxygen source;
- the age needed to explain 2130's measured DO under posterior median aquifer rates;
- independent mixing fractions implied by isotope and chloride **prior means**, which may disagree;
- in-sample posterior fitted and posterior predictive DO at 2130, including the posterior for its true MRT (not held-out CV).

Outputs are saved in `outputs/01_DO_dual_pathway_kinetics/high_DO_*.csv` and, when the fit objects are available, `plots/01_DO_dual_pathway_kinetics/09_high_DO_site2130_audit.png`.

**The 2130 posterior check is in-sample and cannot establish external predictive adequacy.** If concerning, next refit with 2130's DO masked but its tracer observations retained, and measure the held-out predictive probability of 8.8 mg/L. Do not discard 2130 solely because it fits poorly. Do not change the 7.5 mg/L river entry boundary based on site 2801 alone.

### Plotting and source colours

The plotting-only script creates seven standard PNG figures: observed vs fitted DO; δ¹⁸O–chloride two-endmember mixing; per-well LSR fractions; LSR DO deficit vs aquifer rate; LSR pre-water-table vs aquifer DO loss; source-pathway DO decay; and DO residuals versus MRT. An eighth plot of held-out DO predictions is generated if `cv_predictions.csv` exists.

All data colours use **Turbo** from `viridisLite`; **river recharge is Turbo blue and rainfall/LSR is Turbo red**. The continuous LSR-fraction scale is explicitly limited to the Turbo blue–red range. Aquifer loss across both sources uses Turbo green so it cannot be confused with river recharge. Neutral grey is used only for reference lines. Graphs are printed explicitly and saved at 300 dpi. The well-wise intervals stored by script 01 are `q5`–`q95` (90% posterior credible), while `global_posteriors.csv` contains central 95% intervals. The fitted DO chart is in-sample; the optional CV chart alone represents held-out predictions.

### River-entry DO prior sensitivity (experiment 04)

**Do not overwrite the established results from model 01.** Experiment 04 creates separate outputs, and uses the same δ¹⁸O/Cl mixing and aquifer kinetics as 01. It fits six models (three river-entry priors × zero- and first-order aquifer kinetics):

| Scenario | River-entry DO |
|---|---|
| `A_fixed_7p5` | Fixed 7.5 mg/L |
| `B_prior_sd1` | Normal(7.5, 1.0²) mg/L, truncated to 0–14.65 mg/L |
| `C_prior_sd1p5` | Normal(7.5, 1.5²) mg/L, truncated to 0–14.65 mg/L |

**Site 2801 is observed at 10.2 mg/L, not 10.5 mg/L.** `00_setup.R` no longer substitutes 7.5 mg/L for its DO. The baseline 01 fit **masks** 2801 DO as before; in experiment 04, the default `include_proxy_DO=TRUE` includes the actual 2801 observation as ONE river-proximal groundwater sample for **all six scenario fits**, rather than fixing the river boundary to that value. No source file is modified. The reported zero-MRT age at 2801 remains fixed at zero, which is an additional assumption requiring sensitivity analysis. To repeat all six scenarios without using 2801 DO, change `cfg04$include_proxy_DO=FALSE`; this writes to a separate folder.

Run after pulling the latest changes:

```r
source("scripts/04_river_DO_prior_sensitivity.R")
```

Results, individual fit objects, diagnostic reports and posterior correlations are written under `outputs/04_river_DO_prior/2801_included/` (or `2801_masked/`). Run `source("scripts/02_DO_plotting.R")` afterwards to generate the additional Turbo-colour boundary and rate sensitivity plots. **02 does not refit anything**.

The fixed 04 scenario is directly comparable to the uncertain 04 scenarios because the observations are identical. **Do not attribute changes between the old baseline 01 and 04 solely to the river prior**—the inclusion of 2801 DO differs. Posterior predictive tail probabilities in 04 are in-sample checks, not grouped cross-validation. A global boundary DO can mask seasonality and spatial variation; convergence and posterior stability must be checked before adopting the prior in the primary model.


### Grouped cross-validation of broad river DO prior (script 05)

To test the candidate **river-entry DO ~ Normal(7.5, 1.5²) mg/L**, use:

```r
source("scripts/05_river_DO_grouped_CV.R")
```


This is a **standalone CV runner**. It will NOT rerun full-data fits 01 or 04. It runs five matched groups of held-out DO for four variants:

| River boundary | Aquifer kinetics |
|---|---|
| Fixed 7.5 mg/L | Zero order |
| Fixed 7.5 mg/L | First order |
| Normal(7.5, 1.5²) truncated 0–14.65 mg/L | Zero order |
| Normal(7.5, 1.5²) truncated 0–14.65 mg/L | First order |

**Primary: 2801 excluded from the DO likelihood AND from scoring.** Its measured 10.2 mg/L is left intact in the data and its δ¹⁸O and Cl are retained for source mixing. All remaining **40 wells** receive one held-out DO prediction in each model. Site 2130 is **held out** when its fold is evaluated. All four variants use identical group folds and the same reported MRTs, endmember priors, DO observation model and aquifer-rate priors.

The first run performs **20 CV MCMC fits** (four variants × five folds), saves each fit object immediately and writes per-fit predictions/diagnostics. If interrupted, re-run the same command: previously completed fit objects are reused. The script checks a saved design fingerprint (data, source scripts, Stan model, options, fold assignments) to reject incompatible cached fits. For an initial single-fold test, set `cfg05$run_only_fold=1L` inside script 05; set it back to `NULL` for the complete study. **Do not rank models using incomplete folds or failed MCMC diagnostics.**

Outputs appear in `outputs/05_river_DO_grouped_CV/2801_masked/` including:

- `cv_predictions.csv`: held-out mean DO, Student-t posterior predictive 95% limits, LSR fraction and MRT.
- `cv_scores.csv`: RMSE, MAE, bias, predictive 95% coverage, fold counts and convergence flags.
- `cv_fit_status.csv`: R-hat, bulk and tail ESS, E-BFMI, divergences, treedepth.
- `cv_2130_predictions.csv`: genuinely held-out DO for site 2130.
- `cv_group_bootstrap_summary.csv`: paired **cluster bootstrap** (by named well family) of predictive differences, produced only once all 20 fits pass diagnostics. Negative **C-minus-A** favours the broad prior; negative **first-minus-zero** favours first-order kinetics.

**2801 inclusion sensitivity:** Set `cfg05$include_2801_DO=TRUE` in script 05 and rerun. This repeats the same 20 fits with 2801’s true DO **included in training only**, while always scoring the same 40 non-proxy wells. It writes to the distinct `outputs/05_river_DO_grouped_CV/2801_in_training/` folder. This separates the effect of including 2801 from that of changing the river prior.

For additional **Turbo plots**, run:

```r
source("scripts/02_DO_plotting.R")
```


Script 02 automatically picks up available CV outputs and creates plots 11 (held-out observed/predicted DO, 95% predictive intervals, 2130 labelled) and 12 (per-well change in absolute error for broad prior minus fixed). **No MCMC occurs in script 02.**

#### Reviewing CV diagnostics and selectively retrying a single fold

The original script 05 had a **false-negative `all_folds_complete` flag** caused by comparing integer sample indices with CSV-loaded double indices using `identical()`. This is now fixed by explicit integer-normalised set checks. `cv_scores.csv` now reports convergence **per scenario × kinetic order**, instead of marking all four models as failed if a single fold fails.

The first 20-fit run completed 5/5 folds and 40 predictions per model. Only **`C_prior_sd1p5` × first-order × fold 3** failed the additional tail-ESS threshold (**116 vs minimum 200**); it had R-hat 1.01, bulk ESS 431, and no divergent transitions. This does **not** imply poor predictive accuracy or invalidate the other 19 fits, but uncertainty tails for that fold need better sampling.

**No MCMC needed to correct the scores:** pull the updated script, then source `scripts/05_river_DO_grouped_CV.R` with the default `cfg05$refit_failed=FALSE`. It reuses all existing RDS fit objects, rewrites the corrected flags and scores and reports the single failing fold.

To **rerun only the failed fit**, edit the `cfg05` settings near the top of script 05:

```r
run_only_fold=3L,
resume=TRUE,
refit_failed=TRUE,
failed_sampling=3500L,
```

Then source script 05 again. This reads the existing fit/status files and resamples **only the failing model in fold 3**, leaving the other 19 fits untouched. The old fit is preserved at `*.previous_fail.rds`, and sampling-iteration counts are recorded in fit status. Afterwards, reset `run_only_fold=NULL` and `refit_failed=FALSE` for subsequent full report refreshes. The original 20-fold design fingerprint and existing fits remain compatible; changes to the input data or Stan model are not silently reused.

The four current RMSEs are **1.64 (fixed zero), 1.59 (fixed first), 1.57 (broad zero), 1.60 (broad first)** mg/L. These differences are small; do not claim a preferred kinetic order or a CV improvement for the broad river prior without the paired group-bootstrap uncertainty.

### Simple DO-site exclusion in grouped CV

The **only exclusion setting** for script 05 is at the top of
`scripts/05_river_DO_grouped_CV.R`:

```r
include_2801_DO = FALSE,
exclude_DO_sites = "2130", # omit 2130 DO from ALL training and scoring
```

This defaults to a **39-well regional CV**. It masks 2130's DO from every fit (not merely the scored predictions), but **retains its δ¹⁸O and chloride** in the conservative mixing likelihood. The original private CSV and the standard `00_setup.R` remain unchanged. The historical 40-well analysis is preserved.

Outputs are automatically isolated in `outputs/05_river_DO_grouped_CV/2801_masked_no2130/`. To return to the original 40-well setup, change the one setting to `exclude_DO_sites = NULL`; that selects the existing `2801_masked/` output directory and reuses the previous completed fits. All five folds are generated afresh for the 39 scored wells.

Script `02_DO_plotting.R` recognizes both paths, with and without 2130, automatically. With the default excluded-site configuration it labels plots accordingly and makes no claim to predict 2130. All fit/scoring diagnostics and grouped uncertainty comparisons remain separate across these two designs.

**Scientific rationale:** 2130 is a shallow water-race-adjacent observation in the far west that may not represent regional aquifer groundwater. This is a hydrogeological *representativeness exclusion*, not an attempt to manipulate fit residuals. The original 40-well model and its CV are retained for sensitivity reporting.


### Regional model posterior refits after excluding well 2130 (script 04)

**Primary 39-DO-site model:** Download/replace `scripts/04_river_DO_prior_sensitivity.R` from the latest repository version. At the top of script 04:

```r
include_proxy_DO=FALSE, # masks 2801's DO (10.2 mg/L)
exclude_DO_sites="2130", # masks 2130's DO
scenarios=tibble::tibble(
  scenario=c("A_fixed_7p5","C_prior_sd1p5"),
  infer_river_DO=c(0L,1L),
  river_DO_prior_sd=c(0,1.5)
),
```


There is **no** `include_2801_DO` setting in script 04; that setting belongs to script 05. The two excluded DO values are zeroed in Stan data and their likelihood indicators are off. **Both wells' δ¹⁸O and chloride data remain in the mixing model; neither is dropped from the private source CSV.** Their original observed DO is preserved in `d$DO_original` and reported as a diagnostic only.

This version runs four full-data posterior fits (fixed vs broad river prior, each with zero- and first-order kinetics). The earlier intermediate SD=1 prior may be restored for sensitivity, but adds two fits and previously required extra MCMC attention. Outputs are written to **`outputs/04_river_DO_prior/2801_masked_no2130/`**, leaving prior `2801_included/` and `2801_masked/` posterior fits untouched. Compare `comparison_summary.csv` and `fit_status.csv`. The primary selected specification is `A_fixed_7p5` with `aquifer_order=1`; the other three fits provide kinetic/boundary sensitivities.

Run:

```r
source("scripts/04_river_DO_prior_sensitivity.R")
```


The updated plotting-only script `scripts/02_DO_plotting.R` can additionally read this new `2801_masked_no2130/` directory, producing the river-boundary and parameter-comparison figures without refitting the models.
