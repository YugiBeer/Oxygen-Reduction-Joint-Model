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
