# Heretaunga Plains — final two-source Bayesian groundwater DO model

**Current analysis: 39 DO observations**, excluding sites 2801 and 2130 from the dissolved-oxygen likelihood, but retaining their δ¹⁸O and chloride measurements. Site 16758 is excluded from all analysis. The original DO measurements are not altered in the private input CSV.

## Run this final workflow in Positron

Open the repository as your R working directory and ensure the private input file exists at `data/Uwe_HP_data_prepared.csv`.

| Active script | What it does | MCMC? |
|---|---|---|
| [`00_setup.R`](scripts/00_setup.R) | Load/check private data; prepare tracer, age and DO inputs | **No** |
| [`04_river_DO_prior_sensitivity.R`](scripts/04_river_DO_prior_sensitivity.R) | Full-data posterior fits: fixed/broad river DO × zero-/first-order aquifer kinetics | **Yes — 4 fits** |
| [`05_river_DO_grouped_CV.R`](scripts/05_river_DO_grouped_CV.R) | Matched 5-fold held-out DO comparison, with diagnostics and paired group bootstrap; reuses saved fold fits | **Yes — only if necessary** |
| [`02_DO_plotting.R`](scripts/02_DO_plotting.R) | Figures for final 04 parameter summaries and 05 CV (no old 01 results required) | **No** |
| [`06_DO_posterior_distributions.R`](scripts/06_DO_posterior_distributions.R) | **One** preferred-model first-order rate posterior plot; verifies posterior summaries directly from the saved fit | **No** |
| [`07_DO_LSR_MRT_sensitivity.R`](scripts/07_DO_LSR_MRT_sensitivity.R) | One-at-a-time robustness check: recharge temperature, LSR-survival prior and MRT uncertainty; reuses saved broad-prior first-order baseline | **Yes — 7 additional fits, cached** |

Recommended sequence **for a new run** (each `source()` is a separate Positron R-console command):

```r
source("scripts/00_setup.R")                        # optional check, no MCMC
source("scripts/04_river_DO_prior_sensitivity.R")    # fit four full-data models
source("scripts/05_river_DO_grouped_CV.R")           # CV; reuses valid cached fits
source("scripts/02_DO_plotting.R")                   # plots only
source("scripts/06_DO_posterior_distributions.R")   # one preferred-rate posterior figure
source("scripts/07_DO_LSR_MRT_sensitivity.R")        # robustness fits; cached if complete
```

**If the saved 04 and 05 fits and reports already exist, run only scripts 02 and 06 for normal plotting.** Script 07 is a *separate sensitivity study*: the user has completed all seven new fits and reported passing diagnostics. It tests scientific assumptions, not whether the original chains were long enough. Do not rerun 04 or 05. All four original full-data fits passed diagnostics; posterior bulk/tail ESS for key first-order parameters exceeded 3,000.

### Model configuration and scientific findings

- **Primary physical interpretation:** river-entry DO has **Normal(7.5, 1.5²)** mg/L prior, truncated to [0,14.65] (`C_prior_sd1p5`), together with first-order aquifer kinetics. Fixed river DO (`A_fixed_7p5`) is retained as a predictive comparator/sensitivity formulation.
- LSR DO begins near freshwater saturation at 15°C; its water-table oxygen content is inferred. Both pathways share either zero- or first-order aquifer consumption kinetics.
- The selected **physical primary model** is **uncertain river-entry DO + first-order kinetics**, with a reported rate of **0.0632 yr⁻¹** (95% CrI **0.0453–0.0851**) in the earlier script-04 summary. The latest script-07 *baseline drawn from the saved fit* gives **0.0639 yr⁻¹** (95% CrI **0.0458–0.0857**): use script 06 to verify the saved posterior and resolve this small summary discrepancy before publication. The **fixed-boundary first-order** alternative performed slightly better under grouped CV: rate **0.0625 yr⁻¹** (95% CrI **0.0484–0.0812**).
- Primary uncertain-boundary LSR water-table DO: **3.21 mg/L** (95% CrI **1.25–5.55**); inferred LSR deficit **6.79 mg/L** (95% CrI **4.45–8.75**).
- The aquifer rate is similar with fixed versus uncertain river DO boundaries. In script-07's **completed** seven-scenario robustness experiment, changing LSR recharge temperature from 15°C to 10/20°C shifted median rate by <0.2%; changing LSR-survival priors shifted it -1.0% to +3.4%; and MRT log-SD 0.10/0.50 shifted it -4.55% / +18.5% relative to the 0.25 baseline. All eight scenario/baseline sampler checks passed. These are one-at-a-time conditional sensitivities, **not** a test of source-specific travel times.
- Baseline script-07 water-table LSR DO median **3.26 mg/L** (95% CrI **1.31–5.57**), deficit **6.78 mg/L** (95% CrI **4.46–8.72**). At 10/20°C, inferred WT DO stayed ~3.25/3.23 mg/L, while the surface-to-WT deficit shifted to 8.02/5.79 mg/L. This difference is partly definitional, since the deficit is referenced to assumed surface saturation.
- 39-well grouped CV favored first-order kinetics over zero-order with fixed river DO: median paired ΔRMSE (first − zero) **−0.132 mg/L**, 95% grouped-bootstrap interval **[−0.238, −0.0126]**. Broader river DO did not improve first-order CV RMSE.

### Output locations

- Full posterior fits/diagnostics: `outputs/04_river_DO_prior/2801_masked_no2130/`.
- Grouped CV results/bootstrap: `outputs/05_river_DO_grouped_CV/2801_masked_no2130/`.
- CV and parameter-summary plots: `plots/final_39well/`.
- **One preferred-rate posterior plot** and posterior-draw summary CSV: `plots/04_river_DO_prior/2801_masked_no2130/posterior_distributions/`. The previous seven-figure plotter is preserved in `archive/scripts/06_DO_posterior_distributions_legacy.R` (not in the active workflow).
- Script 07 sensitivity summaries and fits: `outputs/07_DO_LSR_MRT_sensitivity/2801_masked_no2130/`; figures: `plots/07_DO_LSR_MRT_sensitivity/2801_masked_no2130/`.

`outputs/`, `plots/` and private `data/*.csv` are gitignored. **GitHub does not contain the user's private measurements or local saved MCMC fits.** The cleanup does not delete or modify these directories locally.

### Code and provenance

The active Stan model is [`stan/04_river_DO_prior_sensitivity.stan`](stan/04_river_DO_prior_sensitivity.stan). [`stan/01_DO_dual_pathway_kinetics.stan`](stan/01_DO_dual_pathway_kinetics.stan) is retained **solely because the unchanged shared setup checks that it exists**. Do not edit `scripts/00_setup.R` simply to remove that historical dependency: script 05 fingerprints its exact contents to protect the existing saved CV fit cache.

Research history, the old 01 baseline runner, 03 high-DO audit, the original script 02 plotting workflow, and superseded protocol text are safely preserved in [`archive/`](archive/README.md), not in the active `scripts/` folder. The archived scripts are **not** part of the current pipeline.

The first-order model choice reflects relative predictive performance; choosing the uncertain river boundary as primary is a physical-modelling decision despite a small CV advantage of the fixed boundary. Acceptable MCMC diagnostics establish numerical convergence, not hydrological model truth. See [`docs/final_39well_model.md`](docs/final_39well_model.md) for scientific interpretation and limitations.

### Known structural limitation

All fitted models assign **the same well-specific effective aquifer MRT** to the river and land-surface recharge fractions. The inferred 39-well mixing fractions do **not** specify where along each flowpath LSR entered. If younger LSR enters progressively, the common-MRT assumption can bias the effective rate and the inferred LSR oxygen deficit; the sign is not established when reported MRT may represent mixed-age groundwater. Sensitivity to MRT *uncertainty* does not solve this structural problem.

### Repository audit, November 2026

The active working scripts are intentionally limited to **00, 02, 04, 05, 06, 07**; historical development code is in `archive/`. Scripts 04/05 and the Stan model are left untouched except for explanatory comments: altering model inputs, setup, or Stan would invalidate the existing script-05 CV fit design checksum. Script 06 is a lightweight **one-figure** plotter; its seven-figure predecessor is archived. Plots and fit data are local-only; the numerical values recorded here come from the user's reported results, not a re-execution of Stan on GitHub.
