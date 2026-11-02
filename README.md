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
| [`06_DO_posterior_distributions.R`](scripts/06_DO_posterior_distributions.R) | Posterior density and joint parameter plots using the saved 04 CmdStanR RDS fits | **No** |

Recommended sequence **for a new run** (each `source()` is a separate Positron R-console command):

```r
source("scripts/00_setup.R")                        # optional check, no MCMC
source("scripts/04_river_DO_prior_sensitivity.R")    # fit four full-data models
source("scripts/05_river_DO_grouped_CV.R")           # CV; reuses valid cached fits
source("scripts/02_DO_plotting.R")                   # plots only
source("scripts/06_DO_posterior_distributions.R")   # plots only
```

**If the saved 04 and 05 fits and reports already exist, run only scripts 02 and 06**. No more calibration iterations are needed: all four final full-data fits passed the sampler checks, and posterior bulk/tail effective sample sizes for the primary parameters exceeded 3,000.

### Model configuration and scientific findings

- River DO entering the aquifer: primary **7.5 mg/L fixed** (`A_fixed_7p5`), sensitivity prior **Normal(7.5, 1.5²)** mg/L truncated to [0, 14.65] (`C_prior_sd1p5`).
- LSR DO begins near freshwater saturation at 15°C; its water-table oxygen content is inferred. Both pathways share either zero- or first-order aquifer consumption kinetics.
- The preferred regional model (based on grouped CV) is **fixed river DO + first-order kinetics**. First-order aquifer reaction rate: **0.0625 yr⁻¹** (95% credible interval **0.0484–0.0812**).
- LSR water-table DO: **3.28 mg/L** (95% CrI **1.42–5.47**); inferred LSR DO deficit: **6.72 mg/L** (95% CrI **4.53–8.58**).
- Posterior rate is stable under a broad river DO prior (~0.0632 yr⁻¹). The LSR deficit remains model-dependent: it does **not** uniquely identify the spatial or physical location of oxygen loss.
- 39-well grouped CV favored first-order kinetics over zero-order with fixed river DO: median paired ΔRMSE (first − zero) **−0.132 mg/L**, 95% grouped-bootstrap interval **[−0.238, −0.0126]**. Broader river DO did not improve first-order CV RMSE.

### Output locations

- Full posterior fits/diagnostics: `outputs/04_river_DO_prior/2801_masked_no2130/`.
- Grouped CV results/bootstrap: `outputs/05_river_DO_grouped_CV/2801_masked_no2130/`.
- CV and parameter-summary plots: `plots/final_39well/`.
- Marginal/joint posterior plots: `plots/04_river_DO_prior/2801_masked_no2130/posterior_distributions/`.

`outputs/`, `plots/` and private `data/*.csv` are gitignored. **GitHub does not contain the user's private measurements or local saved MCMC fits.** The cleanup does not delete or modify these directories locally.

### Code and provenance

The active Stan model is [`stan/04_river_DO_prior_sensitivity.stan`](stan/04_river_DO_prior_sensitivity.stan). [`stan/01_DO_dual_pathway_kinetics.stan`](stan/01_DO_dual_pathway_kinetics.stan) is retained **solely because the unchanged shared setup checks that it exists**. Do not edit `scripts/00_setup.R` simply to remove that historical dependency: script 05 fingerprints its exact contents to protect the existing saved CV fit cache.

Research history, the old 01 baseline runner, 03 high-DO audit, the original script 02 plotting workflow, and superseded protocol text are safely preserved in [`archive/`](archive/README.md), not in the active `scripts/` folder. The archived scripts are **not** part of the current pipeline.

The first-order model choice reflects relative predictive performance; acceptable MCMC diagnostics establish numerical convergence, not hydrological model truth. See [`docs/final_39well_model.md`](docs/final_39well_model.md) for scientific interpretation and limitations.
