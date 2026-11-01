# Heretaunga groundwater oxygen reduction — two-endmember DO model

## Primary analysis: `01_DO_dual_pathway_kinetics`

We model two recharge endmembers—**river recharge** and **land-surface recharge (LSR)**—using only **δ¹⁸O and chloride** to constrain their mixing. Dissolved oxygen (DO) is reactive: each recharge pathway loses oxygen during aquifer residence, and the resulting concentrations are then mixed. The two competing aquifer models are **zero-order** and **first-order** kinetics, fitted independently with uncertainty.

**Sodium and excess air are excluded**, following diagnostic experiments. The model contains no sodium or excess-air parameters. The former model and ablation experiment are preserved for reproducibility, not used in the primary analysis.

### Active files

| Path | Purpose |
|---|---|
| [`scripts/00_setup.R`](scripts/00_setup.R) | **Setup only**: parameters, private data validation, filtered wells, O18/Cl observations and `stan_data`; no MCMC |
| [`scripts/01_DO_dual_pathway_kinetics.R`](scripts/01_DO_dual_pathway_kinetics.R) | **Primary R runner**: zero and first order, diagnostics, outputs, optional grouped CV; sources setup automatically |
| [`stan/01_DO_dual_pathway_kinetics.stan`](stan/01_DO_dual_pathway_kinetics.stan) | **Primary Stan model**: δ¹⁸O/Cl mixing and dual recharge-pathway DO kinetics |
| [`docs/01_DO_dual_pathway_kinetics.md`](docs/01_DO_dual_pathway_kinetics.md) | Scientific protocol, assumptions, outputs, limitations |
| [`data/README.md`](data/README.md) | Private data requirements |

Earlier models `23` and `23a`, their corresponding Stan code and the former `00_check_setup.R` preflight are preserved unchanged in [`archive/`](archive/README.md), alongside older experiments. **Only scripts 00 and 01 and the matching Stan file are active.** Setup is separate from sampling. The earlier `archive/scripts/00_check_setup.R` is a different, obsolete preflight.

### Run from Positron

1. Open the **Oxygen-Reduction-Joint-Model** repository folder as the working directory.
2. Copy the unchanged private source `Uwe_HP_data_prepared.csv` into `data/`; it must not be committed.
3. Ensure R packages `tidyverse`, `posterior`, `cmdstanr` and CmdStan are installed.
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

The initial run fits **two full models only** (four chains per model). The model-runner saves each fitted posterior immediately before producing reports. Results appear in:

```text
outputs/01_DO_dual_pathway_kinetics/
plots/01_DO_dual_pathway_kinetics/
```

First inspect `fit_status.csv` and `sampler_chains.csv`. Then examine `global_posteriors.csv`, `well_posteriors.csv`, `identifiability.csv` and the joint-posterior plot. Do **not** interpret model parameters if full-fit diagnostics fail.

**Grouped held-out DO comparison:** After both full fits converge, set `run_CV = TRUE` in the top-level `cfg` list of `scripts/00_setup.R`. It runs five matched grouped folds for both kinetic orders (10 extra fits), saves CV diagnostics and predictions, and performs a matched group bootstrap only if all folds converge. Source-tracer observations are retained for the held-out wells; their DO observations are masked.

### Fixed scientific assumptions

- Rainfall starts at approximate atmospheric DO saturation at **15 °C** and an assumed atmospheric pressure ratio of 1.
- River-to-aquifer entry DO is **fixed at 7.5 mg/L**.
- River-proximal proxy **2801** is assigned 7.5 mg/L in memory only and its DO is **excluded** from likelihood and held-out scoring. Original private data are unchanged.
- Far-south anomalous well **16758** is excluded.
- LSR enters the water table with an inferred DO concentration between 0 and rainfall saturation. The result is an **LSR oxygen deficit**, not an identifiable LSR kinetic rate without independent vadose travel time.
- Positive measured MRTs receive non-centred, mean-preserving lognormal uncertainty; proxy 2801 has MRT=0. The model assumes MRT effectively represents aquifer residence time for both recharge sources.
- Zero-order rates are reported in **mg/L/year**; first-order rates in **year⁻¹**. These are alternative aquifer models, not simultaneous reactions.
- Arithmetic averages over wells are **not** flow-weighted catchment oxygen budgets.

A converged sampler demonstrates computational reliability, **not** structural identifiability or validation of the endmembers. Check grouped predictive performance and sensitivity to endmembers, MRT and boundaries before physical interpretation or publication.

**Code verification:** The new files were statically checked and aligned with the previously converged no-sodium/no-air model. Compilation and MCMC execution must be verified in the local R/CmdStan environment.

### Obsolete local `0_setup.R` (one zero)

The previous local `scripts/0_setup.R` is **not** the new `scripts/00_setup.R`. It was an incorrectly named copy of the MCMC runner (contained `cmdstan_model()` and `$sample()`). It is not in the current GitHub tree and cannot be removed remotely from a local working copy. Preserve it outside the repository if needed, but **do not run it**. The correct setup file has **two zeros**: `scripts/00_setup.R`.
