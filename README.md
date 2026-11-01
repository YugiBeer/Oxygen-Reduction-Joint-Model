# Heretaunga groundwater oxygen reduction — current workflow

**Start here.** This repository now has one active primary analysis: the final Bayesian comparison of **zero-order versus first-order aquifer oxygen reduction**, with an uncertain land-surface recharge (LSR) oxygen deficit.

Earlier experiments and exploratory scripts **01–22** are retained, unchanged, under [`archive/`](archive/README.md). You **do not** have to run scripts 01–22 before script 23. The full historical README is retained at [`archive/README_legacy.md`](archive/README_legacy.md).

## Active files

| File | Role |
|---|---|
| [`scripts/00_check_setup.R`](scripts/00_check_setup.R) | Check local private input, required columns, packages and CmdStan before fitting |
| [`scripts/23_final_DO_two_pathway_kinetics.R`](scripts/23_final_DO_two_pathway_kinetics.R) | **Primary model runner**: zero- and first-order aquifer kinetics |
| [`scripts/23a_sodium_ablation_diagnostics.R`](scripts/23a_sodium_ablation_diagnostics.R) | Optional matched **with/without sodium** full-fit convergence diagnostic (four Stan fits) |
| [`stan/23_final_DO_two_pathway_kinetics.stan`](stan/23_final_DO_two_pathway_kinetics.stan) | Stan model used by the runner |
| [`docs/final_DO_model_protocol.md`](docs/final_DO_model_protocol.md) | Scientific assumptions, diagnostics, uncertainty and interpretation |
| [`archive/README.md`](archive/README.md) | Guide to earlier studies and recovering archived scripts |

## Run on your computer (Positron)

1. Open the **Oxygen-Reduction-Joint-Model** repository folder as the Positron project / working directory. Check `getwd()`: `scripts/23_final_DO_two_pathway_kinetics.R` and `stan/23_final_DO_two_pathway_kinetics.stan` must be visible relative to it.
2. Place the original **private** `Uwe_HP_data_prepared.csv` in `data/Uwe_HP_data_prepared.csv`. It contains sodium; the older reduced `groundwater.csv` is **not** the correct input. Research data remain local and must **never** be committed.
3. Have R, `tidyverse`, `posterior`, `cmdstanr` and **CmdStan** installed. If needed, see [CmdStanR installation](https://mc-stan.org/cmdstanr/articles/cmdstanr.html).
4. In the R console, run:

```r
source("scripts/00_check_setup.R")
source("scripts/23_final_DO_two_pathway_kinetics.R")
```

Run those commands **from the R console**, not by inserting a `source()` of the runner inside the runner itself.

### Stage A — first run (default; 2 Stan fits)

The runner now defaults to:

```r
run_CV=FALSE,
run_prior_sensitivity=FALSE,
```

This fits **one zero-order and one first-order model**, using four chains for each. It does not depend on results from any prior scripts. Inspect:

- `outputs/final_DO_23_boundary_conditions.csv` — fixed rainfall and river entry boundary conditions
- `outputs/final_DO_23_global_posteriors.csv` — posterior rates (with **different units** for each reaction order) and LSR DO deficit
- `outputs/final_DO_23_fit_status.csv` and `outputs/final_DO_23_full_diagnostics.csv` — convergence and sampler health
- `outputs/final_DO_23_identifiability_summary.csv` and `plots/final_DO_23_identifiability.png` — trade-off between LSR loss and aquifer rate

**Do not infer the better reaction order from two in-sample fits.** Check diagnostics first.

### Optional current diagnostic — sodium inclusion versus omission

Previous full fits had severe R-hat / ESS and E-BFMI problems, most evident in Na–Cl residual and source-fraction parameters, even after non-centred MRT reparameterisation. **Before cross-validation, run a controlled sodium-ablation diagnostic**:

```r
source("scripts/23a_sodium_ablation_diagnostics.R")
```

This makes **four separate full fits**: zero and first order, each with Na included or completely masked. It uses the unchanged Stan model, all the same wells/DO responses, MRT prior, boundaries, priors and sampler settings. Excess air is off in both, so **without Na the retained mixing tracers are δ18O and Cl only**. No reactive nitrogen is fitted.

Read `outputs/sodium_ablation_23a_comparison.csv` for worst active-parameter R-hat, minimum bulk ESS, divergences, treedepth and per-chain E-BFMI; `..._core.csv` and `..._chains.csv` for key posterior parameters and per-chain medians. Full per-parameter diagnostics and four RDS fits are also saved under `outputs/sodium_ablation_23a_*`. These results **do not overwrite** primary `final_DO_23_*` outputs.

**Important qualification:** Sodium-specific nuisance parameters `Na_R`, `Na_LSR`, `sigma_Na`, and `rho_Na_Cl` remain in Stan but are sampled from their priors only when Na observations are masked. They are labelled in the diagnostics and excluded from the ablation's *active-parameter* convergence aggregate, but their per-parameter diagnostics are still saved. The excess-air nuisance parameters are likewise prior-only in both variants. A convergence improvement would isolate the contribution of the sodium likelihood to posterior geometry; it **would not establish** that sodium is chemically unsuitable or validate source mixing. Interpret no kinetic estimates from fits that fail sampler diagnostics.

Script 23 also now exposes `include_sodium=TRUE` as a switch for later controlled work, but manually turning it off would reuse its standard output filenames; **prefer 23a for matched comparison**.

### Stage B — held-out prediction (12 fits total when run)

Once Stage A runs acceptably, change `cfg` at the top of `scripts/23_final_DO_two_pathway_kinetics.R`:

```r
run_CV=TRUE,
run_prior_sensitivity=FALSE,
```

Re-running fits both full models **again**, plus the **10 grouped-CV refits**. Examine `outputs/final_DO_23_cv_scores.csv`, `outputs/final_DO_23_cv_fit_status.csv` and `outputs/final_DO_23_cv_bootstrap_summary.csv`. To test only one fold first, set `run_only_fold=1L`, but do **not** use single-fold scores to rank kinetic orders.

### Stage C — boundary/prior sensitivities (30 fits total)

For a full final audit, set `run_CV=TRUE`, `run_only_fold=NULL` and `run_prior_sensitivity=TRUE`. This runs 2 full models, 10 CV fits and 18 sensitivity fits. Sensitivities include recharge temperature, fixed river DO, MRT uncertainty, LSR survival prior, mixing-endmember prior width and excess air. Check `outputs/final_DO_23_sensitivity_*.csv` before interpreting physical rates.

**The script was statically reviewed but has not yet been executed under R/CmdStan in this editing environment.** Your local full-fit and CV diagnostics determine whether any posterior is suitable for scientific interpretation.

## Fixed boundary conditions and essential caveats

The primary river-to-aquifer DO boundary is **fixed at 7.5 mg/L**. The river-proximal proxy well **2801** receives 7.5 mg/L only in the *model data*, and its DO response is **masked from both fitting and CV**. The original private observation is never overwritten. Outlier **16758** is excluded. Rainfall starts at assumed atmospheric DO saturation (default 15 °C, approximate pressure ratio 1); DO at the LSR water table is inferred.

The model estimates a *surface-to-water-table LSR DO deficit*, **not independently identifiable unsaturated-zone reaction rates**. Conditional LSR rates in the output depend on assumed vadose travel times. The aquifer rate depends on the reported MRT representing aquifer residence and on the tracer mixing formulation. Chloride/mixing sensitivity from script 12 remains a warning against claiming that LSR loss is precisely determined; script 23's correlated Na–Cl observation model does **not** by itself resolve all source-fraction uncertainty. See [the scientific protocol](docs/final_DO_model_protocol.md).

## Keeping a clean local copy

For an existing local clone, first check `git status` and preserve any local edits. Then run `git pull --ff-only` from the repository directory. This will move previously tracked scripts into `archive/`, but it will not remove ignored local data and outputs. If the pull cannot fast-forward because of local edits, resolve or stash those edits **before** pulling; do not force-reset unless you deliberately want to discard local changes.

All `data/*.csv`, `outputs/`, and `plots/` are ignored by Git. The archive is version-controlled and recoverable.
