# Heretaunga groundwater oxygen reduction — two-endmember DO model

## Primary analysis: `01_DO_dual_pathway_kinetics`

We model two recharge endmembers—**river recharge** and **land-surface recharge (LSR)**—using only **δ¹⁸O and chloride** to constrain their mixing. Dissolved oxygen (DO) is reactive: each recharge pathway loses oxygen during aquifer residence, and the resulting concentrations are then mixed. The two competing aquifer models are **zero-order** and **first-order** kinetics, fitted independently with uncertainty.

**Sodium and excess air are excluded**, following diagnostic experiments. The model contains no sodium or excess-air parameters. The former model and ablation experiment are preserved for reproducibility, not used in the primary analysis.

### Active files

| Path | Purpose |
|---|---|
| [`scripts/00_check_setup.R`](scripts/00_check_setup.R) | Verify project root, private input, R packages and CmdStan |
| [`scripts/01_DO_dual_pathway_kinetics.R`](scripts/01_DO_dual_pathway_kinetics.R) | **Primary R runner**: zero and first order, diagnostics, outputs, optional grouped CV |
| [`stan/01_DO_dual_pathway_kinetics.stan`](stan/01_DO_dual_pathway_kinetics.stan) | **Primary Stan model**: δ¹⁸O/Cl mixing and dual recharge-pathway DO kinetics |
| [`docs/01_DO_dual_pathway_kinetics.md`](docs/01_DO_dual_pathway_kinetics.md) | Scientific protocol, assumptions, outputs, limitations |
| [`data/README.md`](data/README.md) | Private data requirements |

Former model `scripts/23_final_DO_two_pathway_kinetics.R`, `stan/23_final_DO_two_pathway_kinetics.stan` and `scripts/23a_sodium_ablation_diagnostics.R` remain as **historical experiments**. Scripts 01–22 from an earlier experimental series are in [`archive/`](archive/README.md).

### Run from Positron

1. Open the **Oxygen-Reduction-Joint-Model** repository folder as the working directory.
2. Copy the unchanged private source `Uwe_HP_data_prepared.csv` into `data/`; it must not be committed.
3. Ensure R packages `tidyverse`, `posterior`, `cmdstanr` and CmdStan are installed.
4. Run:

```r
source("scripts/00_check_setup.R")
source("scripts/01_DO_dual_pathway_kinetics.R")
```

The initial run fits **two full models only** (four chains per model). The model-runner saves each fitted posterior immediately before producing reports. Results appear in:

```text
outputs/01_DO_dual_pathway_kinetics/
plots/01_DO_dual_pathway_kinetics/
```

First inspect `fit_status.csv` and `sampler_chains.csv`. Then examine `global_posteriors.csv`, `well_posteriors.csv`, `identifiability.csv` and the joint-posterior plot. Do **not** interpret model parameters if full-fit diagnostics fail.

**Grouped held-out DO comparison:** After both full fits converge, set `run_CV = TRUE` in the top-level `cfg` list of `scripts/01_DO_dual_pathway_kinetics.R`. It runs five matched grouped folds for both kinetic orders (10 extra fits), saves CV diagnostics and predictions, and performs a matched group bootstrap only if all folds converge. Source-tracer observations are retained for the held-out wells; their DO observations are masked.

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
