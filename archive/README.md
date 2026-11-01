# Archived research workflows (not the active pipeline)

The active workflow is [script 01](../scripts/01_DO_dual_pathway_kinetics.R) and its matching [Stan model](../stan/01_DO_dual_pathway_kinetics.stan), using δ¹⁸O and chloride only. Script 01 performs its own preflight checks. **No archived script must be run first.**

This directory retains **original, unchanged code snapshots** for reproducibility of the research development. The original 96-kB, experiment-by-experiment README is preserved as [`README_legacy.md`](README_legacy.md), including historical numerical prototypes, caveats and instructions.

## Archive index

| Scripts | Historical question |
|---|---|
| 01–03 | Three-tracer endmember mixing, initial DO partition and figures |
| 04–09 | Bayesian recharge DO, confounding, river DO prior, LSR losses and residual tests |
| 10–11 | Groundwater age distributions, gamma/piston assumptions and Bayesian age priors |
| 12–14 | Joint Bayesian mixing and DO, posterior partition, grouped validation |
| 15–17 | Tracer-only held-out checks, pairwise mixing and Na/Cl ratio sensitivity |
| 18–19 | Sodium as additional mixing tracer; joint DO three-vs-four-tracer comparison |
| 20–22 | Paired DO/excess N₂, excess-air ablation and denitrification-progress alternatives |
| 23 | Full joint Na–Cl mixing plus zero-/first-order aquifer DO kinetics; superseded by script 01 |
| 23a | Matched with/without sodium convergence experiment establishing why the simplified script 01 was selected |
| 00_check_setup | Former preflight for script 01; archived because script 01 validates its own environment and data |

Corresponding Stan files (12, 15, 18, 19, 20, 22, 23) are stored under `archive/stan/`.

The former experimental protocol `docs/final_DO_model_protocol.md` is preserved unchanged at [`archive/docs/final_DO_model_protocol.md`](docs/final_DO_model_protocol.md). The active protocol is [`docs/01_DO_dual_pathway_kinetics.md`](../docs/01_DO_dual_pathway_kinetics.md). Files `23_final_DO_two_pathway_kinetics.R`, `23a_sodium_ablation_diagnostics.R` and `00_check_setup.R` were moved **unchanged** from `scripts/` to `archive/scripts/`.

**Important:** Archiving moved these scripts without rewriting their internal relative paths or dependencies. They are **historical snapshots, not directly runnable at their archived paths**. To reproduce a historical experiment, copy the relevant `archive/scripts/*.R` and `archive/stan/*.stan` back to the original `scripts/` and `stan/` paths in a separate working branch or clone, then follow `README_legacy.md`. Some archived scripts source other archived scripts. Do not mix archived model outputs, sample definitions or rate assumptions with script 01 without a deliberate comparison.

The active Stan file has no dependence on archived Stan files. Original code is also permanently available in Git history.
