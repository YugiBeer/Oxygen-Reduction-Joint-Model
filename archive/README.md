# Archived research workflows (not the active pipeline)

The active workflow is [script 23](../scripts/23_final_DO_two_pathway_kinetics.R) and its corresponding [Stan model](../stan/23_final_DO_two_pathway_kinetics.stan). **No archived script must be run first.**

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

Corresponding Stan files (12, 15, 18, 19, 20, 22) are stored under `archive/stan/`.

**Important:** Archiving moved these scripts without rewriting their internal relative paths or dependencies. They are **historical snapshots, not directly runnable at their archived paths**. To reproduce a historical experiment, copy the relevant `archive/scripts/*.R` and `archive/stan/*.stan` back to the original `scripts/` and `stan/` paths in a separate working branch or clone, then follow `README_legacy.md`. Some archived scripts source other archived scripts. Never mix their output assumptions or sample counts with script 23 without a deliberate comparison.

The active Stan file has no dependence on archived Stan files. Original code is also permanently available in Git history.
