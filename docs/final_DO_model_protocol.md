# Final groundwater DO model — reproducible analysis protocol

**Status: implementation committed; CmdStan fits have NOT been executed or convergence-validated in this environment.** This is the final proposed analysis, not a numerical-results report.

## Scientific objective

Estimate (a) oxygen depletion between rainfall saturation and the water table for land-surface recharge (LSR) and (b) the **aquifer oxygen-reduction rate**, comparing **ZERO-order** with **FIRST-order** kinetics. Test the posterior trade-off between these quantities and their robustness to boundary conditions, tracer and MRT assumptions. Avoid interpreting age-correlated apparent DO loss as incontrovertible evidence of the reaction location.

## Physical boundaries and observations

| Boundary / datum | Assumption | Source or handling |
|---|---|---|
| Rain DO at land surface | Saturated at recharge temperature and elevation | Default 15 °C, approximate atmospheric pressure ratio 1; oxygen-saturation polynomial. Explore 10 and 20 °C. |
| LSR at water table | Unknown \(C_{\mathrm{LSR,WT}} = C_{\mathrm{sat}} s\), \(0<s<1\) | Infer \(s\) jointly from tracer mixing, age and DO; default \(s \sim \mathrm{Beta}(1,1)\). |
| River recharge at aquifer entry | **Exactly 7.5 mg/L for primary model** | User-chosen defensible boundary; explore fixed 6.5 and 8.5 mg/L ONLY as separately labelled sensitivities. |
| River-dominated study well 2801 | **Model-only DO assigned 7.5 mg/L** | Original measured DO (10.2 mg/L) is retained unmodified in the private source CSV. Model-only assigned DO is **masked out of the fitting likelihood** and excluded from CV, avoiding circular pseudo-evidence. Its tracer data remain. |
| Aquifer age | Reported MRT is an uncertain proxy for aquifer residence time | Lognormal mean-centred prior, log SD 0.25 baseline, 0.5 sensitivity. No separate LSR-specific MRT inferred. |

The original private `data/Uwe_HP_data_prepared.csv` contains 43 records. Exclude the far-south outlier **16758**: 42 remain. Of these, 41 have recorded DO and MRT, and 40 *independent DO responses* remain after masking 2801. No source observation is overwritten in the original CSV.

Mixing tracers: δ18O, Cl and Na; Na–Cl residual dependence explicitly modelled. Include excess air **only as a sensitivity**. Tracer-endmember prior centres are partially study-derived, requiring a wide-prior sensitivity. LSR δ18O time/spatial variability contributes additional source-fraction-dependent uncertainty.

DO is the **sole reactive response** in this final primary fit. Earlier paired excess-N₂ and denitrification-progress alternatives did not reliably improve held-out DO prediction. Neither is treated as a conservative tracer or stoichiometric DO constraint.

## Aquifer kinetics: two alternatives

For *each source* before final linear mixing:

- **Zero order:** \(C(t)=\max[0,C(0)-k_0 t]\), \(k_0\) in **mg L⁻¹ yr⁻¹**.
- **First order:** \(C(t)=C(0)\exp(-k_1 t)\), \(k_1\) in **yr⁻¹**.

The source concentrations \(C(0)\) are \(7.5\) mg/L for river recharge and \(C_{\mathrm{LSR,WT}}\) for LSR. The **same** modelled aquifer MRT applies to both source pathways, an assumption that must be acknowledged.

Both models share the weakly informative half-normal prior on fractional depletion at \(C_\mathrm{ref}=7.5\) mg/L, sampled as \(k_\mathrm{ref}\ge0\):
\(k_0=7.5\,k_\mathrm{ref}\) versus \(k_1=k_\mathrm{ref}\), with \(k_\mathrm{ref}\sim\mathrm{Normal}^+(0,0.15)\) yr⁻¹. This is a comparable reference-relative prior, not proof of equal model complexity or mechanism.

**Primary identifiable LSR quantity:** \(\Delta C_{\mathrm{LSR}}=C_{\mathrm{sat}}-C_{\mathrm{LSR,WT}}\) mg/L, including full posterior uncertainty. The data do **not** identify separate unsaturated-zone zero- versus first-order rate constants absent vadose travel-time information. *Conditional* rates for assumed 0.25, 1 and 5-year residence times are only derived transformations of the same inferred DO deficit and must be labelled as such.

A robust Student-t(df=4) DO likelihood prevents a single anomalous measurement from dominating the common k, but does not correct missing hydrochemical processes.

## Run in Positron

From the repository root, with CmdStanR, CmdStan, posterior and tidyverse installed, and the original CSV copied to the local gitignored \`data/\` folder:

```r
source("scripts/23_final_DO_two_pathway_kinetics.R")
```

Current defaults: 4 chains, 1400 warmup and 1400 post-warmup draws per chain; **30 Stan fits**: 2 full models, 2×5 grouped held-out folds, and 2×9 prior/boundary sensitivity fits. This could take hours. Use Positron on a suitable machine. To smoke-test, edit settings at the top of the runner to set `run_prior_sensitivity=FALSE` and `run_only_fold=1L` before executing (4 fits). For full posteriors only, set both `run_CV=FALSE` and `run_prior_sensitivity=FALSE` (2 fits). Restore both switches before final reporting. Each full fit is saved as an RDS (unless `save_fits=FALSE`).

## Mandatory outputs and decision checks

| Local output | Question addressed |
|---|---|
| \`outputs/final_DO_23_boundary_conditions.csv\` | Are rainfall and fixed river boundary DO values exactly what we intend? |
| \`outputs/final_DO_23_global_posteriors.csv\` | Zero- or first-order aquifer rate (correct units), LSR DO deficit, water-table DO, credible intervals, per-well average nominal loss components. |
| \`outputs/final_DO_23_full_diagnostics.csv\`, \`outputs/final_DO_23_fit_status.csv\` | Any problematic Rhat, low bulk ESS, divergences, treedepth saturation, low E-BFMI? |
| \`outputs/final_DO_23_identifiability_summary.csv\` | How strongly does posterior LSR DO deficit trade off against aquifer reaction rate? How much has prior LSR uncertainty narrowed? |
| \`plots/final_DO_23_identifiability.png\`, \`outputs/final_DO_23_joint_posterior_draws.csv\` | Shape, direction, multimodality or nonlinearity of the joint posterior trade-off. |
| \`outputs/final_DO_23_cv_scores.csv\`, \`outputs/final_DO_23_cv_bootstrap_summary.csv\` | Which order predicts *the same held-out wells* better? Does a well-group bootstrap show that advantage is robust? |
| \`outputs/final_DO_23_cv_fit_status.csv\` | Are the held-out refits sufficiently converged for comparison? |
| \`outputs/final_DO_23_sensitivity_posteriors.csv\`, \`outputs/final_DO_23_sensitivity_correlations.csv\`, \`outputs/final_DO_23_sensitivity_diagnostics.csv\` | Do rate and LSR-loss conclusions survive rainfall temperature, river boundary, MRT uncertainty, LSR prior, widened tracer prior and excess-air sensitivity? |
| \`outputs/final_DO_23_conditional_LSR_rates.csv\` | Optional *assumed-time-dependent* LSR zero/first-order rates (NOT identified from well DO alone). |

**Do not select an aquifer model if either the full or held-out fits fail sampling diagnostics.** A descriptive matched group bootstrap interval of the CV RMSE difference crossing zero does not support a strong predictive ranking. Even if an order is better, inspect the **posterior rate-versus-LSR-deficit dependence** and sensitivity to the fixed river DO and MRT assumptions before attributing reaction location or rates physically.

The arithmetic average nominal DO losses across sampled wells are **not a flow-weighted catchment oxygen budget**. Uncertainty in groundwater age distributions and source-specific residence times remains outside the present identifiable model. New observations at/near the LSR water table are probably the most effective additional constraints if the trade-off remains high.
