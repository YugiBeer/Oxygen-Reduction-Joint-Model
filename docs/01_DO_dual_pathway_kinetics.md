# 01 — Dual-pathway DO kinetics (two-endmember model)

## Scientific objective

Estimate the oxygen deficit accumulated along land-surface recharge **before the water table** and the oxygen reduction rate accumulated in the **aquifer**, using two recharge endmembers and the minimum defensible conservative mixing model. Compare **zero-order** and **first-order** aquifer reaction models. Do not infer separate zero- and first-order vadose-zone rates, because the vadose travel time is unknown.

## Source mixing

Let \(f_i\) be the LSR fraction in well \(i\), \(1-f_i\) the river-recharge fraction. For conservative tracer \(X\) in \(\{\delta^{18}O,\mathrm{Cl}\}\):

\[
X_i^{\mathrm{pred}}=(1-f_i)X_R+f_i X_{LSR}.
\]

All four endmember tracer concentrations have bounded, informative priors as established in experiment 23. Each well has its own \(f_i\sim\mathrm{Beta}(1,1)\). Isotope residual variance includes a fraction-dependent LSR source-uncertainty term; chloride residual variance includes a measurement-error floor and an inferred residual SD. There is **no Na–Cl covariance**, sodium observation, excess-air observation or prior-only sodium/air parameter. This is the same *active* likelihood as the converged `noNa` variant of experiment 23a; removing nuisance priors should not change its marginal target distribution for common parameters.

## Oxygen boundary conditions and kinetics

- Assumed rainfall DO at land surface: \(C_{\mathrm{sat}}(15^\circ\mathrm{C})\), approximate atmospheric saturation at pressure ratio 1 (approximately 10.03 mg/L).
- River-recharge DO at aquifer entry: \(C_R=7.5\) mg/L, **fixed**.
- LSR at water table: \(C_{LSR,WT}=C_{\mathrm{sat}}s\), where \(s\sim\mathrm{Beta}(1,1)\) in the baseline.
- LSR pre-water-table oxygen deficit: \(\Delta C_{LSR}=C_{\mathrm{sat}}-C_{LSR,WT}\).
- Reported groundwater MRT is assumed to represent a shared aquifer residence time for both recharge sources.

For an entry concentration \(C\) and aquifer residence time \(t\), alternative aquifer models are:

\[
C_0(t)=\max(0,C-k_0 t),\qquad k_0=k_{\mathrm{ref}}C_{\mathrm{ref}},\;C_{\mathrm{ref}}=7.5\ {\rm mg/L}.
\]

\[
C_1(t)=C\exp(-k_1 t),\qquad k_1=k_{\mathrm{ref}}.
\]

The common prior is \(k_{\mathrm{ref}}\sim \mathrm{Normal}^+(0,0.15)\) with bounds [0, 2]. Apply the reaction law **separately to the river and LSR entry concentrations**, then mix surviving DO:

\[
DO_i^{\mathrm{pred}}
=(1-f_i)\,C_R(t_i)+f_i\,C_{LSR,WT}(t_i).
\]

The DO observation likelihood is Student-t with four degrees of freedom and inferred residual scale, including a 0.2 mg/L floor. Report the zero-order rate in mg/L/year and the first-order rate in year⁻¹.

## Groundwater ages

Reported positive MRTs have the unchanged non-centred lognormal uncertainty:

\[
z_i\sim N(0,1),\quad
t_i=\exp\{\log({\rm MRT}_i)-\tfrac{1}{2}\sigma_t^2+\sigma_t z_i\},
\quad \sigma_t=0.25.
\]

This preserves the reported MRT as the arithmetic mean of the age prior. The sole zero-MRT site, 2801, has fixed age zero and its DO observation is masked.

## Data inclusion

Private input: `data/Uwe_HP_data_prepared.csv`. Required columns: `Site_ID`, `Depth`, `DO`, `MRT`, `O18` and `Cl`. Site 16758 is excluded. Site 2801's original DO measurement (10.2 mg/L) is preserved both in the CSV and in memory. Baseline 01 masks its DO from fitting and scoring; experiment 04 can include it as a groundwater observation. With the currently prepared file, 41 wells remain after exclusion and basic DO/MRT filtering, 40 of which contribute DO likelihood observations.

## Execution and separation of responsibilities

`source("scripts/00_setup.R")` performs preparation **only**: loads required packages, sets all configuration values (`cfg`), checks the private input and CmdStan installation, filters the sites, computes oxygen boundary conditions, prepares the reported-MRT non-centred indexing and builds `stan_data`. It does **not** compile Stan, call `$sample()`, produce fit reports, or source the fitting script. This makes it safe to inspect `d`, `cfg`, `Csat`, `proxy_index` and `stan_data` before fitting.

`source("scripts/01_DO_dual_pathway_kinetics.R")` always runs setup again, then compiles/loads the existing `stan/01_DO_dual_pathway_kinetics.stan`, fits both kinetic orders and performs the same diagnostic and posterior-output steps as before. The Stan likelihood and all physical assumptions are unchanged by this refactor. To adjust `run_CV`, priors, boundary conditions or MCMC iterations, edit `cfg` in `scripts/00_setup.R`, **not** in script 01. Outputs are written only by script 01.

The untracked local `scripts/0_setup.R` (one zero) reported by Positron previously contained MCMC calls and is not part of the active workflow. It must not be confused with the new `scripts/00_setup.R` (two zeros).

## Plotting (script 02)

Run `source("scripts/02_DO_plotting.R")` **after** script 01. The plotting workflow never compiles a Stan model or performs MCMC. It sources `00_setup.R` to reuse the original filtered dataset and model settings, then reads the `global_posteriors.csv`, `well_posteriors.csv`, `joint_draws.csv` and (if available) `cv_predictions.csv` written by script 01. It prints each ggplot explicitly and saves it as a 300-dpi PNG in `plots/01_DO_dual_pathway_kinetics/`; no posterior refitting is required.

Plots 01–07 show observed versus fitted DO, δ¹⁸O–Cl endmember mixing, well LSR fractions, LSR–aquifer rate identifiability, per-well DO-loss partition, source-specific DO decay against MRT, and fitted DO residuals against MRT. Plot 08 adds held-out DO prediction if grouped CV outputs exist.

**Turbo palette and semantics:** river recharge is Turbo blue (`turbo(256)[32]`), rainfall/LSR recharge is Turbo red (`turbo(256)[230]`), and the continuous mixing fraction runs between those blue and red colours. Aquifer loss from both source waters is green (`turbo(256)[128]`), not river blue. Grey is used for reference lines only. These colours apply consistently across all plots.

**Uncertainty conventions:** `global_posteriors.csv` provides central **95%** posterior credible intervals. The saved per-well `fit$summary()` provides `q5` and `q95` (**90%** posterior credible intervals); these are correctly labelled as 90% in the plots. Fitted `mean_DO` credible intervals are **not predictive intervals**. Only the optional CV plot uses posterior predictive 95% intervals. The source-pathway decay curves use posterior median parameters and are illustrative, not posterior prediction bands.

## Diagnostics and outputs

The primary R runner saves:

- `outputs/01_DO_dual_pathway_kinetics/fit_order0.rds` and `fit_order1.rds` before downstream reporting.
- `fit_status.csv` and `sampler_chains.csv`: convergence, R-hat, bulk ESS, divergent trajectories, maximum treedepth and E-BFMI.
- `global_posteriors.csv`: posterior summaries and central 95% credible intervals for the mixing endmembers, LSR entry DO, LSR deficit and aquifer rates.
- `well_posteriors.csv`: mixture fractions, nominal losses, predicted DO, and inferred MRT for each well.
- `identifiability.csv`, `joint_draws.csv`, and `plots/01_DO_dual_pathway_kinetics/identifiability.png`: dependence between LSR deficit and aquifer rate.

On setting `run_CV=TRUE`, it additionally produces **matched five-fold grouped held-out DO predictions**, CV sampler checks and scores. A matched group bootstrap of first-minus-zero RMSE and MAE is produced only if all ten held-out model fits converge.

**Model acceptance:** An acceptable R-hat and ESS alongside no divergences, no treedepth saturation and E-BFMI≥0.3 are necessary; they do not establish correctness or chemical identifiability. The choice of aquifer reaction order requires held-out predictive evidence, not in-sample parameter agreement. Sensitivity to the assumed endmembers, recharge temperature, river DO entry boundary, MRT and LSR survival prior remains essential before publication.

## Provenance and limitations

The 23a sodium-ablation experiment found the Na-inclusive variants failed convergence diagnostics (worst R-hat about 1.30–1.41, E-BFMI about 0.1), while both noNa models passed (worst R-hat about 1.01, minimum bulk ESS above 600, minimum E-BFMI above 0.46, no divergences). That motivates excluding Na from the primary formulation. Earlier excess-air tests found limited useful information, so it is excluded as well. Neither finding proves those analytes intrinsically unsuitable as hydrological tracers; the evidence concerns these data and model specifications.

Two pathways share one modelled aquifer rate and an effective MRT. A low aquifer DO concentration can be explained through combinations of LSR water-table deficit, source fractions and aquifer reaction rates. Their partition is **conditional**, not separately measured, and requires sensitivity analysis. Arithmetic mean losses across wells do not represent a flow-weighted groundwater oxygen budget.

## High-DO observation audit (script 03)

The unaltered raw and prepared records give **2801: DO 10.2 mg/L, measured T 10.2 °C, MRT 0** and **2130: DO 8.8 mg/L, measured T 16.5 °C, MRT 8 years**. These are the only two measured DO observations above the fixed river boundary of 7.5 mg/L. Using the model's freshwater saturation polynomial, both are around 91% oxygen-saturated at their measured temperatures, so a high DO value alone is not evidence of a faulty measurement. Groundwater temperature may not equal recharge-equilibration temperature, and the calculation uses the assumed pressure ratio.

Script `03_high_DO_audit.R` analyses both observations and uses *existing* CmdStanR saved fit objects for site 2130, with no sampling or model modification. It compares 2130's measured DO with the most oxygen-rich theoretically possible source at its reported MRT under each fitted aquifer rate (even allowing zero LSR pre-water-table DO loss). At posterior median rates 0.236 mg/L/y (zero order) and 0.0625 y⁻¹ (first order), the optimistic maxima are approximately **8.15 mg/L and 6.09 mg/L**, respectively, both below measured 8.8 mg/L. Explaining 8.8 mg/L at those rate medians would require aquifer residence no greater than roughly **5.2 years** or **2.1 years**, respectively, even if the recharge entered fully saturated and was 100% from the highest-DO endmember. The saved-draw check tests whether posterior uncertainty in site 2130's true MRT and rate changes that assessment.

The two conservative tracer prior means imply very different simple source fractions at site 2130: approximately 0.68 LSR using δ¹⁸O versus 0.20 using chloride. These are illustrative *prior-endmember calculations*, not the fitted posterior mixing fractions; the mismatch motivates checking spatial/temporal endmember variation and the assumption of a single conservative two-source line. The masked DO for proxy 2801 must not be used to score predictive model fit. Site 2130 *was included* in model calibration, so in-sample posterior predictive tail probabilities are diagnostics, not independent validation.

## Experiment 04: river-entry DO uncertainty

We separate **river-entry DO boundary uncertainty** from the original 2801 pseudo-boundary practice. The measured DO of river-proximal groundwater site **2801 is 10.2 mg/L**, and `00_setup.R` no longer replaces it with 7.5. `01_DO_dual_pathway_kinetics.R` still masks 2801 DO in the original fixed-boundary model; the previously fitted 01 results remain comparable.

The new scripts `04_river_DO_prior_sensitivity.R` and `stan/04_river_DO_prior_sensitivity.stan` compare, on **exactly the same input observations and mask**, a fixed river-entry boundary of 7.5 mg/L with Normal priors centred at 7.5 mg/L, SD **1.0** and **1.5 mg/L**, truncated to [0, 14.65] mg/L. The cap approximates cold freshwater saturation rather than the distinct 15°C rainfall boundary (10.03 mg/L). The fixed model has NO unused boundary parameter: `river_DO_free` is a zero-length Stan vector unless the boundary is inferred. All two-endmember mixing, LSR survival, aquifer kinetics, and non-centred MRT assumptions are otherwise unchanged.

The default experiment includes the original 2801 DO observation as one ordinary groundwater response at a fixed reported MRT of zero, **not** as a direct measurement or fixed prior for river entry. In all scenario fits, the DO mask is identical. To assess dependence on that assumption, `cfg04$include_proxy_DO=FALSE` repeats the matched six-fit suite with 2801 DO masked and writes to a different results folder.

Output files include `scenario_design.csv`, `fit_status.csv`, `sampler_chains.csv`, `global_posteriors.csv`, `comparison_summary.csv`, `posterior_correlations.csv`, and `high_DO_site_predictions.csv`, plus independently saved full-fit RDS objects. The existing script 02 automatically generates two additional Turbo-colour figures if the experiment results are present: **river (blue) versus LSR (red) DO entry concentrations**, and **LSR deficit (red) versus aquifer reaction rates (green)**. No plot execution launches MCMC.

A widening posterior or strong correlation of river-entry DO with aquifer rate or LSR deficit would indicate insufficient separability; a high in-sample predictive tail probability for 2801 alone is not scientific validation. The direct effect of **including 2801** must not be conflated with the effect of **uncertain river DO**. Do not adopt the uncertain river boundary as primary until full diagnostics, predictive validation, and boundary/age sensitivity have been reviewed.


## Five-fold grouped held-out DO validation (script 05)

A rigorous validation suite compares fixed river DO at 7.5 mg/L to a truncated Normal(7.5, 1.5²) global river-entry DO prior under the same zero- and first-order kinetic laws. The model and chemical observation processes are reused from `stan/04_river_DO_prior_sensitivity.stan`. Four variants and five **identical, group-held-out DO folds** require 20 distinct Stan fits. Source isotope and chloride and reported MRT of held-out wells are included in the mixing/age model, but DO likelihood indicators **and DO values** for the held-out wells are zeroed. This is held-out DO prediction conditional on tracer and age information, not a fully external validation of source endmembers.

The primary design masks site 2801's measured 10.2 mg/L DO from **all** likelihoods and never scores it, while retaining its δ¹⁸O/Cl chemistry. All 40 remaining DO samples are scored once per model, including 2130. The alternative "2801_in_training" suite includes the actual 2801 DO in each training fit and scores exactly the **same** 40 wells. No 2801 pseudo-observation at 7.5 mg/L is created.

The fitted conditional DO means are averaged over the retained posterior draws as point predictions. Predictive 95% limits come directly from Stan-generated `replicated_DO[i]` draws, using the same heavy-tailed Student-t observation model as the training likelihood. Scores are RMSE, MAE, bias and 95% predictive coverage; MCMC checks include all active parameters and aggregate LSR loss, with R-hat/ESS, E-BFMI, divergence and treedepth checks. Model rankings and group-resampled 95% comparison intervals are withheld until all 20 fits have completed and passed diagnostics.

Paired comparisons: (i) C-minus-A at each kinetic order asks whether the broad river prior improves prediction; (ii) first-minus-zero within each river DO assumption asks which kinetic law predicts better. The paired bootstrap resamples pre-declared named sample groups, not individual rows, preserving group dependence. Outputs are separate from experiment 04 and the baseline; fits are individually cached with a saved input/design fingerprint to avoid accidentally combining different experiments. Script 02 can graph the held-out results using Turbo with blue river, red LSR and well point colours based on source fraction.

The CV does not by itself identify the physical location of oxygen loss; structural sensitivity and source-water monitoring remain important.