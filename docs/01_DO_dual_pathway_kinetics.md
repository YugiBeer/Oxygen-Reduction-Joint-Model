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

Private input: `data/Uwe_HP_data_prepared.csv`. Required columns: `Site_ID`, `Depth`, `DO`, `MRT`, `O18` and `Cl`. Site 16758 is excluded. Site 2801's original DO measurement is not changed in the CSV; in-memory value 7.5 mg/L is not a fitted or scored DO observation. With the currently prepared file, 41 wells remain after exclusion and basic DO/MRT filtering, 40 of which contribute DO likelihood observations.

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
