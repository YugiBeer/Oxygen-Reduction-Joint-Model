# Final Heretaunga DO model — 39-well regional analysis

## Scope and exclusions

Data from `data/Uwe_HP_data_prepared.csv` are preprocessed by script 00. Site 16758 is excluded; site 2801 (river-proximal, measured DO 10.2 mg/L, assumed MRT 0) and site 2130 (shallow, water-race adjacent, measured DO 8.8 mg/L) are **masked from the DO likelihood**, while their oxygen-18 and chloride tracer information remains. Site 2130 is excluded on hydrogeological representativeness grounds; the older 40-well experiment is retained in local historical outputs and the Git history, not used as the regional primary.

## Model structure

Bayesian two-endmember mixing with per-well source fractions inferred from δ¹⁸O and chloride. Rainfall/LSR DO begins at freshwater saturation (assumed 15°C and atmospheric pressure ratio 1); the water-table LSR entry DO is inferred through a survival fraction. River-entry DO is estimated under a Normal(7.5, 1.5²) mg/L prior truncated [0,14.65] (**physical primary**); the fixed 7.5 mg/L boundary is a **sensitivity/predictive comparator**. Each source's DO is then exposed to the same aquifer consumption kinetics for its groundwater mean residence time before mixing. The two independently fitted aquifer kinetics are:
- Zero order: `max(0, C_entry - k0*t)`; k0 in mg/L/year.
- First order: `C_entry*exp(-k1*t)`; k1 in year^-1.

Positive measured MRT values use non-centred lognormal uncertainty (log-SD 0.25). Source mixture priors, a common LSR entry DO and a shared aquifer reaction rate are important structural assumptions. The model does not estimate an independent vadose-zone reaction rate.

## Model selection and principal estimates

Five-fold grouped CV (39 held-out DO observations per model, 20 matched fits) compares fixed and broad river DO × both aquifer kinetic orders, using grouped bootstrap paired differences. In the fixed river-DO scenario, first-order RMSE improvement over zero-order is median 0.132 mg/L (95% bootstrap improvement 0.0126–0.238 mg/L); the corresponding MAE interval overlaps zero. Broader river DO does not improve first-order predictive RMSE. First-order kinetics are preferred. Although fixed river DO has slightly better grouped CV RMSE, the final *physical primary* formulation uses the **uncertain river DO boundary** to avoid treating 7.5 mg/L as exactly known. The fixed-boundary first-order model is retained as the predictive comparator. This distinction should be explicit in reporting.

Earlier script-04 CSV summaries for the uncertain-river/first-order **primary**:
- `river_DO`: 7.66 mg/L (95% CrI 5.75–9.74)
- `aquifer_rate_actual`: 0.0632 yr^-1 (95% CrI 0.0453–0.0851)
- `LSR_DO_at_water_table`: 3.21 mg/L (95% CrI 1.25–5.55)
- `LSR_DO_deficit`: 6.79 mg/L (95% CrI 4.45–8.75)

The **latest script-07 baseline**, read from the saved primary CmdStanR fit object, instead reports `k1 = 0.0639` yr^-1 (95% CrI 0.0458–0.0857), LSR WT DO 3.26 (1.31–5.57) mg/L, and deficit 6.78 (4.46–8.72) mg/L. These are similar but **not identical**. Script 06 now reads and summarises that saved RDS directly and warns if the old CSV disagrees; reconcile provenance before manuscript submission. The fixed-river first-order **comparator** was k1=0.0625 yr^-1 (0.0484–0.0812), LSR WT DO 3.28 (1.42–5.47), and deficit 6.72 (4.53–8.58). Broad-river/first-order posterior river DO versus rate correlation was previously +0.544, indicating some compensation; the fixed-boundary comparison's LSR deficit versus rate correlation was -0.397.

## Numerical convergence

All four full-data fits passed R-hat (max 1.01), bulk-ESS (minimum 713 across monitored population/latent parameters), no divergences/treedepth failures and E-BFMI (minimum ~0.444). For the preferred model's four target parameters, bulk ESS was 3,628–5,070 and tail ESS 3,013–4,373. Longer MCMC chains are not presently warranted.

## Scientific limitations

The LSR deficit is inferred, **not** a direct measurement of unsaturated-zone oxygen reduction. The apparent spatial attribution of DO consumption is conditional on recharge endmembers, groundwater age assumptions, a common water-table LSR entry concentration, spatial/temporal river-boundary variability and the two-source structure. The inferred mean loss is not a spatially/volumetrically weighted aquifer oxygen budget. Validation is for the defined 39-well regional dataset; further independent measurements are preferable to additional free model parameters.

## Active workflow

Use `scripts/00_setup.R` (preflight), `scripts/04_river_DO_prior_sensitivity.R` (four full-data posterior fits), `scripts/05_river_DO_grouped_CV.R` (matched CV with cache), `scripts/02_DO_plotting.R` (plot 04/05 CSVs), `scripts/06_DO_posterior_distributions.R` (one primary first-order rate posterior figure from the saved RDS), and `scripts/07_DO_LSR_MRT_sensitivity.R` (separate assumption-sensitivity study with cached fits). Current RDS/CSV outputs remain local and gitignored. Older protocols and scripts are archived rather than deleted.

## Completed targeted LSR/MRT sensitivity study (script 07; user-reported outputs)

This study isolates three previously untested assumption categories, **one factor at a time**, using the selected first-order aquifer kinetics and broad river-DO prior (mean 7.5, SD 1.5 mg/L). It **reuses** the validated `04` posterior `fit_C_prior_sd1p5_order1.rds` as the baseline and writes new results separately under `outputs/07_DO_LSR_MRT_sensitivity/2801_masked_no2130/`.

- **Surface DO saturation** at assumed LSR recharge temperatures of 10, 15 (baseline) and 20°C. Changing this alters both the land-surface DO reference and the prior induced on water-table DO.
- **LSR survival prior** on the water-table fraction: Beta(1,1) baseline, Beta(2,2), Beta(2,5) and Beta(5,2). The alternatives are deliberately contrasting prior-sensitivity probes and are not independently justified hydrological priors.
- **MRT uncertainty** log-standard deviation 0.10, 0.25 baseline and 0.50. The model still assumes the same well-specific aquifer residence time for both recharge pathways.

The user completed all seven additional fits and reported that **all eight (baseline plus seven new fits) passed** diagnostics: max R-hat 1.01, zero divergences and treedepth hits, minimum E-BFMI 0.460, minimum bulk ESS 768, and minimum tail ESS 669. These have not been independently re-run in this environment.

| Assumption | Scenario | Median k1 (yr^-1) | Change versus baseline | Median LSR WT DO (mg/L) | Median LSR deficit (mg/L) |
|---|---|---:|---:|---:|---:|
| Baseline | 15°C; Beta(1,1); MRT log-SD 0.25 | 0.0639 | 0% | 3.26 | 6.78 |
| Saturation temperature | 10°C | 0.0638 | -0.120% | 3.25 | 8.02 |
| Saturation temperature | 20°C | 0.0639 | +0.0363% | 3.23 | 5.79 |
| LSR survival | Beta(2,2) | 0.0638 | -0.0318% | 3.42 | 6.61 |
| LSR survival | Beta(2,5) | 0.0632 | -0.997% | 3.01 | 7.03 |
| LSR survival | Beta(5,2) | 0.0660 | +3.38% | 4.28 | 5.75 |
| MRT log-SD | 0.10 | 0.0610 | -4.55% | 3.13 | 6.91 |
| MRT log-SD | 0.50 | 0.0757 | +18.5% | 3.90 | 6.13 |

Therefore **the estimated first-order aquifer rate is much less sensitive to the tested LSR DO boundary and survival-prior choices than to the breadth of assumed groundwater-age uncertainty**. The LSR DO deficit responds to the assumed saturation concentration, while the inferred WT DO stays near 3.2–3.3 mg/L across the 10–20°C comparisons. The alternative survival priors are deliberate prior-stress tests, not independently constrained recharge priors. These are one-at-a-time conditional sensitivities. They do **not** demonstrate the validity of the assumed common-MRT age for river and progressively added LSR, or determine the sign of bias from source-specific ages.
