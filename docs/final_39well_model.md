# Final Heretaunga DO model — 39-well regional analysis

## Scope and exclusions

Data from `data/Uwe_HP_data_prepared.csv` are preprocessed by script 00. Site 16758 is excluded; site 2801 (river-proximal, measured DO 10.2 mg/L, assumed MRT 0) and site 2130 (shallow, water-race adjacent, measured DO 8.8 mg/L) are **masked from the DO likelihood**, while their oxygen-18 and chloride tracer information remains. Site 2130 is excluded on hydrogeological representativeness grounds; the older 40-well experiment is retained in local historical outputs and the Git history, not used as the regional primary.

## Model structure

Bayesian two-endmember mixing with per-well source fractions inferred from δ¹⁸O and chloride. Rainfall/LSR DO begins at freshwater saturation (assumed 15°C and atmospheric pressure ratio 1); the water-table LSR entry DO is inferred through a survival fraction. River-entry DO is fixed to 7.5 mg/L (primary), or estimated under a Normal(7.5, 1.5²) mg/L prior truncated [0,14.65] (sensitivity). Each source's DO is then exposed to the same aquifer consumption kinetics for its groundwater mean residence time before mixing. The two independently fitted aquifer kinetics are:
- Zero order: `max(0, C_entry - k0*t)`; k0 in mg/L/year.
- First order: `C_entry*exp(-k1*t)`; k1 in year^-1.

Positive measured MRT values use non-centred lognormal uncertainty (log-SD 0.25). Source mixture priors, a common LSR entry DO and a shared aquifer reaction rate are important structural assumptions. The model does not estimate an independent vadose-zone reaction rate.

## Model selection and principal estimates

Five-fold grouped CV (39 held-out DO observations per model, 20 matched fits) compares fixed and broad river DO × both aquifer kinetic orders, using grouped bootstrap paired differences. In the fixed river-DO scenario, first-order RMSE improvement over zero-order is median 0.132 mg/L (95% bootstrap improvement 0.0126–0.238 mg/L); the corresponding MAE interval overlaps zero. Broader river DO does not improve first-order predictive RMSE. First order / fixed river DO is preferred for its CV performance and parsimony, not considered proven uniquely correct.

The preferred model posterior has:
- `aquifer_rate_actual`: median 0.0625 yr^-1 (95% CrI 0.0484–0.0812)
- `LSR_DO_at_water_table`: 3.28 mg/L (95% CrI 1.42–5.47)
- `LSR_DO_deficit`: 6.72 mg/L (95% CrI 4.53–8.58)
- `river_DO`: 7.5 mg/L fixed (not inferred)

In the broad river DO / first-order sensitivity: reaction rate ~0.0632 yr^-1, LSR deficit ~6.79 mg/L. The preferred model's posterior correlation of LSR DO deficit vs aquifer rate is approximately -0.397. Broad river DO / first-order shows correlation +0.544 between inferred river DO and aquifer reaction rate, indicating some parameter trade-off.

## Numerical convergence

All four full-data fits passed R-hat (max 1.01), bulk-ESS (minimum 713 across monitored population/latent parameters), no divergences/treedepth failures and E-BFMI (minimum ~0.444). For the preferred model's four target parameters, bulk ESS was 3,628–5,070 and tail ESS 3,013–4,373. Longer MCMC chains are not presently warranted.

## Scientific limitations

The LSR deficit is inferred, **not** a direct measurement of unsaturated-zone oxygen reduction. The apparent spatial attribution of DO consumption is conditional on recharge endmembers, groundwater age assumptions, a common water-table LSR entry concentration, spatial/temporal river-boundary variability and the two-source structure. The inferred mean loss is not a spatially/volumetrically weighted aquifer oxygen budget. Validation is for the defined 39-well regional dataset; further independent measurements are preferable to additional free model parameters.

## Active workflow

Use `scripts/00_setup.R` (preflight), `scripts/04_river_DO_prior_sensitivity.R` (four full-data posterior fits), `scripts/05_river_DO_grouped_CV.R` (matched CV with cache), `scripts/02_DO_plotting.R` (plot final 04/05 CSVs) and `scripts/06_DO_posterior_distributions.R` (plot saved posterior draws). Current RDS/CSV outputs remain local and gitignored. Older protocols and scripts are archived rather than deleted.
