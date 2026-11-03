# Tracer-only river/LSR groundwater mixing fractions

This is a **separate Bayesian model** that estimates the **water mixing fractions** of river recharge and land-surface recharge (LSR) from **δ¹⁸O and chloride only**. It does **not** use dissolved oxygen, mean residence time (MRT), excess N₂, nitrate, or the preassigned `source` category when inferring fractions.

## Run

From the repository root in Positron, with the original `data/Uwe_HP_data_prepared.csv` available locally:

```r
source("scripts/07_export_tracer_only_fractions.R")
```

The **first run fits a new four-chain Stan model** (it cannot extract tracer-only posterior fractions from the existing DO-conditioned fit). Subsequent runs reuse the saved tracer-only fit if the input and model checksums and sampling configuration still match. This is completely separate from the existing DO fit/CV workflow.

The script creates:

- `data/Uwe_HP_data_prepared_with_tracer_fractions.csv` — shareable copy of the complete original data with eight appended fields; **original file unchanged**
- `outputs/08_tracer_only_mixing/tracer_only_fit.rds` — full posterior, useful for downstream uncertainty propagation
- `outputs/08_tracer_only_mixing/tracer_only_diagnostics.csv` — convergence and data-quality checks
- `outputs/08_tracer_only_mixing/tracer_only_method.txt` — reproducibility and interpretation notes

Both `data/*.csv` and `outputs/` are gitignored, so **no raw sample data or derived data will be pushed to GitHub**.

## Added CSV fields

| Column | Meaning |
|---|---|
| `LSR_fraction_median` | Posterior median proportion of **water** originating from LSR, 0–1 |
| `LSR_fraction_lower95`, `LSR_fraction_upper95` | Equal-tailed 95% Bayesian credible interval for LSR water fraction |
| `river_fraction_median` | Posterior median proportion of water originating from river recharge; `1 - LSR_fraction_median` |
| `river_fraction_lower95`, `river_fraction_upper95` | Equal-tailed river-fraction interval; endpoints are the complements of the opposite LSR bounds |
| `mixing_tracers_used` | `O18 + Cl`, `O18 only`, `Cl only`, or `none` |
| `mixing_status` | `two tracers`, `one tracer - caution`, `no available tracer` or `excluded (16758)` |

The input data contain 43 rows (with a repeated `Site_ID` at one site). The script operates by **input row**, not by joining on `Site_ID`, and preserves the original row order and every existing column. Site **16758** is excluded from fitting but retained in the output with missing fractions. Sites **2801 and 2130** remain included because their conservative tracer observations are valid. The script also allows a sample with unknown MRT or missing DO, since neither variable is part of tracer-only inference. One-tracer samples receive fractions but are flagged; these are generally less constrained than two-tracer estimates.

## Scientific method and limitations

The standalone model at `stan/05_tracer_only_mixing.stan` duplicates the **δ¹⁸O and chloride submodel** and default endmember prior specifications of the main DO model, `stan/04_river_DO_prior_sensitivity.stan`: river/LSR endmember isotope and chloride concentrations are uncertain; each row has its own `f_LSR`, with a Beta(1,1) prior. The model also estimates tracer residual variation, including the LSR isotope-dispersion term. Source fractions are volumetric water fractions, **not** fractional contributions to measured DO or excess N₂.

The posterior medians and credible intervals summarise marginal uncertainties. For MRT, excess N₂ or denitrification calculations, mixing source water concentrations may be reasonable, but apparent tracer ages and noble-gas excesses require appropriate tracer-specific transformations, equilibration assumptions and endmember uncertainties. **Do not treat a measured apparent MRT as a simple linear age mixture without considering the dating method.** For rigorous downstream uncertainty propagation, use the joint posterior draws from the saved tracer-only fit rather than substituting independent marginal medians.

The export is withheld if convergence diagnostics fail (R-hat, effective sample size, divergences, E-BFMI or maximum treedepth hits). The source input is never overwritten.
