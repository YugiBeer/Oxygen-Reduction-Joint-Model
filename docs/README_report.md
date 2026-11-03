# Heretaunga DO model report

The report source is [`docs/heretaunga_do_model_report.tex`](heretaunga_do_model_report.tex).
Its five DO–MRT figures are produced from **saved posterior draws** by
`scripts/05_DO_kinetic_trajectories and summary.R`. The report never embeds
screenshot copies, and the plotting script does not run Stan or refit any model.

## Rebuild in Positron

Pull the latest GitHub changes, open the repository root as your R working directory,
and ensure the existing private `data/Uwe_HP_data_prepared.csv` and saved preferred-model
RDS/diagnostic files are still available locally.

In the R console:

```r
source("scripts/05_DO_kinetic_trajectories and summary.R")
```

The script writes these **five numbered PNG files** into `plots/final_39well/`:

| Script plot | Exact output filename | What it displays |
|---|---|---|
| **1** | `01_DO_mixing_adjusted_decay_2801_masked_no2130.png` | Observed DO normalised by inferred initial mixed DO, against the common fitted first-order survival curve and its 95% credible band |
| **2** | `02_DO_observed_vs_fitted_mixed_2801_masked_no2130.png` | Measured versus individually fitted mixed DO, with posterior predictive intervals |
| **3** | `03_DO_river_reference_raw_samples_2801_masked_no2130.png` | Raw measured DO versus a hypothetical 100%-river first-order reference, which is not a fit to mixed groundwater |
| **4** | `04_DO_river_contribution_2801_masked_no2130.png` | Inferred river-attributable contribution to each mixed sample versus the unweighted river-only reference |
| **5** | `05_DO_fitted_components_2801_masked_no2130.png` | River and LSR contributions stacked to the fitted total DO, alongside measured DO |

The separate four-parameter model summary is written as
`06_DO_primary_posterior_summary_2801_masked_no2130.md`.

**Interpretation:** Plot 1's mixing-adjusted observed points use initial concentrations
inferred from the *same* joint model as the rate `k`. They are explanatory, not
independent observations of DO survival or a separate estimate of `k`.
Their point intervals reflect uncertainty in the inferred normalising concentration.
Plot 2 uses 95% posterior **predictive** intervals. Plot 5 uses posterior **means**
for additive source components and 95% **credible** intervals for the fitted
mean; its uncertainty whiskers are not posterior predictive intervals.
The shared-MRT assumption remains untested by these plots.

## Compile report with all five figures

From the repository root, in a terminal:

```bash
bash docs/build_heretaunga_report.sh
```

The build script checks **all five** numbered PNGs, then runs `pdflatex` twice and
writes `docs/heretaunga_do_model_report.pdf`. If a PNG is missing or empty,
the script stops and asks you to regenerate the figures.

The `outputs/` and `plots/` folders and private input data are intentionally
gitignored. Pulling the repository updates **R source, LaTeX, build scripts
and documentation**, not the locally cached posterior draws or generated PNGs.

A separately maintained GitHub-compiled PDF is just a **preview** with
explicit missing-output notices, because GitHub Actions cannot access
the gitignored local figures. The full report is the version built locally
after generating the PNGs.

No changes to `scripts/00_setup.R`, Stan, existing fitted RDS files, the
calibration likelihood or grouped CV are needed for this figure update.
