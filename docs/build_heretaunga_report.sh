#!/usr/bin/env bash
# Compile with the CURRENT plot outputs, never with uploaded screenshots.
set -euo pipefail
here="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
root="$(cd -- "$here/.." && pwd)"
figures=(
  "13_DO_observed_vs_fitted_mixed_2801_masked_no2130.png"
  "14_DO_river_reference_raw_samples_2801_masked_no2130.png"
  "15_DO_river_contribution_2801_masked_no2130.png"
)
for name in "${figures[@]}"; do
  file="$root/plots/final_39well/$name"
  if [[ ! -s "$file" ]]; then
    echo "Missing generated figure: $file" >&2
    echo 'Run source("scripts/05_DO_kinetic_trajectories and summary.R") in Positron from repo root.' >&2
    exit 2
  fi
done
command -v pdflatex >/dev/null || { echo "pdflatex is required for this report" >&2; exit 3; }
cd "$here"
pdflatex -interaction=nonstopmode -halt-on-error heretaunga_do_model_report.tex
pdflatex -interaction=nonstopmode -halt-on-error heretaunga_do_model_report.tex
printf 'Created %s/heretaunga_do_model_report.pdf with all three script-generated plots.\n' "$here"
