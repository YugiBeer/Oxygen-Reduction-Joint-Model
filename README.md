# Oxygen Reduction — Joint Model

A standalone two-script workflow for the 30 October 2026 Heretaunga groundwater analysis. It is independent of the earlier **fixed-rainfall-DO** rate calculations.

## Files

- `scripts/01_three_tracer_mixing.R` — source fractions from δ¹⁸O, chloride and excess air; three-tracer agreement screen and rainfall-endmember sensitivity.
- `scripts/02_joint_oxygen_partition.R` — simultaneous fitting of source-associated oxygen deficit and aquifer reduction, comparing first- and zero-order models using leave-one-well-out prediction and a conditional bootstrap.
- `data/groundwater.csv` — **not committed**. The source research repository is private and this repository is public; use the groundwater.csv provided in the separately prepared workflow ZIP.

## Run in Positron

Install `tidyverse` if needed (`install.packages('tidyverse')`). Open this repository as your working directory, place the prepared dataset at `data/groundwater.csv`, then run:

```r
source("scripts/01_three_tracer_mixing.R")
source("scripts/02_joint_oxygen_partition.R")
```

Results are written to `outputs/` and figures to `plots/`. These directories are Git-ignored.

## Assumptions

- Mixing uses site 2801 as the river **tracer endpoint**, and 1032, 4294 and the 108 m HB Prison sample as the land-surface **tracer endpoint**.
- For each well, the median of the clipped rainfall fractions from δ¹⁸O, chloride and excess air is used. Reliability requires all three tracers and a range of **unclipped** fractions ≤ 0.35. No minimum river recharge fraction.
- 29 of 42 records pass the production screen (five rainfall-dominated).
- River DO = 7.5 mg/L; **assumed** pre-infiltration rainfall DO reference = 9.5 mg/L.
- For the first-order model:
  
  DO_pred = [fR * 7.5 + fP * (9.5 - L)] * exp(-k * MRT)

  L and k are fitted together. For zero-order reduction, the post-mixing concentration is reduced by k * MRT, bounded below by zero.

Only **effective** rainfall-source DO (9.5 - L) is constrained by the data. An inferred source-associated deficit cannot, by itself, be assigned uniquely to infiltration or the unsaturated zone. The first-order rate is an apparent residence-time-dependent rate, not independent proof of saturated-zone reaction.

Excess air is inferred from dissolved gases; its treatment as a conservative mixing tracer is an **assumption needing investigation**. This analysis does not model excess N₂ or denitrification kinetics.

## Verification

The prepared standalone model and selection were numerically cross-checked against the earlier analysis. R itself was not executed in the preparation environment; please run both scripts in Positron to validate the outputs.