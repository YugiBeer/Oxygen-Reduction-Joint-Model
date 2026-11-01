# Local, private data

Place the source research file `Uwe_HP_data_prepared.csv` in this folder on your computer before running `scripts/00_check_setup.R` and `scripts/23_final_DO_two_pathway_kinetics.R`.

The file must contain `Site_ID`, `Depth`, `DO`, `MRT`, `O18`, `Cl`, `Na`, `xs_Air`, and `xs_Air_uncertainty`. It is **not** committed to this public repository. Keep the original measurements unchanged; script 23 masks proxy site 2801 DO and excludes site 16758 in memory.
