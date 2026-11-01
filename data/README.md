# Local private data

Copy the unchanged original prepared research file to:

`data/Uwe_HP_data_prepared.csv`

The **primary two-endmember DO model** (`scripts/01_DO_dual_pathway_kinetics.R`) requires `Site_ID`, `Depth`, `DO`, `MRT`, `O18` and `Cl`. The CSV may contain other measured analytes, but **sodium and excess air are not used** in the new Stan model.

Site **16758** is excluded as a far-south outlier. River-proximal site **2801** has a fixed, model-only DO value of 7.5 mg/L and its DO is masked from the likelihood; its original measured DO is not altered. Input data remain local and must never be committed to this repository.
