# Pre-flight checks for the PRIMARY 01_DO_dual_pathway_kinetics workflow.
# Run from the Oxygen-Reduction-Joint-Model repository root:
# source("scripts/00_check_setup.R")
# No model fits, modifications to original data or output files.

required_paths <- c(
  "scripts/01_DO_dual_pathway_kinetics.R",
  "stan/01_DO_dual_pathway_kinetics.stan",
  "data/Uwe_HP_data_prepared.csv"
)
missing_paths <- required_paths[!file.exists(required_paths)]
if (length(missing_paths)) {
  stop("Missing project files relative to ",getwd(),":\n  ",
       paste(missing_paths,collapse="\n  "),
       "\nSet the working directory to the repository root.")
}
required_packages <- c("tidyverse","posterior","cmdstanr")
not_installed <- required_packages[!vapply(
  required_packages,requireNamespace,logical(1),quietly=TRUE
)]
if (length(not_installed))
  stop("Missing R packages: ",paste(not_installed,collapse=", "))
cmdstan_path <- tryCatch(
  cmdstanr::cmdstan_path(),
  error=function(e) stop("CmdStan not configured: ",conditionMessage(e))
)

raw <- readr::read_csv(
  "data/Uwe_HP_data_prepared.csv",
  show_col_types=FALSE,
  col_types=readr::cols(
    Site_ID=readr::col_character(),
    .default=readr::col_guess()
  )
)
required_columns <- c("Site_ID","Depth","DO","MRT","O18","Cl")
missing_columns <- setdiff(required_columns,names(raw))
if(length(missing_columns))
  stop("Missing private input columns: ",
       paste(missing_columns,collapse=", "))
if (sum(raw$Site_ID=="2801",na.rm=TRUE)!=1L)
  stop("Expected exactly one site 2801 proxy record.")
if (!any(raw$Site_ID=="16758",na.rm=TRUE))
  stop("Site 16758 missing; verify input version.")
keep <- !is.na(raw$Site_ID) & raw$Site_ID!="16758" &
  is.finite(raw$DO) & raw$DO>=0 &
  is.finite(raw$MRT) & raw$MRT>=0
if (sum(keep)<35L)
  stop("Too few usable records: ",sum(keep))
if (!all(raw$MRT[keep & raw$Site_ID!="2801"]>0))
  stop("Zero reported MRT outside proxy site 2801.")
message(
  "PREFLIGHT OK | root: ",normalizePath(getwd()),
  " | CmdStan: ",cmdstan_path,
  " | original rows: ",nrow(raw),
  " | usable two-pot samples: ",sum(keep),
  " | fitted DO observations: ",sum(keep)-1L,
  " | delta18O observations: ",sum(is.finite(raw$O18[keep])),
  " | Cl observations: ",sum(is.finite(raw$Cl[keep]))
)
message('Run: source("scripts/01_DO_dual_pathway_kinetics.R")')
