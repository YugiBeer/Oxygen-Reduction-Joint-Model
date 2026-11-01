# Pre-flight checks for the ACTIVE script 23 workflow.
# Run from the Oxygen-Reduction-Joint-Model repository root:
# source("scripts/00_check_setup.R")
# Does not fit the model, modify source data or call another script.

required_paths <- c(
  "scripts/23_final_DO_two_pathway_kinetics.R",
  "stan/23_final_DO_two_pathway_kinetics.stan",
  "data/Uwe_HP_data_prepared.csv"
)
missing_paths <- required_paths[!file.exists(required_paths)]
if (length(missing_paths)) {
  stop("Missing files relative to getwd() = ", getwd(), ":\n  ",
       paste(missing_paths, collapse = "\n  "),
       "\nOpen the repository root in Positron. Keep private data under data/.")
}

packages <- c("tidyverse", "posterior", "cmdstanr")
missing_pkgs <- packages[!vapply(
  packages, requireNamespace, logical(1), quietly = TRUE
)]
if (length(missing_pkgs)) {
  stop("Missing R packages: ", paste(missing_pkgs, collapse = ", "),
       ". Install these packages before running script 23.")
}

cmdstan_path <- tryCatch(cmdstanr::cmdstan_path(), error = function(e) {
  stop("CmdStan is not configured: ", conditionMessage(e))
})

raw <- readr::read_csv(
  "data/Uwe_HP_data_prepared.csv",
  show_col_types = FALSE,
  col_types = readr::cols(
    Site_ID = readr::col_character(),
    .default = readr::col_guess()
  )
)
required_columns <- c(
  "Site_ID", "Depth", "DO", "MRT", "O18",
  "Cl", "Na", "xs_Air", "xs_Air_uncertainty"
)
missing_columns <- setdiff(required_columns, names(raw))
if (length(missing_columns)) {
  stop("Private CSV missing columns: ",
       paste(missing_columns, collapse = ", "))
}
if (sum(raw$Site_ID == "2801", na.rm = TRUE) != 1L) {
  stop("Expected exactly one 2801 proxy record in the original CSV.")
}
if (!any(raw$Site_ID == "16758", na.rm = TRUE)) {
  stop("Expected site 16758 in input; script 23 excludes it explicitly.")
}

keep <- !is.na(raw$Site_ID) & raw$Site_ID != "16758" &
  is.finite(raw$DO) & raw$DO >= 0 &
  is.finite(raw$MRT) & raw$MRT >= 0
if (sum(keep) < 35L) {
  stop("Too few usable rows after excluding 16758 and missing DO/MRT: ",
       sum(keep))
}

message(
  "PREFLIGHT OK | project: ", normalizePath(getwd()),
  " | CmdStan: ", cmdstan_path,
  " | original rows: ", nrow(raw),
  " | usable rows: ", sum(keep),
  " | independent DO responses after masking 2801: ", sum(keep) - 1L
)
message("Now run: source(\"scripts/23_final_DO_two_pathway_kinetics.R\")")
