# 00_setup.R — primary DO two-endmember model setup (NO MCMC)
#
# Responsibilities: configure model / sampling options, validate CmdStan and
# the private input, prepare the two-pot O18 + Cl mixing observations and
# DO boundary conditions, and expose cfg, d, stan_data and other inputs.
# It does NOT compile the model, run Stan sampling, write analysis outputs,
# or source script 01. Source it alone to inspect the input and assumptions.
#
# Run from repository root:
#   source("scripts/00_setup.R")
# Execution of MCMC is ONLY in scripts/01_DO_dual_pathway_kinetics.R.

suppressPackageStartupMessages(library(tidyverse))

cfg <- list(
  input = "data/Uwe_HP_data_prepared.csv",
  stan = "stan/01_DO_dual_pathway_kinetics.stan",
  excluded_site = "16758",
  river_proxy_site = "2801",
  river_DO_boundary = 7.5,
  recharge_temperature_C = 15,
  atmospheric_pressure_ratio = 1,
  C_ref = 7.5,
  survival_prior_a = 1,
  survival_prior_b = 1,
  endmember_prior_scale = 1,
  log_mrt_SD = 0.25,
  seed = 20261031L,
  chains = 4L,
  parallel_chains = 4L,
  warmup = 1400L,
  sampling = 1400L,
  adapt_delta = 0.98,
  max_treedepth = 13L,
  max_rhat = 1.05,
  min_bulk_ESS = 400,
  run_CV = FALSE,
  folds = 5L,
  run_only_fold = NULL
)

if (!file.exists(cfg$stan))
  stop("Set working directory to repo root: missing ",cfg$stan)
if (!file.exists(cfg$input))
  stop("Private input missing: ",cfg$input)
if (!requireNamespace("cmdstanr", quietly=TRUE) ||
    !requireNamespace("posterior", quietly=TRUE))
  stop("Install cmdstanr, posterior and CmdStan.")
tryCatch(cmdstanr::cmdstan_path(),
         error=function(e) stop("CmdStan unavailable: ",conditionMessage(e)))

outdir <- "outputs/01_DO_dual_pathway_kinetics"
plotdir <- "plots/01_DO_dual_pathway_kinetics"
output <- function(name) file.path(outdir, name)

do_sat_mgL <- function(T, pressure_ratio=1) {
  # Freshwater approximation at 0-30 C; pressure correction approximate.
  (14.652 - 0.41022*T + 0.007991*T^2 -
     0.000077774*T^3)*pressure_ratio
}
stopifnot(
  cfg$recharge_temperature_C>=0,
  cfg$recharge_temperature_C<=30,
  cfg$atmospheric_pressure_ratio>0,
  cfg$river_DO_boundary>0,
  cfg$C_ref>0,
  cfg$log_mrt_SD>0,
  cfg$survival_prior_a>0,
  cfg$survival_prior_b>0
)
Csat <- do_sat_mgL(
  cfg$recharge_temperature_C,
  cfg$atmospheric_pressure_ratio
)

raw <- readr::read_csv(
  cfg$input,
  show_col_types=FALSE,
  col_types=readr::cols(
    Site_ID=readr::col_character(),
    .default=readr::col_guess()
  )
)
needed <- c("Site_ID","Depth","DO","MRT","O18","Cl")
missing <- setdiff(needed, names(raw))
if (length(missing))
  stop("Missing input columns: ",paste(missing,collapse=", "))
stopifnot(
  nrow(raw)>=40L,
  sum(raw$Site_ID==cfg$river_proxy_site,na.rm=TRUE)==1L,
  any(raw$Site_ID==cfg$excluded_site,na.rm=TRUE)
)

d <- raw %>%
  filter(
    Site_ID!=cfg$excluded_site,
    is.finite(DO), DO>=0,
    is.finite(MRT), MRT>=0
  ) %>%
  mutate(
    sample_index=row_number(),
    cv_group=case_when(
      str_detect(Site_ID,"Eastbourne")~"Eastbourne",
      str_detect(Site_ID,"Gilbertson")~"Gilbertson",
      str_detect(Site_ID,"Flaxmere")~"Flaxmere",
      TRUE~Site_ID
    )
  )
stopifnot(nrow(d)>=35L)
proxy_index <- which(d$Site_ID==cfg$river_proxy_site)
if (length(proxy_index)!=1L)
  stop("Missing or duplicated river proxy after filtering.")
if (any(d$MRT==0 & d$Site_ID!=cfg$river_proxy_site))
  stop("Unexpected zero-MRT sample: review age assumptions before fitting.")
# Preserve the ORIGINAL observed DO at every site, including proxy 2801.
# Site 2801 has measured DO=10.2 mg/L in the current source dataset.
# Never overwrite it with the river-entry boundary (7.5 mg/L).
# Baseline model 01 masks 2801 in its DO likelihood; the sensitivity
# model 04 can explicitly include its measured DO as a groundwater sample.
d$DO_original <- d$DO

positive_mrt <- d$MRT>0
mrt_index <- as.integer(ifelse(
  positive_mrt,cumsum(positive_mrt),0L
))
stopifnot(
  identical(sort(unique(mrt_index[mrt_index>0L])),
            seq_len(sum(positive_mrt)))
)
as_obs <- function(x) ifelse(is.finite(x),as.numeric(x),0)

stan_data <- list(
  N=nrow(d),
  has_O18=as.integer(is.finite(d$O18)),
  has_Cl=as.integer(is.finite(d$Cl)),
  has_DO=as.integer(seq_len(nrow(d))!=proxy_index),
  O18=as_obs(d$O18),
  Cl=as_obs(d$Cl),
  DO_obs=as.numeric(d$DO_original),
  Cl_sd=pmax(1,0.03*as_obs(d$Cl)),
  mrt_reported=as.numeric(d$MRT),
  N_mrt_positive=as.integer(sum(positive_mrt)),
  mrt_index=mrt_index,
  aq_order=0L,
  Csat_rain=Csat,
  C_ref=cfg$C_ref,
  endmember_prior_scale=cfg$endmember_prior_scale,
  log_mrt_SD=cfg$log_mrt_SD,
  river_DO_boundary=cfg$river_DO_boundary,
  beta_survival_a=cfg$survival_prior_a,
  beta_survival_b=cfg$survival_prior_b
)
eligible <- which(stan_data$has_DO==1L)
if (!length(eligible) ||
    !sum(stan_data$has_O18) ||
    !sum(stan_data$has_Cl))
  stop("Insufficient DO, delta18O or chloride data.")

message("01 TWO-POT DO MODEL: ",nrow(d)," wells; ",
        length(eligible)," fitted DO responses; ",
        sum(stan_data$has_O18)," isotope and ",
        sum(stan_data$has_Cl)," chloride observations.")
message("Fixed river DO ",cfg$river_DO_boundary,
        " mg/L; rainfall saturation DO ",signif(Csat,4),
        " mg/L. 2801 original DO ",d$DO_original[proxy_index],
        " mg/L preserved unchanged; masked in baseline 01. ",
        "Sodium and excess air NOT modelled.")

