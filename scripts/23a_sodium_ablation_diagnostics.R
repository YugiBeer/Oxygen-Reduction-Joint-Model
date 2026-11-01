# 23a_sodium_ablation_diagnostics.R
#
# Controlled full-fit experiment: compare sodium included versus sodium
# completely masked, independently for aquifer orders 0 and 1.
#
# Same Stan model, prior distributions, input wells, four chains, sampler
# settings, random seeds by order, oxygen boundaries and MRT treatment.
# No cross-validation or sensitivity fitting in this diagnostic experiment.
# DO NOT change script 23's primary model based solely on MCMC diagnostics.
#
# No-Na analysis uses O18 + chloride only (excess air off), NOT three tracers.
# Na-specific nuisance parameters remain in Stan under proper priors but are
# PRIOR-ONLY when Na is masked. They are separately flagged in the summaries.
#
# Output goes to outputs/sodium_ablation_23a_* (never final_DO_23_*).
# The original private data are never modified or committed.

suppressPackageStartupMessages(library(tidyverse))

cfg <- list(
  input = "data/Uwe_HP_data_prepared.csv",
  stan = "stan/23_final_DO_two_pathway_kinetics.stan",
  river_DO = 7.5,
  proxy_site = "2801",
  excluded_site = "16758",
  recharge_T = 15,
  pressure_ratio = 1,
  C_ref = 7.5,
  endmember_prior_scale = 1,
  log_mrt_SD = 0.25,
  include_excess_air = FALSE,
  seed = 20261031L,
  chains = 4L,
  parallel_chains = 4L,
  warmup = 1400L,
  sampling = 1400L,
  adapt_delta = 0.98,
  max_treedepth = 13L,
  max_rhat = 1.05,
  min_bulk_ESS = 400
)

if (!file.exists(cfg$stan)) {
  stop("Set the working directory to the repository root; missing ", cfg$stan)
}
if (!file.exists(cfg$input)) {
  stop("Missing PRIVATE input file: ", cfg$input)
}
if (!requireNamespace("cmdstanr", quietly = TRUE) ||
    !requireNamespace("posterior", quietly = TRUE)) {
  stop("Install cmdstanr, posterior and CmdStan before running.")
}
cmdstanr::cmdstan_path()
dir.create("outputs", showWarnings = FALSE, recursive = TRUE)

as_obs <- function(x) ifelse(is.finite(x), as.numeric(x), 0)
as_sigma <- function(x, default) {
  ifelse(is.finite(x) & x > 0, as.numeric(x), default)
}
do_sat_mgL <- function(T, pressure_ratio = 1) {
  (14.652 - 0.41022 * T + 0.007991 * T^2 -
     0.000077774 * T^3) * pressure_ratio
}

raw <- readr::read_csv(
  cfg$input, show_col_types = FALSE,
  col_types = readr::cols(
    Site_ID = readr::col_character(),
    .default = readr::col_guess()
  )
)
required <- c(
  "Site_ID", "Depth", "DO", "MRT", "O18", "Cl", "Na",
  "xs_Air", "xs_Air_uncertainty"
)
if (length(setdiff(required, names(raw)))) {
  stop("Missing columns: ", paste(setdiff(required, names(raw)),
                                  collapse = ", "))
}
stopifnot(sum(raw$Site_ID == cfg$proxy_site, na.rm = TRUE) == 1L)
stopifnot(any(raw$Site_ID == cfg$excluded_site, na.rm = TRUE))

d <- raw %>%
  filter(
    Site_ID != cfg$excluded_site,
    is.finite(DO), DO >= 0,
    is.finite(MRT), MRT >= 0
  ) %>%
  mutate(sample_index = row_number())
stopifnot(nrow(d) >= 35L)

proxy <- which(d$Site_ID == cfg$proxy_site)
stopifnot(length(proxy) == 1L)
d$DO[proxy] <- cfg$river_DO  # model-only; excluded from likelihood
positive_mrt <- d$MRT > 0
mrt_index <- as.integer(ifelse(
  positive_mrt, cumsum(positive_mrt), 0L
))
stopifnot(
  identical(
    sort(unique(mrt_index[mrt_index > 0L])),
    seq_len(sum(positive_mrt))
  )
)
Csat <- do_sat_mgL(cfg$recharge_T, cfg$pressure_ratio)

make_data <- function(use_na) {
  list(
    N = nrow(d),
    N_mrt_positive = as.integer(sum(positive_mrt)),
    mrt_index = mrt_index,
    has_O18 = as.integer(is.finite(d$O18)),
    has_Cl = as.integer(is.finite(d$Cl)),
    has_Na = if (use_na) as.integer(is.finite(d$Na))
             else rep.int(0L, nrow(d)),
    has_air = if (cfg$include_excess_air)
                as.integer(is.finite(d$xs_Air))
              else rep.int(0L, nrow(d)),
    has_DO = replace(rep.int(1L, nrow(d)), proxy, 0L),
    O18 = as_obs(d$O18),
    Cl = as_obs(d$Cl),
    Na = if (use_na) as_obs(d$Na) else rep(0, nrow(d)),
    air = if (cfg$include_excess_air)
            as_obs(d$xs_Air)
          else rep(0, nrow(d)),
    DO_obs = as.numeric(d$DO),
    Cl_sd = pmax(1, 0.03 * as_obs(d$Cl)),
    Na_sd = if (use_na) pmax(0.5, 0.03 * as_obs(d$Na))
            else rep(0.5, nrow(d)),
    air_sd = as_sigma(d$xs_Air_uncertainty, 0.15),
    mrt_reported = as.numeric(d$MRT),
    aq_order = 0L,
    Csat_rain = Csat,
    C_ref = cfg$C_ref,
    endmember_prior_scale = cfg$endmember_prior_scale,
    log_mrt_SD = cfg$log_mrt_SD,
    river_DO_boundary = cfg$river_DO,
    beta_survival_a = 1,
    beta_survival_b = 1
  )
}

mod <- cmdstanr::cmdstan_model(cfg$stan, quiet = TRUE)
prefix <- "outputs/sodium_ablation_23a_"
stored <- list(
  sampler = list(), parameters = list(),
  core = list(), chains = list(), comparison = list()
)
save_tables <- function() {
  for (kind in names(stored)) {
    if (length(stored[[kind]])) {
      readr::write_csv(
        dplyr::bind_rows(stored[[kind]]),
        paste0(prefix, kind, ".csv")
      )
    }
  }
}

core_vars <- c(
  "LSR_DO_deficit", "LSR_DO_at_water_table", "survival_LSR",
  "aquifer_rate_actual", "k_ref", "sigma_DO", "mean_f_LSR",
  "sigma_Cl", "sigma_Na", "rho_Na_Cl", "Cl_LSR", "Na_LSR"
)
chain_vars <- c(
  "LSR_DO_deficit", "aquifer_rate_actual", "mean_f_LSR",
  "sigma_Cl", "sigma_Na", "rho_Na_Cl", "Cl_LSR", "Na_LSR"
)
prior_only <- c("air_R", "air_LSR", "sigma_air")
prior_only_no_na <- c(
  "Na_R", "Na_LSR", "sigma_Na", "rho_Na_Cl"
)

for (use_na in c(TRUE, FALSE)) {
  version <- if (use_na) "withNa" else "noNa"
  x <- make_data(use_na)

  message("SODIUM ABLATION: ", version, "; N=", x$N,
          "; Na measurements used=", sum(x$has_Na),
          "; DO observations used=", sum(x$has_DO))

  for (order in 0:1) {
    x$aq_order <- as.integer(order)
    message("Sampling: ", version, " aquifer order ", order)

    fit <- mod$sample(
      data = x,
      seed = cfg$seed + 101L * order,
      chains = cfg$chains,
      parallel_chains = cfg$parallel_chains,
      iter_warmup = cfg$warmup,
      iter_sampling = cfg$sampling,
      adapt_delta = cfg$adapt_delta,
      max_treedepth = cfg$max_treedepth,
      refresh = 250
    )

    # Preserve MCMC result before any downstream table generation.
    saved_fit <- paste0(prefix, version, "_order", order, ".rds")
    fit$save_object(file = saved_fit)

    diagnostics <- fit$diagnostic_summary()
    sampler <- tibble(
      variant = version,
      aquifer_order = order,
      chain = seq_len(cfg$chains),
      num_divergent = as.integer(diagnostics$num_divergent),
      num_max_treedepth = as.integer(diagnostics$num_max_treedepth),
      ebfmi = as.numeric(diagnostics$ebfmi)
    )
    stored$sampler[[length(stored$sampler) + 1L]] <- sampler

    sm <- tibble::as_tibble(fit$summary()) %>%
      mutate(
        variant = version, aquifer_order = order,
        prior_only_in_fit = variable %in% c(
          prior_only, if (!use_na) prior_only_no_na
        ),
        .before = 1
      )
    stored$parameters[[length(stored$parameters) + 1L]] <- sm

    active <- sm %>%
      filter(!prior_only_in_fit, is.finite(rhat),
             is.finite(ess_bulk))
    if (!nrow(active)) stop("No active variables with finite diagnostics.")
    bad_bfmi <- sum(!is.finite(sampler$ebfmi) |
                    sampler$ebfmi < 0.3)
    worst_rhat <- max(active$rhat)
    min_ess <- min(active$ess_bulk)
    divergent <- sum(sampler$num_divergent)
    treedepth <- sum(sampler$num_max_treedepth)

    stored$comparison[[length(stored$comparison) + 1L]] <- tibble(
      variant = version,
      aquifer_order = order,
      n_wells = x$N,
      n_na_observations = sum(x$has_Na),
      max_rhat_active = worst_rhat,
      min_bulk_ess_active = min_ess,
      n_divergences = divergent,
      n_max_treedepth = treedepth,
      min_ebfmi = min(sampler$ebfmi),
      n_chains_ebfmi_below_0p3 = bad_bfmi,
      passes_current_sampler_gates =
        worst_rhat <= cfg$max_rhat &&
        min_ess >= cfg$min_bulk_ESS &&
        divergent == 0 && treedepth == 0 && bad_bfmi == 0
    )

    core <- sm %>%
      filter(variable %in% core_vars) %>%
      select(
        variant, aquifer_order, variable, mean, median, sd,
        rhat, ess_bulk, ess_tail, prior_only_in_fit
      ) %>%
      mutate(
        kinetic_rate_units = if_else(
          variable == "aquifer_rate_actual",
          if (order == 0) "mg/L/year" else "year^-1",
          NA_character_
        )
      )
    stored$core[[length(stored$core) + 1L]] <- core

    draws <- tibble::as_tibble(posterior::as_draws_df(
      fit$draws(variables = chain_vars)
    ))
    chain_medians <- draws %>%
      group_by(.chain) %>%
      summarise(
        across(all_of(chain_vars), median),
        .groups = "drop"
      ) %>%
      mutate(variant = version, aquifer_order = order, .before = 1)
    stored$chains[[length(stored$chains) + 1L]] <- chain_medians

    save_tables()
    print(tail(dplyr::bind_rows(stored$comparison), 1))
  }
}

cat("\nFULL MATCHED SODIUM ABLATION COMPARISON\n")
print(dplyr::bind_rows(stored$comparison), n = Inf)
cat("\nOutput tables: ", prefix, "{comparison,sampler,parameters,core,chains}.csv\n",
    "Fits: ", prefix, "{withNa,noNa}_order{0,1}.rds\n",
    "IMPORTANT: Na-specific parameter posteriors are PRIOR-ONLY in noNa fits.\n",
    "Passing MCMC diagnostics is necessary but does not validate the mixing model.\n",
    sep = "")
