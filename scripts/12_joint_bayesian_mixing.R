# 12_joint_bayesian_mixing.R
#
# Fully joint Bayesian fit: endmember uncertainty, tracer-informed LSR
# fractions, uncertainty around reported groundwater MRTs, uncertain LSR
# DO, and apparent age-dependent reduction k.
#
# Conservative mixing tracers: chloride, delta-18O, excess AIR.
# Excess N2 is a REACTIVE response, optionally included via an empirical
# association with age-associated oxygen loss. It is NEVER used as an
# ordinary conservative mixing tracer.
#
# Requires cmdstanr + CmdStan. See README.
# Local input only: data/groundwater.csv (NEVER committed).
#
# WARNING: the endmember priors are partially informed by the same
# groundwater samples included in this study (empirical Bayes). They are
# not independent endmember measurements. Review and widen them before
# treating posterior source fractions as definitive.
library(tidyverse)
setwd(here::here("Oxygen-Reduction-Joint-Model"))

# ==================================================================
# USER SETTINGS
# ==================================================================
cfg <- list(
  mrts = c(0.25, 0.50),    # ASSUMED SD of log(true MRT), sensitivity
  use_n2 = c(FALSE, TRUE),# xsN2 as an optional REACTIVE response
  default_gamma = FALSE,  # gamma alpha priors optional; see Stan model
  vary_lsr_do = FALSE,    # TRUE adds per-well hierarchical LSR DO
  gamma_prior_median = 10,
  tracer_prior_scale = 1, # 2 doubles source-endmember uncertainty
  seed = 20261030L,
  warmup = 1100L,
  sampling = 1100L,
  chains = 4L,
  parallel_chains = 4L,
  adapt_delta = 0.94,
  max_treedepth = 12L,
  input = "data/groundwater.csv",
  stan_file = "stan/12_joint_bayesian_mixing.stan"
)

if (!requireNamespace("cmdstanr", quietly = TRUE)) {
  stop("Install cmdstanr (see README) and CmdStan before running.")
}
if (!requireNamespace("posterior", quietly = TRUE)) {
  stop("Install the posterior R package.")
}
if (!file.exists(cfg$input)) {
  stop("Missing private data/groundwater.csv. Run the local data-prep workflow.")
}
if (!file.exists(cfg$stan_file)) {
  stop("Stan file missing: ", cfg$stan_file)
}
tryCatch(
  cmdstanr::cmdstan_path(),
  error = function(e) {
    stop("CmdStan is not installed/configured: ", conditionMessage(e))
  }
)

dir.create("outputs", showWarnings = FALSE, recursive = TRUE)
dir.create("plots", showWarnings = FALSE, recursive = TRUE)

raw <- readr::read_csv(
  cfg$input, show_col_types = FALSE,
  col_types = cols(Site_ID = col_character(),
                   .default = col_guess())
)
stopifnot(
  all(c("Site_ID", "Depth", "DO", "MRT", "O18", "Cl",
        "xs_Air", "xs_Air_uncertainty", "xs_N2",
        "xs_N2_sig") %in% names(raw))
)

d <- raw %>%
  filter(is.finite(DO), DO >= 0,
         is.finite(MRT), MRT >= 0) %>%
  mutate(sample_index = row_number())
stopifnot(
  nrow(d) >= 20,
  all(is.finite(d$DO)),
  all(is.finite(d$MRT))
)
message(
  "Using ", nrow(d), " groundwater samples with DO + MRT. ",
  "Original script 01 used 29 screened samples; this joint model ",
  "does NOT impose the earlier tracer-agreement screen."
)

# Dummy numeric observations are ignored whenever has_* is 0.
observed <- function(x) {
  ifelse(is.finite(x), x, 0)
}
observed_sd <- function(x, fallback) {
  ifelse(is.finite(x) & x > 0, x, fallback)
}

stan_data_base <- list(
  N = nrow(d),
  has_o18 = as.integer(is.finite(d$O18)),
  has_cl = as.integer(is.finite(d$Cl)),
  has_air = as.integer(is.finite(d$xs_Air)),
  has_n2 = as.integer(is.finite(d$xs_N2)),
  o18 = observed(d$O18),
  cl = observed(d$Cl),
  air = observed(d$xs_Air),
  do_obs = as.numeric(d$DO),
  n2_obs = observed(d$xs_N2),
  cl_sd = pmax(1, 0.03 * observed(d$Cl)),
  air_sd = observed_sd(d$xs_Air_uncertainty, .15),
  n2_sd = observed_sd(d$xs_N2_sig, .3),
  mrt_reported = as.numeric(d$MRT)
)

# Run the three main sensitivity scenarios, with same underlying
# measurements and tracer priors. Other combinations can be added.
scenarios <- tribble(
  ~scenario,       ~log_mrt_SD, ~use_n2,
  "MRT025",         cfg$mrts[1], FALSE,
  "MRT025_N2",      cfg$mrts[1], TRUE,
  "MRT050",         cfg$mrts[2], FALSE
)

model <- cmdstanr::cmdstan_model(
  cfg$stan_file,
  quiet = TRUE,
  force_recompile = FALSE
)

summaries <- list()
fraction_summaries <- list()
age_summaries <- list()
fitted_DO <- list()
diagnostics <- list()

for (i in seq_len(nrow(scenarios))) {
  s <- scenarios[i, ]
  message("Fitting joint model: ", s$scenario)
  x <- c(
    stan_data_base,
    list(
      log_mrt_SD = as.numeric(s$log_mrt_SD),
      use_n2 = as.integer(s$use_n2),
      use_gamma = as.integer(cfg$default_gamma),
      vary_lsr_do = as.integer(cfg$vary_lsr_do),
      gamma_prior_median = cfg$gamma_prior_median,
      tracer_prior_scale = cfg$tracer_prior_scale
    )
  )
  fit <- model$sample(
    data = x,
    seed = cfg$seed + i * 101L,
    chains = cfg$chains,
    parallel_chains = cfg$parallel_chains,
    iter_warmup = cfg$warmup,
    iter_sampling = cfg$sampling,
    adapt_delta = cfg$adapt_delta,
    max_treedepth = cfg$max_treedepth,
    refresh = 200
  )

  # Posterior summaries include Rhat and effective sample sizes.
  select_globals <- c(
    "river_DO", "central_LSR_entry_DO",
    "central_LSR_DO_loss", "sample_mean_LSR_DO_loss",
    "k", "sigma_DO",
    "o18_R", "o18_LSR", "cl_R", "cl_LSR",
    "air_R", "air_LSR", "sigma_o18", "sigma_cl", "sigma_air",
    "tau_LSR", "alpha", "n2_R", "n2_LSR",
    "beta_n2", "sigma_n2"
  )
  globals <- fit$summary(
    variables = select_globals
  ) %>%
    mutate(scenario = s$scenario, .before = 1)
  summaries[[i]] <- globals
  diagnostics[[i]] <- fit$diagnostic_summary() %>%
    mutate(scenario = s$scenario, .before = 1)

  site_output <- function(variable) {
    fit$summary(variables = variable) %>%
      mutate(
        sample_index = as.integer(
          str_match(variable, "\\[(\\d+)\\]")[, 2]
        )
      ) %>%
      inner_join(
        d %>% select(sample_index, Site_ID, Depth, MRT, DO),
        by = "sample_index"
      ) %>%
      mutate(scenario = s$scenario, .before = 1)
  }
  fraction_summaries[[i]] <- site_output("f_LSR")
  age_summaries[[i]] <- site_output("true_MRT")
  fitted_DO[[i]] <- site_output("mean_DO")

  # Preserve all posterior draws locally, never in public GitHub.
  fit$save_object(
    file = file.path(
      "outputs",
      paste0("joint_bayes_mixing_", s$scenario, ".rds")
    )
  )
}

all_parameters <- bind_rows(summaries)
all_fractions <- bind_rows(fraction_summaries)
all_ages <- bind_rows(age_summaries)
all_DO <- bind_rows(fitted_DO)
all_diagnostics <- bind_rows(diagnostics)

write_csv(all_parameters, "outputs/joint_bayes_mixing_parameters.csv")
write_csv(all_fractions, "outputs/joint_bayes_mixing_fractions.csv")
write_csv(all_ages, "outputs/joint_bayes_mixing_MRT.csv")
write_csv(all_DO, "outputs/joint_bayes_mixing_fitted_DO.csv")
write_csv(all_diagnostics, "outputs/joint_bayes_mixing_diagnostics.csv")

# Compare old 29-well median source fractions with the FULL joint
# posterior, using Site_ID AND Depth to identify duplicate site IDs.
prior_mix_file <- "outputs/three_tracer_mixing.csv"
if (file.exists(prior_mix_file)) {
  older <- read_csv(
    prior_mix_file, show_col_types = FALSE,
    col_types = cols(Site_ID = col_character(),
                     .default = col_guess())
  )
  if (all(c("Site_ID", "Depth",
            "rainfall_contribution") %in% names(older))) {
    comparison <- all_fractions %>%
      left_join(
        older %>% select(
          Site_ID, Depth, mixing_reliable,
          old_LSR_fraction = rainfall_contribution
        ),
        by = c("Site_ID", "Depth")
      )
    write_csv(
      comparison,
      "outputs/joint_bayes_mixing_vs_old_fractions.csv"
    )
  }
}

# The site summary table uses cmdstanr/posterior default q5/q95:
# these are 90% intervals. Do not describe them as 95% intervals.
pp <- all_parameters %>%
  filter(variable %in% c("k", "central_LSR_DO_loss",
                         "river_DO")) %>%
  ggplot(aes(median, scenario, xmin = q5, xmax = q95)) +
  geom_segment(aes(x = q5, xend = q95,
                   y = scenario, yend = scenario)) +
  geom_point(size = 2.5) +
  facet_wrap(~ variable, scales = "free_x") +
  labs(
    title = "Joint Bayesian parameter sensitivity",
    subtitle = "Posterior median and central 90% intervals",
    x = "Parameter value", y = "Scenario"
  ) +
  theme_bw(base_size = 12)
print(pp)
ggsave(
  "plots/joint_bayes_mixing_parameter_intervals.png",
  pp, width = 12, height = 5, dpi = 300
)

cat("\nGLOBAL POSTERIORS:\n")
print(
  all_parameters %>%
    filter(variable %in% c(
      "k", "central_LSR_DO_loss", "sample_mean_LSR_DO_loss",
      "river_DO", "cl_R", "cl_LSR",
      "o18_R", "o18_LSR", "air_R", "air_LSR"
    )),
  n = Inf
)
cat("\nMCMC DIAGNOSTICS:\n")
print(all_diagnostics, n = Inf)
if (any(
  all_parameters$rhat > 1.05,
  na.rm = TRUE
)) {
  warning(
    "At least one parameter has Rhat > 1.05. ",
    "Do not interpret posterior estimates without more sampling."
  )
}
cat(
  "\nINTERPRETATION:\n",
  "1. Excess AIR is a mixing tracer; excess N2 is not.\n",
  "2. 41 samples are included if the local source data match the study.\n",
  "3. MRT priors and source endmember priors are assumed, not verified.\n",
  "4. The N2 model is an empirical association, not denitrification\n",
  "   stoichiometry. Use the no-N2 model as the primary result.\n",
  "5. No part of this fit proves oxygen consumption happened\n",
  "   exclusively in the saturated regional aquifer.\n",
  "6. Exports are local and ignored by the repository's .gitignore.\n",
  sep = ""
)
