# 04_bayesian_recharge_DO.R
# Bayesian sensitivity test: can uncertain recharge DO explain observed DO
# without a first-order groundwater-age term?
#
# Run AFTER scripts/01_three_tracer_mixing.R from the project directory.
# Requires tidyverse only (MCMC and predictive comparisons use base R).
#
# CRITICAL: This model assumes the distribution of *initial* recharge DO
# is independent of groundwater age, conditional on mixing fraction.
# A positive k does NOT locate reduction in the saturated aquifer.
#
# The well-specific river/rainfall DO distributions are Gaussian, and their
# random effects are analytically integrated out. This is a convenient
# approximation, NOT a physically bounded dissolved-gas model; extreme
# latent negative DO values are possible. Revisit with truncated source
# distributions if needed.
library(tidyverse)

# ======================================================================
# USER SETTINGS
# ======================================================================
cfg <- list(
  river_DO_prior_mean = 7.5,          # Fernhill long-term river mean (mg/L)
  river_mean_prior_sd = 0.7,          # ASSUMED uncertainty in mean, not SD
  river_variability_prior_sd = 1.9,   # Temporal variability scale (mg/L)
  rainfall_DO_upper = 9.5,           # Uniform prior on population mean
  rainfall_variability_prior_sd = 2.5,
  model_discrepancy_prior_sd = 1.5,
  DO_measurement_sd = 0.2,
  k_prior_scale = 0.10,              # Half-normal scale, year^-1
  source_DO_upper_bound = 12,
  n_chains = 4,
  n_posterior = 2200,                # Retained draws per chain
  warmup = 1800,
  thin = 2,
  n_folds = 5,
  cv_draws = 1000,
  cv_warmup = 1200,
  cv_thin = 2,
  make_plots = TRUE,
  plot_font_size = 13,
  plot_dpi = 300,
  seed = 20261030
)
# Broad-prior stress test: not presented as "the" correct prior.
wide_prior <- modifyList(cfg, list(
  river_mean_prior_sd = 1.5,
  river_variability_prior_sd = 3.5,
  rainfall_variability_prior_sd = 3.5
))

# ======================================================================
# DATA: SAME 29-WELL COHORT AS THE JOINT MODEL
# ======================================================================
file <- "outputs/three_tracer_mixing.csv"
if (!file.exists(file)) stop("Run scripts/01_three_tracer_mixing.R first.")
dir.create("outputs", showWarnings = FALSE)
dir.create("plots", showWarnings = FALSE)

d <- read_csv(
  file, show_col_types = FALSE,
  col_types = cols(Site_ID = col_character(), .default = col_guess())
) %>%
  filter(mixing_reliable, is.finite(DO), DO > 0,
         is.finite(MRT), MRT >= 0,
         is.finite(rainfall_contribution),
         is.finite(river_contribution)) %>%
  transmute(
    Site_ID, site, DO, MRT,
    fP = rainfall_contribution,
    fR = river_contribution
  )
stopifnot(nrow(d) == 29, n_distinct(d$Site_ID) == 29,
          all(abs(d$fR + d$fP - 1) < 1e-8))

# ======================================================================
# STATISTICAL MODEL (LATENT SOURCE DO ANALYTICALLY MARGINALISED)
# ======================================================================
#  R_i ~ Normal(mu_R, sigma_R)
#  P_i ~ Normal(mu_P, sigma_P)
#  y_i = (fR_i * R_i + fP_i * P_i) * exp(-k*MRT_i) + epsilon_i
#  epsilon_i ~ Normal(0, sqrt(sigma_extra^2 + DO_measurement_sd^2))
#
#  mu_R ~ Normal(7.5, 0.7), bounded 0..12
#  mu_P ~ Uniform(0, 9.5)
#  sigma_R ~ HalfNormal(1.9)
#  sigma_P ~ HalfNormal(2.5)
#  sigma_extra ~ HalfNormal(1.5)
#  k ~ HalfNormal(0.1) in the age-dependent model, otherwise exactly zero.
#
# After marginalisation:
#  E[y_i] = (fR_i*mu_R + fP_i*mu_P) exp(-k*MRT_i)
#  Var[y_i] = (fR_i^2 sigma_R^2 + fP_i^2 sigma_P^2)
#              exp(-2k*MRT_i) + sigma_extra^2 + measurement_sd^2

unpack <- function(z, use_k) {
  c(
    mu_R = z[1], mu_P = z[2],
    sigma_R = exp(z[3]), sigma_P = exp(z[4]),
    sigma_extra = exp(z[5]),
    k = if (use_k) exp(z[6]) else 0
  )
}

log_posterior <- function(z, data, use_k, prior) {
  p <- unpack(z, use_k)
  if (any(!is.finite(p)) ||
      p["mu_R"] <= 0 ||
      p["mu_R"] >= prior$source_DO_upper_bound ||
      p["mu_P"] <= 0 ||
      p["mu_P"] >= prior$rainfall_DO_upper ||
      any(p[c("sigma_R", "sigma_P", "sigma_extra")] >= 30) ||
      p["k"] >= 3) return(-Inf)

  # Priors on SDs are half-normals. Add log-Jacobians for log(SD).
  lp <- dnorm(
    p["mu_R"], prior$river_DO_prior_mean,
    prior$river_mean_prior_sd, log = TRUE
  ) +
    dnorm(p["sigma_R"], 0, prior$river_variability_prior_sd,
          log = TRUE) + z[3] +
    dnorm(p["sigma_P"], 0, prior$rainfall_variability_prior_sd,
          log = TRUE) + z[4] +
    dnorm(p["sigma_extra"], 0, prior$model_discrepancy_prior_sd,
          log = TRUE) + z[5]

  if (use_k) {
    lp <- lp +
      dnorm(p["k"], 0, prior$k_prior_scale, log = TRUE) + z[6]
  }

  atten <- exp(-p["k"] * data$MRT)
  expected <- (data$fR * p["mu_R"] + data$fP * p["mu_P"]) * atten
  variance <- (data$fR^2 * p["sigma_R"]^2 +
               data$fP^2 * p["sigma_P"]^2) *
    atten^2 + p["sigma_extra"]^2 + prior$DO_measurement_sd^2
  lp + sum(dnorm(data$DO, expected, sqrt(variance), log = TRUE))
}

# Componentwise adaptive Metropolis-Hastings in transformed parameter space.
# Adaptation takes place ONLY during warmup; retained draws are stationary.
run_chain <- function(data, use_k, prior, seed,
                      n_draw = prior$n_posterior,
                      warmup = prior$warmup, thin = prior$thin,
                      chain = 1L) {
  set.seed(seed)
  muP_start <- c(3.0, 5.0, 8.0, 4.5)[(chain - 1L) %% 4L + 1L]
  k_start <- c(0.03, 0.07, 0.12, 0.05)[(chain - 1L) %% 4L + 1L]
  z <- c(
    prior$river_DO_prior_mean, muP_start,
    log(prior$river_variability_prior_sd),
    log(prior$rainfall_variability_prior_sd),
    log(1.0)
  )
  if (use_k) z <- c(z, log(k_start))

  step <- c(0.27, 0.60, 0.22, 0.22, 0.20)
  if (use_k) step <- c(step, 0.23)
  npar <- length(z)
  current <- log_posterior(z, data, use_k, prior)
  n_iter <- warmup + n_draw * thin
  retained <- matrix(NA_real_, nrow = n_draw, ncol = 6)
  colnames(retained) <- c("mu_R", "mu_P", "sigma_R",
                          "sigma_P", "sigma_extra", "k")
  accepted <- integer(npar)
  counter <- 0L

  for (iteration in seq_len(n_iter)) {
    for (j in sample.int(npar)) {
      proposal <- z
      proposal[j] <- z[j] + rnorm(1, 0, step[j])
      candidate <- log_posterior(proposal, data, use_k, prior)
      if (is.finite(candidate) &&
          log(runif(1)) < candidate - current) {
        z <- proposal
        current <- candidate
        accepted[j] <- accepted[j] + 1L
      }
    }

    if (iteration <= warmup && iteration %% 100L == 0L) {
      accept_rate <- accepted / 100
      step <- pmax(
        0.015, pmin(
          2,
          step * exp(pmax(-0.25, pmin(0.25, accept_rate - 0.34)))
        )
      )
      accepted[] <- 0L
    }

    if (iteration > warmup &&
        (iteration - warmup) %% thin == 0L) {
      counter <- counter + 1L
      retained[counter, ] <- unpack(z, use_k)
    }
  }

  as_tibble(retained) %>%
    mutate(chain = chain, iteration = seq_len(n_draw))
}

fit_chains <- function(data, use_k, prior, offset = 0L) {
  map_dfr(seq_len(prior$n_chains), function(chain) {
    run_chain(
      data, use_k, prior,
      seed = prior$seed + offset + 1009L * chain,
      chain = chain
    )
  })
}

# Basic between-chain R-hat; this is not split/rank-normalised R-hat.
rhat_basic <- function(draws, parameter) {
  splits <- split(draws[[parameter]], draws$chain)
  n <- min(vapply(splits, length, integer(1)))
  if (length(splits) < 2L || n < 20L) return(NA_real_)
  means <- vapply(splits, function(x) mean(x[seq_len(n)]), numeric(1))
  v <- vapply(splits, function(x) var(x[seq_len(n)]), numeric(1))
  W <- mean(v)
  B <- n * var(means)
  sqrt(((n - 1) / n * W + B / n) / W)
}

summarize_draws <- function(draws, model_name, prior_name) {
  draws %>%
    select(mu_R, mu_P, sigma_R, sigma_P, sigma_extra, k) %>%
    pivot_longer(everything(), names_to = "parameter",
                 values_to = "value") %>%
    group_by(parameter) %>%
    summarise(
      median = median(value),
      lower_95 = unname(quantile(value, 0.025)),
      upper_95 = unname(quantile(value, 0.975)),
      .groups = "drop"
    ) %>%
    mutate(
      model = model_name, prior_scenario = prior_name,
      .before = 1
    )
}

# ======================================================================
# FIT FULL MODELS AND PRIOR SENSITIVITY
# ======================================================================
cat("Fitting 4-chain Bayesian models. This may take a few minutes in R.\n")
no_k <- fit_chains(d, FALSE, cfg, offset = 0L)
with_k <- fit_chains(d, TRUE, cfg, offset = 100000L)
wide_k <- fit_chains(d, TRUE, wide_prior, offset = 200000L)

summary_table <- bind_rows(
  summarize_draws(no_k, "Recharge only (k=0)", "Baseline"),
  summarize_draws(with_k, "Recharge + age effect", "Baseline"),
  summarize_draws(wide_k, "Recharge + age effect", "Wide priors")
)
write_csv(summary_table, "outputs/bayesian_posterior_summary.csv")

convergence <- tibble(
  scenario = c("Baseline", "Wide priors"),
  k_Rhat_basic = c(
    rhat_basic(with_k, "k"), rhat_basic(wide_k, "k")
  )
)
write_csv(convergence, "outputs/bayesian_convergence.csv")
if (any(convergence$k_Rhat_basic > 1.05)) {
  warning("k may not have converged: increase warmup and n_posterior.")
}

# ======================================================================
# HELD-OUT PREDICTION: SAME FIVE FOLDS FOR BOTH MODELS
# ======================================================================
posterior_predictive <- function(draws, holdout, prior) {
  n <- nrow(draws)
  m <- nrow(holdout)
  attenuation <- exp(-outer(draws$k, holdout$MRT))
  expected <- (
    tcrossprod(draws$mu_R, holdout$fR) +
      tcrossprod(draws$mu_P, holdout$fP)
  ) * attenuation
  variance <- (
    tcrossprod(draws$sigma_R^2, holdout$fR^2) +
      tcrossprod(draws$sigma_P^2, holdout$fP^2)
  ) * attenuation^2 +
    matrix(
      draws$sigma_extra^2 + prior$DO_measurement_sd^2,
      nrow = n, ncol = m
    )

  map_dfr(seq_len(m), function(j) {
    logdens <- dnorm(
      holdout$DO[j], expected[, j],
      sqrt(variance[, j]), log = TRUE
    )
    mx <- max(logdens)
    tibble(
      Site_ID = holdout$Site_ID[j],
      site = holdout$site[j],
      MRT = holdout$MRT[j],
      rainfall_fraction = holdout$fP[j],
      observed_DO = holdout$DO[j],
      predicted_DO = mean(expected[, j]),
      log_predictive_density = mx + log(mean(exp(logdens - mx)))
    )
  })
}

set.seed(cfg$seed)
fold <- sample(rep(seq_len(cfg$n_folds), length.out = nrow(d)))
folds <- split(seq_len(nrow(d)), fold)

cv_results <- map_dfr(c(FALSE, TRUE), function(use_k) {
  map_dfr(seq_along(folds), function(j) {
    test_ids <- folds[[j]]
    train <- d[-test_ids, , drop = FALSE]
    held <- d[test_ids, , drop = FALSE]
    cat(
      "CV:", if (use_k) "age-dependent" else "recharge-only",
      "fold", j, "of", length(folds), "\n"
    )

    # Two independent chains for each CV fold.
    temp_cfg <- modifyList(cfg, list(
      n_chains = 2L,
      n_posterior = cfg$cv_draws,
      warmup = cfg$cv_warmup,
      thin = cfg$cv_thin
    ))
    fit <- fit_chains(
      train, use_k, temp_cfg,
      offset = 300000L + 10000L * as.integer(use_k) + j * 101L
    )

    posterior_predictive(fit, held, cfg) %>%
      mutate(
        model = if (use_k) "Age-dependent" else "Recharge-only",
        fold = j, .before = 1
      )
  })
})
write_csv(cv_results, "outputs/bayesian_CV_predictions.csv")

cv_summary <- cv_results %>%
  group_by(model) %>%
  summarise(
    n = n(),
    RMSE_mg_L = sqrt(mean((observed_DO - predicted_DO)^2)),
    MAE_mg_L = mean(abs(observed_DO - predicted_DO)),
    expected_log_predictive_density = sum(log_predictive_density),
    .groups = "drop"
  )
write_csv(cv_summary, "outputs/bayesian_CV_summary.csv")

# ======================================================================
# PLOTS: POSTERIOR k AND HELD-OUT PREDICTIONS
# ======================================================================
if (isTRUE(cfg$make_plots)) {
  theme_set(theme_bw(base_size = cfg$plot_font_size))

  posterior_k <- bind_rows(
    with_k %>% mutate(prior_scenario = "Baseline"),
    wide_k %>% mutate(prior_scenario = "Wide")
  )
  p_k <- ggplot(
    posterior_k, aes(k, colour = prior_scenario, fill = prior_scenario)
  ) +
    geom_density(alpha = 0.2, linewidth = 0.9) +
    labs(
      title = "Posterior age-dependent DO reduction rate",
      subtitle = "Hierarchical uncertain river and rainfall DO model",
      x = expression(k~(year^{-1})),
      y = "Posterior density",
      colour = "Prior", fill = "Prior"
    )
  print(p_k)
  ggsave(
    "plots/bayesian_posterior_k.png", p_k,
    width = 9, height = 6, dpi = cfg$plot_dpi
  )

  p_cv <- ggplot(
    cv_results, aes(observed_DO, predicted_DO)
  ) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                colour = "grey45") +
    geom_point(aes(colour = rainfall_fraction), size = 3) +
    facet_wrap(~ model) +
    coord_equal() +
    scale_colour_viridis_c(labels = scales::percent) +
    labs(
      title = "Held-out DO predictions",
      subtitle = "Five-fold cross-validation; no observation used in own fit",
      x = "Observed DO (mg/L)",
      y = "Predicted DO (mg/L)",
      colour = "Rainfall contribution"
    )
  print(p_cv)
  ggsave(
    "plots/bayesian_CV_observed_predicted.png", p_cv,
    width = 11, height = 6, dpi = cfg$plot_dpi
  )

  p_resid <- cv_results %>%
    mutate(residual = observed_DO - predicted_DO) %>%
    ggplot(aes(MRT, residual)) +
    geom_hline(yintercept = 0, colour = "grey45",
               linetype = "dashed") +
    geom_point(aes(colour = rainfall_fraction), size = 3) +
    facet_wrap(~ model) +
    scale_colour_viridis_c(labels = scales::percent) +
    labs(
      title = "Held-out prediction errors versus groundwater age",
      x = "Mean residence time (years)",
      y = "Observed minus predicted DO (mg/L)",
      colour = "Rainfall contribution"
    )
  print(p_resid)
  ggsave(
    "plots/bayesian_CV_residuals.png", p_resid,
    width = 11, height = 6, dpi = cfg$plot_dpi
  )
}

cat("\nPosterior summary:\n")
print(summary_table)
cat("\nConvergence check (basic, non-rank R-hat):\n")
print(convergence)
cat("\nFive-fold predictive comparison:\n")
print(cv_summary)
cat(
  "\nINTERPRETATION: A posterior k above zero shows the need for an",
  "\nage-associated predictor GIVEN assumed age-independent source DO.",
  "\nIt does not establish the physical location of oxygen consumption.\n"
)
