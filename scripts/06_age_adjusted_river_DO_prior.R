# 06_age_adjusted_river_DO_prior.R
#
# Does attenuation of river DO fluctuations with groundwater age change
# the inference of an age-associated first-order DO depletion coefficient?
#
# The long-term river DO reference is 7.5 mg/L. The reported ~1.9 mg/L
# variation comes from DAILY MEANS, NOT diurnal oscillations. We do NOT
# interpret 1.9 mg/L as a measured diurnal amplitude.
#
# The river DO PRIOR MEAN is unchanged with MRT. Only the variance changes.
# We split the river daily-mean variance into an uncertain persistent part
# and an assumed annual/seasonal part. An additional unmeasured diurnal
# component is explicitly labelled a sensitivity assumption.
#
# IMPORTANT: A mean residence time does not define travel-time dispersion.
# The attenuation below ASSUMES an exponential travel-time distribution:
#   H(period, MRT) = [1 + (2*pi*MRT/period)^2]^(-1/2)
# which is the harmonic amplitude attenuation for an exponential response.
# Piston flow would not attenuate the harmonic amplitude at all.
#
# Gaussian source-DO random effects are integrated out for speed. This is
# an approximation, and latent source DO is not physically truncated.
#
# Run after scripts/01_three_tracer_mixing.R from project root.
library(tidyverse)
setwd(here::here("Oxygen-Reduction-Joint-Model"))

# ======================================================================
# USER SETTINGS: CHANGE / REVIEW ALL ASSUMPTIONS HERE
# ======================================================================
cfg <- list(
  river_mean = 7.5,                # mg/L
  river_mean_prior_sd = 0.7,      # assumed uncertainty in river mean
  river_daily_mean_sd = 1.9,      # HALF-NORMAL PRIOR SCALE from daily means
  rainfall_upper = 9.5,          # upper bound on effective rain DO mean
  rainfall_variability_sd = 2.5, # half-normal SD prior scale
  residual_variability_sd = 1.5, # half-normal SD prior scale
  DO_measurement_sd = 0.2,       # assumed instrumental / analytical SD
  k_prior_scale = 0.10,          # half-normal scale (year^-1)
  diurnal_sd_assumed = 1.0,      # mg/L; NOT measured from daily means
  year_days = 365.25,
  seasonal_period_years = 1,
  persistent_variance_fractions = c(
    "100% persistent" = 1.00,
    "50% persistent" = 0.50,
    "25% persistent" = 0.25
  ),
  n_chains = 4L,
  warmup = 1400L,
  retained_per_chain = 1800L,
  thin = 2L,
  n_folds = 5L,
  cv_n_chains = 2L,
  cv_warmup = 900L,
  cv_retained = 800L,
  seed = 20261030L,
  make_plots = TRUE,
  plot_font_size = 12,
  plot_dpi = 300
)

dir.create("outputs", showWarnings = FALSE, recursive = TRUE)
dir.create("plots", showWarnings = FALSE, recursive = TRUE)
file <- "outputs/three_tracer_mixing.csv"
if (!file.exists(file)) stop("Run scripts/01_three_tracer_mixing.R first.")
d <- read_csv(
  file, show_col_types = FALSE,
  col_types = cols(Site_ID = col_character(), .default = col_guess())
) %>%
  filter(mixing_reliable,
         is.finite(DO), DO > 0,
         is.finite(MRT), MRT >= 0,
         is.finite(rainfall_contribution),
         is.finite(river_contribution)) %>%
  transmute(
    Site_ID, DO, MRT,
    fR = river_contribution,
    fP = rainfall_contribution
  )
stopifnot(nrow(d) == 29L, n_distinct(d$Site_ID) == 29L,
          all(abs(d$fR + d$fP - 1) < 1e-8))

# ======================================================================
# FREQUENCY-DEPENDENT RIVER DO VARIABILITY
# ======================================================================
# For an exponential travel-time distribution, the amplitude of an
# input sinusoid of period T is attenuated by:
# H = 1 / sqrt(1 + (2*pi*MRT/T)^2).
# The component of variance associated with this sinusoid is H^2 times
# its input variance.
harmonic_attenuation <- function(MRT_years, period_years) {
  1 / sqrt(1 + (2 * pi * MRT_years / period_years)^2)
}

# q is the fraction of DAILY-MEAN variance assumed to remain as an
# age-independent persistent floor. The remaining fraction is assigned
# to a one-year seasonal component. Both are SENSITIVITY ASSUMPTIONS.
#
# An additional diurnal component is included with assumed (not measured)
# input SD. It is not counted within the 1.9 mg/L daily-mean SD.
river_DO_variance <- function(MRT_years, sigma_R, q, prior = cfg) {
  h_year <- harmonic_attenuation(
    MRT_years, prior$seasonal_period_years
  )
  h_day <- harmonic_attenuation(
    MRT_years, 1 / prior$year_days
  )
  sigma_R^2 * (q + (1 - q) * h_year^2) +
    prior$diurnal_sd_assumed^2 * h_day^2
}

# ======================================================================
# MARGINAL BAYESIAN MODEL
# R_i ~ N(mu_R, river_DO_variance(MRT_i))
# P_i ~ N(mu_P, sigma_P^2)
# y_i | R_i,P_i ~ N([fR_i R_i+fP_i P_i] exp(-k*MRT_i),
#                     sqrt(sigma_extra^2+measurement_sd^2))
# k=0 in recharge-only model.
# ======================================================================
unpack <- function(z, use_age) {
  c(
    mu_R = z[1], mu_P = z[2],
    sigma_R = exp(z[3]), sigma_P = exp(z[4]),
    sigma_extra = exp(z[5]),
    k = if (use_age) exp(z[6]) else 0
  )
}

log_posterior <- function(z, data, use_age, q, prior = cfg) {
  p <- unpack(z, use_age)
  if (any(!is.finite(p)) ||
      p["mu_R"] <= 0 || p["mu_R"] >= 12 ||
      p["mu_P"] <= 0 || p["mu_P"] >= prior$rainfall_upper ||
      any(p[c("sigma_R", "sigma_P", "sigma_extra")] >= 30) ||
      p["k"] >= 3) return(-Inf)

  # Half-normal priors on positive SDs and k; exp transforms require
  # log-Jacobian correction. Normal/Uniform prior constants that are
  # independent of parameters cancel in within-model sampling.
  lp <- dnorm(
    p["mu_R"], prior$river_mean, prior$river_mean_prior_sd,
    log = TRUE
  ) +
    dnorm(p["sigma_R"], 0, prior$river_daily_mean_sd,
          log = TRUE) + z[3] +
    dnorm(p["sigma_P"], 0, prior$rainfall_variability_sd,
          log = TRUE) + z[4] +
    dnorm(p["sigma_extra"], 0, prior$residual_variability_sd,
          log = TRUE) + z[5]

  if (use_age) {
    lp <- lp + dnorm(p["k"], 0, prior$k_prior_scale,
                     log = TRUE) + z[6]
  }

  att <- exp(-p["k"] * data$MRT)
  expected <- (
    data$fR * p["mu_R"] + data$fP * p["mu_P"]
  ) * att
  prior_river_variance <- river_DO_variance(
    data$MRT, p["sigma_R"], q, prior
  )
  variance <- (
    data$fR^2 * prior_river_variance +
      data$fP^2 * p["sigma_P"]^2
  ) * att^2 + p["sigma_extra"]^2 + prior$DO_measurement_sd^2
  lp + sum(dnorm(data$DO, expected, sqrt(variance), log = TRUE))
}

# Componentwise adaptive Metropolis. Proposal tuning stops at warmup.
run_chain <- function(data, use_age, q, chain, seed,
                      warmup, n_retained, thin, prior = cfg) {
  set.seed(seed)
  start_P <- c(3, 5, 8, 6)[(chain - 1L) %% 4L + 1L]
  start_k <- c(.03, .07, .12, .05)[(chain - 1L) %% 4L + 1L]
  z <- c(prior$river_mean, start_P,
         log(prior$river_daily_mean_sd),
         log(prior$rainfall_variability_sd),
         log(1))
  if (use_age) z <- c(z, log(start_k))
  step <- c(.27, .60, .22, .22, .20)
  if (use_age) step <- c(step, .23)
  current <- log_posterior(z, data, use_age, q, prior)
  if (!is.finite(current)) stop("Invalid starting state.")
  accepted <- integer(length(z))

  kept <- matrix(NA_real_, nrow = n_retained, ncol = 6)
  colnames(kept) <- c(
    "mu_R", "mu_P", "sigma_R", "sigma_P", "sigma_extra", "k"
  )
  write_row <- 0L
  for (it in seq_len(warmup + n_retained * thin)) {
    for (j in sample.int(length(z))) {
      trial <- z
      trial[j] <- z[j] + rnorm(1, 0, step[j])
      proposed_lp <- log_posterior(trial, data, use_age, q, prior)
      if (is.finite(proposed_lp) &&
          log(runif(1)) < proposed_lp - current) {
        z <- trial
        current <- proposed_lp
        accepted[j] <- accepted[j] + 1L
      }
    }

    if (it <= warmup && it %% 100L == 0L) {
      step <- pmax(
        .015, pmin(
          2, step * exp(pmax(-.25, pmin(.25, accepted / 100 - .34)))
        )
      )
      accepted[] <- 0L
    }

    if (it > warmup && (it - warmup) %% thin == 0L) {
      write_row <- write_row + 1L
      kept[write_row, ] <- unpack(z, use_age)
    }
  }
  as_tibble(kept) %>% mutate(chain = chain)
}

fit_chains <- function(data, use_age, q, seed,
                       n_chains = cfg$n_chains,
                       warmup = cfg$warmup,
                       n_retained = cfg$retained_per_chain) {
  map_dfr(seq_len(n_chains), function(ch) {
    run_chain(data, use_age, q, ch,
              seed + ch * 1009L,
              warmup = warmup, n_retained = n_retained,
              thin = cfg$thin)
  })
}

summarise_posterior <- function(draws, scenario, model) {
  draws %>%
    select(mu_R, mu_P, sigma_R, sigma_P, sigma_extra, k) %>%
    pivot_longer(cols = everything(), names_to = "parameter") %>%
    group_by(parameter) %>%
    summarise(
      median = median(value),
      lower95 = unname(quantile(value, .025)),
      upper95 = unname(quantile(value, .975)),
      .groups = "drop"
    ) %>%
    mutate(scenario = scenario, model = model, .before = 1)
}

basic_rhat <- function(draws, name) {
  chains <- split(draws[[name]], draws$chain)
  n <- min(lengths(chains))
  if (n < 20 || length(chains) < 2) return(NA_real_)
  W <- mean(vapply(chains, function(x) var(x[seq_len(n)]),
                   numeric(1)))
  B <- n * var(vapply(chains, function(x) mean(x[seq_len(n)]),
                      numeric(1)))
  if (W <= 0) return(NA_real_)
  sqrt(((n - 1) * W / n + B / n) / W)
}

posterior_predict <- function(draws, held, q) {
  n <- nrow(draws)
  m <- nrow(held)
  att <- exp(-tcrossprod(draws$k, held$MRT))
  expected <- (
    tcrossprod(draws$mu_R, held$fR) +
      tcrossprod(draws$mu_P, held$fP)
  ) * att

  annual_h2 <- harmonic_attenuation(
    held$MRT, cfg$seasonal_period_years
  )^2
  diel_h2 <- harmonic_attenuation(
    held$MRT, 1 / cfg$year_days
  )^2
  r_var <- tcrossprod(
    draws$sigma_R^2,
    held$fR^2 * (q + (1 - q) * annual_h2)
  ) + matrix(
    cfg$diurnal_sd_assumed^2 * held$fR^2 * diel_h2,
    nrow = n, ncol = m, byrow = TRUE
  )
  variance <- (
    r_var + tcrossprod(draws$sigma_P^2, held$fP^2)
  ) * att^2 + matrix(
    draws$sigma_extra^2 + cfg$DO_measurement_sd^2,
    nrow = n, ncol = m
  )
  list(mean = expected, variance = variance)
}

# ======================================================================
# MATCHED 5-FOLD CV FOR ALL SIX MODELS
# ======================================================================
set.seed(cfg$seed)
fold <- sample(rep(seq_len(cfg$n_folds), length.out = nrow(d)))
folds <- split(seq_len(nrow(d)), fold)

scenario <- tibble(
  label = names(cfg$persistent_variance_fractions),
  q = as.numeric(cfg$persistent_variance_fractions)
)
specs <- map_dfr(seq_len(nrow(scenario)), function(i) {
  tibble(
    label = scenario$label[i],
    q = scenario$q[i],
    use_age = c(FALSE, TRUE)
  )
}) %>%
  mutate(model = if_else(
    use_age, "Age-dependent", "Recharge-only"
  ))

posterior_list <- vector("list", nrow(specs))
posterior_summary <- vector("list", nrow(specs))
CV_results <- vector("list", nrow(specs))
convergence <- vector("list", nrow(specs))

for (i in seq_len(nrow(specs))) {
  s <- specs[i, ]
  cat(
    "Age-adjusted river DO prior:", s$label,
    "|", s$model, "\n"
  )
  draws <- fit_chains(
    d, s$use_age, s$q,
    seed = cfg$seed + 100000L * i
  )
  posterior_list[[i]] <- draws
  posterior_summary[[i]] <- summarise_posterior(
    draws, s$label, s$model
  )
  convergence[[i]] <- tibble(
    scenario = s$label,
    model = s$model,
    k_basic_Rhat = if (s$use_age) basic_rhat(draws, "k") else NA_real_
  )

  CV_results[[i]] <- map_dfr(seq_along(folds), function(j) {
    ids <- folds[[j]]
    train <- d[-ids, , drop = FALSE]
    held <- d[ids, , drop = FALSE]

    train_draws <- fit_chains(
      train, s$use_age, s$q,
      seed = cfg$seed + 700000L + 30000L * i + 1100L * j,
      n_chains = cfg$cv_n_chains,
      warmup = cfg$cv_warmup,
      n_retained = cfg$cv_retained
    )
    pp <- posterior_predict(train_draws, held, s$q)

    map_dfr(seq_len(nrow(held)), function(u) {
      l <- dnorm(held$DO[u], mean = pp$mean[, u],
                 sd = sqrt(pp$variance[, u]), log = TRUE)
      mx <- max(l)
      tibble(
        scenario = s$label, model = s$model, fold = j,
        Site_ID = held$Site_ID[u], MRT = held$MRT[u],
        rainfall_fraction = held$fP[u],
        observed_DO = held$DO[u],
        predicted_DO = mean(pp$mean[, u]),
        log_predictive_density = mx + log(mean(exp(l - mx)))
      )
    })
  })
}

summary_tab <- bind_rows(posterior_summary)
cv_tab <- bind_rows(CV_results)
convergence_tab <- bind_rows(convergence)
performance <- cv_tab %>%
  group_by(scenario, model) %>%
  summarise(
    n = n(),
    RMSE_mg_L = sqrt(mean((observed_DO - predicted_DO)^2)),
    MAE_mg_L = mean(abs(observed_DO - predicted_DO)),
    sum_log_predictive_density = sum(log_predictive_density),
    .groups = "drop"
  )

write_csv(summary_tab, "outputs/age_adjusted_river_prior_posteriors.csv")
write_csv(performance, "outputs/age_adjusted_river_prior_CV_summary.csv")
write_csv(cv_tab, "outputs/age_adjusted_river_prior_CV_predictions.csv")
write_csv(convergence_tab, "outputs/age_adjusted_river_prior_convergence.csv")

if (any(convergence_tab$k_basic_Rhat > 1.05, na.rm = TRUE)) {
  warning("Some k chains may not have converged; increase samples.")
}

# ======================================================================
# FIGURES
# ======================================================================
if (isTRUE(cfg$make_plots)) {
  theme_set(theme_bw(base_size = cfg$plot_font_size))
  k_draws <- map_dfr(
    which(specs$use_age), function(i) {
      posterior_list[[i]] %>%
        mutate(scenario = specs$label[i])
    }
  )
  p_k <- ggplot(k_draws, aes(k, colour = scenario)) +
    geom_density(linewidth = .9) +
    labs(
      title = "First-order DO age coefficient",
      subtitle = "Uncertain river DO with age-adjusted variance",
      x = expression(k~(year^{-1})), y = "Posterior density",
      colour = "Persistent variance"
    )
  print(p_k)
  ggsave("plots/age_adjusted_river_prior_k.png", p_k,
         width = 9, height = 6, dpi = cfg$plot_dpi)

  p_cv <- ggplot(performance, aes(scenario, RMSE_mg_L,
                                  shape = model, group = model)) +
    geom_point(size = 3) +
    geom_line() +
    labs(
      title = "Prediction accuracy under river DO variance attenuation",
      x = "Persistent fraction of river daily-mean DO variance",
      y = "Five-fold held-out DO RMSE (mg/L)",
      shape = "Model"
    )
  print(p_cv)
  ggsave("plots/age_adjusted_river_prior_CV.png", p_cv,
         width = 10, height = 6, dpi = cfg$plot_dpi)
}

cat("\nPosterior age coefficient by scenario:\n")
print(
  summary_tab %>%
    filter(parameter == "k", model == "Age-dependent")
)
cat("\nFive-fold CV performance:\n")
print(performance)
cat("\nBasic R-hat diagnostic:\n")
print(convergence_tab)

cat(
  "\nLIMITATIONS:",
  "\n* The 1.9 mg/L empirical SD is of daily means, NOT the diurnal cycle.",
  "\n* The variance split and diurnal component are assumed; sensitivity only.",
  "\n* An exponential travel-time distribution is imposed, not inferred.",
  "\n* MRT=0 for site 2801 is not proof of zero travel time.",
  "\n* A nonzero k does not locate oxygen reduction within the aquifer.",
  "\n* Source DO distributions are Gaussian, not physically truncated.",
  "\n* Mixing fractions are fixed, and 5-fold CV is not spatial CV.\n"
)
