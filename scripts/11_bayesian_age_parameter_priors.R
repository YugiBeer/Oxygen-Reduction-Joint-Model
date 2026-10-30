# 11_bayesian_age_parameter_priors.R
#
# Bayesian uncertainty analysis for the GAMMA age-distribution shape
# alpha, apparent first-order oxygen-reduction coefficient k, effective
# land-surface recharge (LSR) DO, mean river DO, and residual SD.
#
# DO_i = [(1-fLSR_i)*R + fLSR_i*P] *
#        (1 + k*MRT_i/alpha)^(-alpha) + epsilon_i
# Piston flow is fitted separately as exp(-k*MRT_i).
#
# We do NOT introduce a river-recharge oxygen-loss parameter.
# LSR means LAND-SURFACE recharge (rainfall through soil/unsaturated zone).
# This analysis does not infer where age-associated reduction occurred.
#
# Important: ONLY the 29 wells passing the three-tracer QC screen are
# used. Their mixing fractions and reported MRTs are treated as given.
# We are testing prior uncertainty, NOT independently dating groundwater.
#
# IMPORTANT PRIOR ASSUMPTIONS:
#   R       ~ Normal(7.5, 0.7), truncated (0,12) mg/L
#   P       ~ Uniform(0,9.5) mg/L (effective LSR DO at aquifer entry)
#   k       ~ HalfNormal(0.2) year^-1
#   sigma   ~ HalfNormal(1.5) mg/L
#   ln(alpha) ~ Normal(log(median), 1.3), alpha in (0.25,500)
#
# The 0.7 river-mean uncertainty, k prior, residual dispersion, gamma
# shape priors and DO measurement SD=0.2 are ASSUMPTIONS, not measured
# errors. The gamma prior medians 2/10/50 are sensitivity scenarios.
#
# Uses simple adaptive component-wise Metropolis-Hastings, four chains
# with warmup-only adaptation. No extra MCMC packages required.
#
# Run after scripts/01_three_tracer_mixing.R in Positron.
library(tidyverse)
setwd(here::here("Oxygen-Reduction-Joint-Model"))

# ===============================================================
# SETTINGS
# ===============================================================
cfg <- list(
  river_DO_prior_mean = 7.5,      # mg/L, monitored river reference
  river_DO_prior_SD = 0.7,        # mg/L, assumed uncertainty in mean
  max_LSR_entry_DO = 9.5,        # mg/L, reference before LSR loss
  k_prior_SD = 0.2,              # year^-1, HalfNormal
  extra_residual_SD_prior = 1.5, # mg/L, HalfNormal
  measurement_SD = 0.2,          # mg/L, assumed
  gamma_prior_medians = c(2, 10, 50),
  gamma_log_prior_SD = 1.3,
  gamma_shape_bounds = c(0.25, 500),
  full_n_chains = 4L,
  full_warmup = 1700L,
  full_keep_per_chain = 1800L,
  cv_n_chains = 2L,
  cv_warmup = 1000L,
  cv_keep_per_chain = 1100L,
  thin = 2L,
  n_cv_folds = 5L,
  exclude_2801 = TRUE,
  seed = 20261030L,
  plot_size = 12L,
  plot_dpi = 300L
)

dir.create("outputs", showWarnings = FALSE, recursive = TRUE)
dir.create("plots", showWarnings = FALSE, recursive = TRUE)
input <- "outputs/three_tracer_mixing.csv"
if (!file.exists(input)) {
  stop("Run scripts/01_three_tracer_mixing.R first.")
}

d <- read_csv(
  input, show_col_types = FALSE,
  col_types = cols(
    Site_ID = col_character(),
    .default = col_guess()
  )
) %>%
  filter(
    mixing_reliable,
    is.finite(DO), DO > 0,
    is.finite(MRT), MRT >= 0,
    is.finite(rainfall_contribution),
    is.finite(river_contribution)
  ) %>%
  transmute(
    Site_ID, DO, MRT,
    fLSR = rainfall_contribution,
    fRiver = river_contribution
  ) %>%
  mutate(cv_group = case_when(
    str_detect(Site_ID, "Eastbourne") ~ "Eastbourne",
    str_detect(Site_ID, "Gilbertson") ~ "Gilbertson",
    str_detect(Site_ID, "Flaxmere") ~ "Flaxmere",
    TRUE ~ Site_ID
  ))
stopifnot(
  nrow(d) == 29L,
  n_distinct(d$Site_ID) == 29L,
  sum(d$fLSR >= .5) == 5L,
  all(abs(d$fLSR + d$fRiver - 1) < 1e-8)
)

# ===============================================================
# BAYESIAN LIKELIHOOD AND PRIOR
# ===============================================================
# Sampling state z:
#   R, P, log(k), log(sigma), [log(alpha) if gamma model]
# Jacobians +log(k) and +log(sigma) transform the physical-space
# HalfNormal priors. alpha is specified as normal IN log(alpha),
# so its prior needs no further Jacobian in the sampling coordinates.
#
# The source-weighted input DO is held constant through the age
# distribution. In reality the river and LSR components can have
# different age distributions, and redox kinetics may be nonuniform.

attenuation <- function(MRT, k, alpha) {
  if (is.infinite(alpha)) {
    exp(-k * MRT)
  } else {
    exp(-alpha * log1p(k * MRT / alpha))
  }
}

parameters <- function(z, piston) {
  list(
    R = z[1],
    P = z[2],
    k = exp(z[3]),
    sigma = exp(z[4]),
    alpha = if (piston) Inf else exp(z[5])
  )
}

log_posterior <- function(z, data, piston, prior_median) {
  if (any(!is.finite(z))) return(-Inf)
  R <- z[1]
  P <- z[2]
  logk <- z[3]
  logs <- z[4]
  if (R <= 0 || R >= 12 || P <= 0 ||
      P >= cfg$max_LSR_entry_DO ||
      logk <= -9 || logk >= 1 ||
      logs <= -8 || logs >= 3) return(-Inf)

  k <- exp(logk)
  sigma <- exp(logs)

  if (piston) {
    alpha <- Inf
    prior_alpha <- 0
  } else {
    logalpha <- z[5]
    if (logalpha <= log(cfg$gamma_shape_bounds[1]) ||
        logalpha >= log(cfg$gamma_shape_bounds[2])) {
      return(-Inf)
    }
    alpha <- exp(logalpha)
    prior_alpha <- dnorm(
      logalpha, log(prior_median),
      cfg$gamma_log_prior_SD, log = TRUE
    )
  }

  att <- attenuation(data$MRT, k, alpha)
  mu <- (data$fRiver * R + data$fLSR * P) * att
  sd_total <- sqrt(sigma^2 + cfg$measurement_SD^2)

  prior <- dnorm(
    R, cfg$river_DO_prior_mean, cfg$river_DO_prior_SD,
    log = TRUE
  ) +
    # Constant Uniform P prior omitted because it cancels in MH.
    dnorm(k, 0, cfg$k_prior_SD, log = TRUE) +
    logk +
    dnorm(sigma, 0, cfg$extra_residual_SD_prior,
          log = TRUE) +
    logs + prior_alpha

  sum(dnorm(data$DO, mu, sd_total, log = TRUE)) + prior
}

starting_state <- function(piston, chain, prior_median) {
  i <- (chain - 1L) %% 4L + 1L
  z <- c(
    R = cfg$river_DO_prior_mean + (i - 2.5) * .12,
    P = c(8, 6, 7, 8.5)[i],
    logk = log(c(.055, .07, .09, .11)[i]),
    logsigma = log(c(.7, .9, 1.2, 1.5)[i])
  )
  if (!piston) {
    z <- c(z, logalpha = log(
      pmin(300, pmax(.5,
        prior_median * c(1/3, 1, 2, 4)[i]))
    ))
  }
  z
}

sample_one_chain <- function(data, piston, prior_median,
                             chain, seed, warmup, draws, thin) {
  set.seed(seed)
  z <- starting_state(piston, chain, prior_median)
  step <- c(.19, .38, .17, .15)
  if (!piston) step <- c(step, .24)
  lp <- log_posterior(z, data, piston, prior_median)
  if (!is.finite(lp)) stop("Invalid MCMC start")

  accepted <- numeric(length(z))
  out <- matrix(
    NA_real_, nrow = draws, ncol = length(z),
    dimnames = list(NULL, names(z))
  )
  written <- 0L
  total <- warmup + draws * thin

  for (it in seq_len(total)) {
    for (j in sample.int(length(z))) {
      trial <- z
      trial[j] <- trial[j] + rnorm(1, sd = step[j])
      ltrial <- log_posterior(
        trial, data, piston, prior_median
      )
      if (is.finite(ltrial) &&
          log(runif(1)) < ltrial - lp) {
        z <- trial
        lp <- ltrial
        accepted[j] <- accepted[j] + 1
      }
    }

    if (it <= warmup && it %% 100L == 0L) {
      step <- pmax(.025, pmin(2,
        step * exp(pmax(-.25,
                        pmin(.25, accepted / 100 - .34)))
      ))
      accepted[] <- 0
    }

    if (it > warmup && (it - warmup) %% thin == 0L) {
      written <- written + 1L
      out[written, ] <- z
    }
  }
  as_tibble(out) %>% mutate(chain = chain)
}

fit_posterior <- function(data, piston, prior_median,
                          n_chains, warmup, draws, seed) {
  bind_rows(lapply(seq_len(n_chains), function(chain) {
    sample_one_chain(
      data, piston, prior_median, chain,
      seed + 104729L * chain,
      warmup, draws, cfg$thin
    )
  }))
}

# Basic split-Rhat diagnostic on untransformed sampling state. This
# checks between-chain agreement but NOT effective sample size.
split_rhat <- function(draws, variable) {
  chain_data <- split(draws[[variable]], draws$chain)
  n <- min(lengths(chain_data))
  n <- n - n %% 2L
  if (n < 40L) return(NA_real_)

  halves <- unlist(lapply(chain_data, function(x) {
    list(x[seq_len(n/2)],
         x[(n/2 + 1):n])
  }), recursive = FALSE)
  within <- mean(vapply(halves, var, numeric(1)))
  between <- (n/2) * var(vapply(
    halves, mean, numeric(1)
  ))
  if (!is.finite(within) || within <= 0) return(NA_real_)
  sqrt(
    (((n/2)-1) / (n/2) * within + between / (n/2)) / within
  )
}

describe_posterior <- function(draws, piston) {
  value <- function(name) {
    if (name == "river_DO") return(draws$R)
    if (name == "LSR_entry_DO") return(draws$P)
    if (name == "LSR_DO_loss") return(cfg$max_LSR_entry_DO - draws$P)
    if (name == "k") return(exp(draws$logk))
    if (name == "extra_residual_SD") return(exp(draws$logsigma))
    if (name == "alpha") return(exp(draws$logalpha))
    if (name == "age_coefficient_of_variation") {
      return(exp(-draws$logalpha/2))
    }
    stop("Unknown posterior field")
  }
  fields <- c(
    "river_DO", "LSR_entry_DO", "LSR_DO_loss",
    "k", "extra_residual_SD"
  )
  if (!piston) fields <- c(
    fields, "alpha", "age_coefficient_of_variation"
  )
  map_dfr(fields, function(f) {
    v <- value(f)
    tibble(
      parameter = f,
      median = median(v),
      lower95 = unname(quantile(v, .025)),
      upper95 = unname(quantile(v, .975))
    )
  })
}

predict_from_draws <- function(draws, held, piston) {
  # Posterior mean expectation of DO (used for squared-error CV).
  # The likelihood also contains residual uncertainty; this is not
  # a posterior predictive confidence interval.
  z <- as.matrix(draws[, c("R", "P", "logk", "logsigma")])
  alpha <- if (piston) {
    rep(Inf, nrow(z))
  } else {
    exp(draws$logalpha)
  }
  as.numeric(vapply(seq_len(nrow(held)), function(i) {
    k <- exp(z[, "logk"])
    att <- if (piston) {
      exp(-k * held$MRT[i])
    } else {
      exp(-alpha * log1p(k * held$MRT[i] / alpha))
    }
    mean(
      (held$fRiver[i] * z[, "R"] +
         held$fLSR[i] * z[, "P"]) * att
    )
  }, numeric(1)))
}

# ===============================================================
# POSTERIORS: PISTON; GAMMA WITH 3 ALTERNATIVE ALPHA PRIORS;
# AND EXCLUDING SITE 2801 UNDER INTERMEDIATE PRIOR.
# ===============================================================
run_specs <- tibble(
  model = c(
    "Piston Bayesian",
    "Gamma prior median 2",
    "Gamma prior median 10",
    "Gamma prior median 50",
    "Piston Bayesian",
    "Gamma prior median 10"
  ),
  piston = c(TRUE, FALSE, FALSE, FALSE, TRUE, FALSE),
  gamma_prior_median = c(
    NA_real_, 2, 10, 50, NA_real_, 10
  ),
  selection = c(
    rep("All 29", 4), rep("Excluding 2801", 2)
  )
)
all_draws <- list()
all_summaries <- list()
all_diagnostics <- list()

for (i in seq_len(nrow(run_specs))) {
  spec <- run_specs[i, ]
  sub <- if (spec$selection == "All 29") {
    d
  } else {
    d %>% filter(Site_ID != "2801")
  }
  median_alpha <- if (spec$piston) {
    10
  } else spec$gamma_prior_median

  message(
    "Posterior ", spec$model, " / ", spec$selection,
    " (", nrow(sub), " samples)"
  )
  draws <- fit_posterior(
    sub, spec$piston, median_alpha,
    cfg$full_n_chains, cfg$full_warmup,
    cfg$full_keep_per_chain,
    seed = cfg$seed + 500L * i
  )
  all_draws[[i]] <- draws
  all_summaries[[i]] <- describe_posterior(
    draws, spec$piston
  ) %>%
    mutate(
      model = spec$model,
      selection = spec$selection,
      .before = 1
    )

  rhats <- map_dfr(
    c("R", "P", "logk", "logsigma",
      if (!spec$piston) "logalpha"),
    function(name) {
      tibble(
        parameter = name,
        split_Rhat = split_rhat(draws, name)
      )
    }
  )
  all_diagnostics[[i]] <- rhats %>%
    mutate(
      model = spec$model,
      selection = spec$selection,
      .before = 1
    )
}

posterior <- bind_rows(all_summaries)
diagnostics <- bind_rows(all_diagnostics)
write_csv(
  posterior,
  "outputs/bayesian_age_parameter_posterior_summary.csv"
)
write_csv(
  diagnostics,
  "outputs/bayesian_age_parameter_convergence.csv"
)

# Print a warning if chains show significant nonconvergence.
if (any(diagnostics$split_Rhat > 1.05, na.rm = TRUE)) {
  warning(
    "Split Rhat > 1.05: extend warmup and sampling before interpreting."
  )
}

# Prior dependence of gamma shape and its correlation with apparent k.
k_alpha_association <- map_dfr(c(2, 3, 4, 6), function(i) {
  spec <- run_specs[i, ]
  draws <- all_draws[[i]]
  tibble(
    selection = spec$selection,
    gamma_prior_median = spec$gamma_prior_median,
    posterior_logalpha_logk_correlation = cor(
      draws$logalpha, draws$logk
    )
  )
})
write_csv(
  k_alpha_association,
  "outputs/bayesian_age_parameter_correlations.csv"
)

# ===============================================================
# GROUPED HELD-OUT PREDICTION: PISTON VS BAYESIAN GAMMA.
# Keep related named well sets in the same validation fold.
# All 29 wells; posterior uncertainty in all parameters is averaged.
# ===============================================================
set.seed(cfg$seed)
grouped <- split(seq_len(nrow(d)), d$cv_group)
order_groups <- sample(seq_along(grouped))
order_groups <- order_groups[
  order(-lengths(grouped)[order_groups])
]
folds <- vector("list", cfg$n_cv_folds)
for (g in order_groups) {
  j <- which.min(lengths(folds))
  folds[[j]] <- c(folds[[j]], grouped[[g]])
}
stopifnot(
  all(lengths(folds) > 0),
  identical(sort(unlist(folds)), seq_len(nrow(d)))
)

cv_specs <- run_specs[1:4, , drop = FALSE]
cv_predictions <- map_dfr(seq_len(nrow(cv_specs)), function(i) {
  spec <- cv_specs[i, ]
  median_alpha <- if (spec$piston) {
    10
  } else spec$gamma_prior_median

  map_dfr(seq_along(folds), function(j) {
    held <- d[folds[[j]], , drop = FALSE]
    train <- d[-folds[[j]], , drop = FALSE]
    message(
      "CV ", spec$model, " fold ", j, "/",
      length(folds), " (held out ", nrow(held), ")"
    )
    draws <- fit_posterior(
      train, spec$piston, median_alpha,
      cfg$cv_n_chains, cfg$cv_warmup,
      cfg$cv_keep_per_chain,
      seed = cfg$seed + 10000L * i + 100L * j
    )
    tibble(
      model = spec$model,
      fold = j,
      Site_ID = held$Site_ID,
      observed_DO = held$DO,
      predicted_DO = predict_from_draws(
        draws, held, spec$piston
      ),
      LSR_fraction = held$fLSR,
      MRT = held$MRT
    )
  })
})
cv_summary <- cv_predictions %>%
  group_by(model) %>%
  summarise(
    n = n(),
    CV_RMSE_mg_L = sqrt(mean(
      (observed_DO - predicted_DO)^2
    )),
    CV_MAE_mg_L = mean(abs(observed_DO - predicted_DO)),
    .groups = "drop"
  )
write_csv(
  cv_predictions,
  "outputs/bayesian_age_parameter_CV_predictions.csv"
)
write_csv(
  cv_summary,
  "outputs/bayesian_age_parameter_CV_summary.csv"
)

# ===============================================================
# VISUAL DIAGNOSTICS
# ===============================================================
theme_set(theme_bw(base_size = cfg$plot_size))
k_intervals <- posterior %>%
  filter(
    selection == "All 29",
    parameter == "k"
  ) %>%
  mutate(model = forcats::fct_reorder(model, median))
p_k <- ggplot(
  k_intervals,
  aes(median, model, xmin = lower95, xmax = upper95)
) +
  geom_segment(aes(x = lower95, xend = upper95,
                   y = model, yend = model), linewidth = 1) +
  geom_point(size = 2.5) +
  labs(
    title = "Age-model prior sensitivity of apparent k",
    subtitle = "Bayesian posterior medians and 95% credible intervals",
    x = expression(k~(year^{-1})),
    y = NULL
  )
print(p_k)
ggsave(
  "plots/bayesian_age_parameter_k_intervals.png",
  p_k, width = 9, height = 5.5, dpi = cfg$plot_dpi
)

alpha_data <- map_dfr(2:4, function(i) {
  tibble(
    model = run_specs$model[i],
    logalpha = all_draws[[i]]$logalpha
  )
})
p_alpha <- ggplot(
  alpha_data, aes(logalpha, colour = model)
) +
  geom_density(linewidth = 1) +
  labs(
    title = "Gamma shape posterior depends on the prior",
    subtitle = "Sampled gamma shape is not independently measured",
    x = "log(gamma shape alpha)",
    y = "Posterior density",
    colour = "Gamma prior"
  )
print(p_alpha)
ggsave(
  "plots/bayesian_age_parameter_alpha_prior_sensitivity.png",
  p_alpha, width = 9, height = 5.5, dpi = cfg$plot_dpi
)

p_cv <- ggplot(
  cv_summary,
  aes(forcats::fct_reorder(model, CV_RMSE_mg_L),
      CV_RMSE_mg_L)
) +
  geom_col() +
  coord_flip() +
  labs(
    title = "Bayesian gamma age models: grouped held-out prediction",
    x = NULL,
    y = "Grouped five-fold CV RMSE (mg/L)"
  )
print(p_cv)
ggsave(
  "plots/bayesian_age_parameter_CV.png",
  p_cv, width = 9, height = 5.5, dpi = cfg$plot_dpi
)

cat("\nBAYESIAN POSTERIOR PARAMETER SUMMARIES:\n")
print(posterior, n = Inf)
cat("\nGAMMA SHAPE VS K POSTERIOR ASSOCIATION:\n")
print(k_alpha_association)
cat("\nSPLIT RHAT (simple convergence check):\n")
print(diagnostics, n = Inf)
cat("\nGROUPED FIVE-FOLD CV:\n")
print(cv_summary, n = Inf)
cat(
  "\nINTERPRETATION:\n",
  "* A gamma shape prior is not an independent age-tracer constraint.\n",
  "* Changes in gamma prior median should be compared with changes in\n",
  "  its posterior; poorly identified alpha leads to prior-sensitive k.\n",
  "* A weakly predictive gamma model is not proof of piston flow.\n",
  "* LSR effective DO is estimated jointly, but a site-specific LSR\n",
  "  oxygen loss process is NOT estimated in this experiment.\n",
  "* No river-recharge DO loss parameter is fitted.\n",
  "* A positive apparent k is an AGE-ASSOCIATED signal, not proof\n",
  "  of oxygen loss solely within the saturated aquifer.\n",
  "* This model assumes independent normal measurement plus residual\n",
  "  errors and one gamma age distribution for both recharge sources.\n",
  sep = ""
)
