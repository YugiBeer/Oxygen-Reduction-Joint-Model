# 08_LSR_vs_aquifer_oxygen_loss.R
#
# Main scientific test: can LAND-SURFACE RECHARGE (LSR) oxygen loss,
# including between-well variability, improve on an apparent first-order
# groundwater-age DO depletion model?
#
# No riverbed / river-recharge oxygen-loss parameter is fitted here.
# Monitored river DO: prior mean 7.5 mg/L, uncertainty of mean 0.7 mg/L
# (the latter assumed, not measured). The river daily-mean variability
# of about 1.9 mg/L is age-filtered with an ASSUMED exponential travel
# time distribution. Its variance partition is a sensitivity assumption.
#
# Model A: effective LSR DO fixed at 9.5 mg/L; k estimated.
# Model B: common effective LSR DO in [0,9.5]; k estimated.
# Model C: effective LSR DO varies across wells within [0,9.5]
#          using a logit-normal hierarchical distribution; k estimated.
# Model D: same variable-LSR model as C, but k=0.
#
# k is an AGE-ASSOCIATED coefficient; DO+MRT data alone do not prove
# oxygen consumption happened inside the saturated regional aquifer.
#
# This script uses base R Metropolis-Hastings and tidyverse only.
# It does not publish or modify private groundwater input data.
library(tidyverse)
setwd(here::here("Oxygen-Reduction-Joint-Model"))

# ================================================================
# USER SETTINGS
# ================================================================
cfg <- list(
  river_mean_prior = 7.5,          # DO mg/L (river record mean)
  river_mean_prior_SD = 0.7,       # assumed uncertainty of mean (mg/L)
  river_daily_mean_SD = 1.9,       # observed temporal variation of DAILY means
  persistent_river_variance = 0.25,# fraction assumed to survive all ages
  seasonal_period_years = 1,
  assumed_diurnal_input_SD = 1,    # mg/L, NOT measured from daily means
  year_days = 365.25,
  LSR_initial_DO_reference = 9.5,  # mg/L; near-saturation reference
  LSR_spatial_logit_SD_prior = 1.5,# half-normal prior on logit SD tau
  residual_SD_prior = 1.5,         # half-normal scale, mg/L
  measurement_SD = 0.2,           # mg/L; assumed
  k_prior_SD = 0.10,              # year^-1; half-normal prior
  LSR_quadrature_nodes = 21L,      # midpoint-normal integration nodes
  n_chains = 4L,
  warmup = 1100L,
  kept_per_chain = 1100L,
  thin = 2L,
  cv_n_chains = 2L,
  cv_warmup = 550L,
  cv_kept_per_chain = 500L,
  cv_folds = 5L,
  exclude_2801_sensitivity = TRUE,
  seed = 20261030L,
  plot_font_size = 12,
  plot_dpi = 300
)

# ================================================================
# DATA: identical 29-well cohort as script 01
# ================================================================
dir.create("outputs", showWarnings = FALSE, recursive = TRUE)
dir.create("plots", showWarnings = FALSE, recursive = TRUE)
input <- "outputs/three_tracer_mixing.csv"
if (!file.exists(input)) stop("Run scripts/01_three_tracer_mixing.R first.")

d <- read_csv(input, show_col_types = FALSE,
              col_types = cols(Site_ID = col_character(),
                               .default = col_guess())) %>%
  filter(mixing_reliable,
         is.finite(DO), DO > 0,
         is.finite(MRT), MRT >= 0,
         is.finite(rainfall_contribution),
         is.finite(river_contribution)) %>%
  transmute(
    Site_ID, DO, MRT,
    fR = river_contribution,
    fLSR = rainfall_contribution
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
  all(abs(d$fR + d$fLSR - 1) < 1e-8)
)

models <- c(
  A = "A - Age loss only",
  B = "B - Age + common LSR loss",
  C = "C - Age + varying LSR loss",
  D = "D - Varying LSR loss only (k=0)"
)
nodes <- qnorm((seq_len(cfg$LSR_quadrature_nodes) - .5) /
                 cfg$LSR_quadrature_nodes)

# ================================================================
# PRIORS AND LIKELIHOOD
# ================================================================
# River source DO: R_i ~ N(mu_R, v_R(MRT_i)), with mu_R
#   having its own N(7.5, 0.7) prior.
# Effective LSR DO:
#   A: P_i = 9.5;
#   B: P_i = P, P ~ Uniform(0,9.5);
#   C/D: P_i = 9.5 * logistic(eta + tau*Z_i),
#        Z_i ~ N(0,1), tau ~ HalfNormal(1.5);
#        central P=9.5*logistic(eta) ~ Uniform(0,9.5).
# Measurement DO_i = [fR_i R_i + fLSR_i P_i] exp(-k*MRT_i)
#                    + epsilon_i
# epsilon_i ~ N(0, sigma_extra^2 + measurement_SD^2).
#
# For C/D the well-specific, physically bounded LSR distribution is
# integrated out by fixed normal-quantile quadrature. No well gets its
# own freely fitted LSR loss.
river_variance <- function(t) {
  h_year2 <- 1 / (1 + (2 * pi * t /
    cfg$seasonal_period_years)^2)
  h_day2 <- 1 / (1 + (2 * pi * t * cfg$year_days)^2)
  cfg$river_daily_mean_SD^2 * (
    cfg$persistent_river_variance +
      (1 - cfg$persistent_river_variance) * h_year2
  ) + cfg$assumed_diurnal_input_SD^2 * h_day2
}

decode <- function(z, model) {
  j <- 1L
  R <- z[j]; j <- j + 1L
  if (model == "A") {
    Pcentral <- cfg$LSR_initial_DO_reference
    eta <- NA_real_
  } else {
    eta <- z[j]; j <- j + 1L
    Pcentral <- cfg$LSR_initial_DO_reference * plogis(eta)
  }
  if (model == "D") {
    k <- 0
  } else {
    k <- exp(z[j]); j <- j + 1L
  }
  sigma_extra <- exp(z[j]); j <- j + 1L
  tau <- if (model %in% c("C", "D")) exp(z[j]) else 0
  list(R = R, Pcentral = Pcentral, eta = eta,
       k = k, sigma_extra = sigma_extra, tau = tau)
}

LSR_nodes <- function(p, model) {
  if (model %in% c("A", "B")) return(p$Pcentral)
  cfg$LSR_initial_DO_reference * plogis(p$eta + p$tau * nodes)
}

# Per-observation log marginal likelihood. Integrate both the river DO
# random effect (analytically) and the LSR DO random effect (quadrature).
log_density <- function(p, data, model) {
  P <- LSR_nodes(p, model)
  a <- exp(-p$k * data$MRT)
  river_part <- data$fR * p$R * a
  lsr_part <- data$fLSR * a
  mu <- outer(lsr_part, P, "*") +
    matrix(river_part, nrow = nrow(data), ncol = length(P))
  variance <- data$fR^2 * river_variance(data$MRT) * a^2 +
    p$sigma_extra^2 + cfg$measurement_SD^2
  logdens <- -.5 * (
    log(2 * pi * matrix(variance, nrow = nrow(data),
                      ncol = length(P))) +
      (matrix(data$DO, nrow = nrow(data), ncol = length(P)) -
         mu)^2 /
      matrix(variance, nrow = nrow(data), ncol = length(P))
  )
  apply(logdens, 1L, function(x) {
    mx <- max(x)
    mx + log(mean(exp(x - mx)))
  })
}

log_posterior <- function(z, data, model) {
  p <- decode(z, model)
  if (any(!is.finite(z)) ||
      !is.finite(p$R + p$Pcentral + p$k +
                   p$sigma_extra + p$tau) ||
      p$R <= 0 || p$R >= 12 ||
      p$k >= 3 || p$sigma_extra >= 15 || p$tau >= 30) {
    return(-Inf)
  }

  logprior <- dnorm(
    p$R, cfg$river_mean_prior, cfg$river_mean_prior_SD,
    log = TRUE
  )
  # Uniform prior on central effective LSR DO with a logit transform.
  if (model != "A") {
    if (p$Pcentral <= 0 ||
        p$Pcentral >= cfg$LSR_initial_DO_reference) return(-Inf)
    logprior <- logprior +
      log(p$Pcentral) +
      log(cfg$LSR_initial_DO_reference - p$Pcentral) -
      log(cfg$LSR_initial_DO_reference)
  }
  # Half-normal priors in physical space; add transformation Jacobians.
  if (model != "D") {
    logprior <- logprior +
      dnorm(p$k, 0, cfg$k_prior_SD, log = TRUE) +
      log(p$k)
  }
  logprior <- logprior +
    dnorm(p$sigma_extra, 0, cfg$residual_SD_prior,
          log = TRUE) + log(p$sigma_extra)

  if (model %in% c("C", "D")) {
    logprior <- logprior +
      dnorm(p$tau, 0, cfg$LSR_spatial_logit_SD_prior,
            log = TRUE) + log(p$tau)
  }
  logprior + sum(log_density(p, data, model))
}

starting_state <- function(model, chain) {
  z <- cfg$river_mean_prior + (chain - 1) * .1
  if (model != "A") {
    start_P <- c(8.5, 6, 4, 9)[(chain - 1L) %% 4L + 1L]
    z <- c(z, qlogis(start_P / cfg$LSR_initial_DO_reference))
  }
  if (model != "D") {
    start_k <- c(.04, .07, .09, .12)[(chain - 1L) %% 4L + 1L]
    z <- c(z, log(start_k))
  }
  z <- c(z, log(c(.6, 1, 1.4, 1.8)[
    (chain - 1L) %% 4L + 1L]))
  if (model %in% c("C", "D")) {
    start_tau <- c(.2, .7, 1.5, 3)[(chain - 1L) %% 4L + 1L]
    z <- c(z, log(start_tau))
  }
  z
}

derive <- function(z, model) {
  p <- decode(z, model)
  Pnodes <- LSR_nodes(p, model)
  c(
    river_mean_DO = p$R,
    LSR_mean_entry_DO = mean(Pnodes),
    LSR_mean_loss = cfg$LSR_initial_DO_reference - mean(Pnodes),
    LSR_between_well_SD = sqrt(mean(
      (Pnodes - mean(Pnodes))^2
    )),
    apparent_k = p$k,
    residual_SD = p$sigma_extra,
    LSR_logit_SD = p$tau
  )
}

# Metropolis component updates; proposal adaptation occurs in warmup only.
sample_chain <- function(data, model, chain, seed, warmup,
                         n_draws, thin) {
  set.seed(seed)
  z <- starting_state(model, chain)
  step <- c(.22)
  if (model != "A") step <- c(step, .45)
  if (model != "D") step <- c(step, .22)
  step <- c(step, .20)
  if (model %in% c("C", "D")) step <- c(step, .20)
  stopifnot(length(step) == length(z))
  lp <- log_posterior(z, data, model)
  if (!is.finite(lp)) stop("Invalid initial state in model ", model)
  accept <- integer(length(z))

  keep <- matrix(NA_real_, nrow = n_draws, ncol = 7)
  colnames(keep) <- c(
    "river_mean_DO", "LSR_mean_entry_DO", "LSR_mean_loss",
    "LSR_between_well_SD", "apparent_k", "residual_SD", "LSR_logit_SD"
  )
  keep_z <- vector("list", n_draws)
  out <- 0L
  for (it in seq_len(warmup + n_draws * thin)) {
    for (j in sample.int(length(z))) {
      proposed <- z
      proposed[j] <- proposed[j] + rnorm(1, sd = step[j])
      cand_lp <- log_posterior(proposed, data, model)
      if (is.finite(cand_lp) && log(runif(1)) < cand_lp - lp) {
        z <- proposed
        lp <- cand_lp
        accept[j] <- accept[j] + 1L
      }
    }
    if (it <= warmup && it %% 100L == 0L) {
      step <- pmax(.015, pmin(2,
        step * exp(pmax(-.25, pmin(.25, accept / 100 - .34)))
      ))
      accept[] <- 0L
    }
    if (it > warmup && (it - warmup) %% thin == 0L) {
      out <- out + 1L
      keep[out, ] <- derive(z, model)
      keep_z[[out]] <- z
    }
  }
  list(
    summary = as_tibble(keep) %>% mutate(chain = chain),
    z = keep_z
  )
}

fit_chains <- function(data, model, n_chains, warmup,
                       n_draws, seed) {
  chains <- map(seq_len(n_chains), function(ch) {
    sample_chain(
      data, model, ch, seed + ch * 1103L,
      warmup, n_draws, cfg$thin
    )
  })
  list(
    summary = bind_rows(map(chains, "summary")),
    z = unlist(map(chains, "z"), recursive = FALSE)
  )
}

summarise_draws <- function(draws, model, selection) {
  draws %>%
    select(-chain) %>%
    pivot_longer(everything(), names_to = "parameter") %>%
    group_by(parameter) %>%
    summarise(
      median = median(value),
      lower95 = unname(quantile(value, .025)),
      upper95 = unname(quantile(value, .975)),
      .groups = "drop"
    ) %>%
    mutate(model = models[[model]],
           sample = selection, .before = 1)
}

basic_Rhat <- function(draws, field) {
  cc <- split(draws[[field]], draws$chain)
  n <- min(lengths(cc))
  if (n < 20L) return(NA_real_)
  v <- mean(vapply(cc, function(x) var(x[seq_len(n)]),
                   numeric(1)))
  if (!is.finite(v) || v < 1e-12) return(NA_real_)
  B <- n * var(vapply(cc, function(x) mean(x[seq_len(n)]),
                      numeric(1)))
  sqrt(((n - 1) / n * v + B / n) / v)
}

# ================================================================
# REPRODUCIBLE GROUPED 5-FOLD VALIDATION
# ================================================================
set.seed(cfg$seed)
groups <- split(seq_len(nrow(d)), d$cv_group)
ord <- sample(seq_along(groups))
ord <- ord[order(-lengths(groups)[ord])]
folds <- vector("list", cfg$cv_folds)
for (g in ord) {
  smallest <- which.min(lengths(folds))
  folds[[smallest]] <- c(folds[[smallest]], groups[[g]])
}
stopifnot(
  all(lengths(folds) > 0),
  identical(sort(unlist(folds)), seq_len(nrow(d)))
)

logmeanexp <- function(x) {
  mm <- max(x)
  mm + log(mean(exp(x - mm)))
}

predict_held_out <- function(fit, held, model) {
  z <- fit$z
  n <- length(z)
  all_means <- matrix(NA_real_, nrow = n, ncol = nrow(held))
  all_lpd <- matrix(NA_real_, nrow = n, ncol = nrow(held))

  for (i in seq_len(n)) {
    p <- decode(z[[i]], model)
    P <- LSR_nodes(p, model)
    a <- exp(-p$k * held$MRT)
    all_means[i, ] <- (
      held$fR * p$R + held$fLSR * mean(P)
    ) * a
    all_lpd[i, ] <- log_density(p, held, model)
  }
  tibble(
    Site_ID = held$Site_ID,
    MRT = held$MRT,
    LSR_fraction = held$fLSR,
    observed_DO = held$DO,
    predicted_DO = colMeans(all_means),
    log_predictive_density = apply(all_lpd, 2, logmeanexp)
  )
}

run_analysis <- function(dd, selection) {
  cat("\nLSR vs aquifer models, selection:", selection,
      "n =", nrow(dd), "\n")
  full <- vector("list", length(models))
  names(full) <- names(models)
  summaries <- list()
  diagnostics <- list()
  predictions <- list()

  for (i in seq_along(models)) {
    model <- names(models)[i]
    cat("Full posterior:", models[[i]], "\n")
    fit <- fit_chains(
      dd, model, cfg$n_chains, cfg$warmup,
      cfg$kept_per_chain, cfg$seed + 50000L * i
    )
    full[[i]] <- fit
    summaries[[i]] <- summarise_draws(fit$summary, model, selection)
    diagnostics[[i]] <- tibble(
      sample = selection, model = models[[model]],
      k_Rhat_basic = basic_Rhat(fit$summary, "apparent_k"),
      LSR_SD_Rhat_basic = basic_Rhat(fit$summary,
                                    "LSR_between_well_SD")
    )

    # CV held-out well samples never enter their training posterior.
    predictions[[i]] <- map_dfr(seq_along(folds), function(j) {
      held_idx <- intersect(
        folds[[j]], which(d$Site_ID %in% dd$Site_ID)
      )
      if (length(held_idx) == 0) return(tibble())
      held <- d[held_idx, , drop = FALSE]
      train <- dd %>% filter(!Site_ID %in% held$Site_ID)
      cat("  CV:", model, "fold", j, "ntrain", nrow(train), "\n")
      trainfit <- fit_chains(
        train, model,
        cfg$cv_n_chains, cfg$cv_warmup, cfg$cv_kept_per_chain,
        cfg$seed + 200000L + i * 10000L + j * 101L
      )
      predict_held_out(trainfit, held, model) %>%
        mutate(
          sample = selection, model = models[[model]],
          fold = j, .before = 1
        )
    })
  }

  list(
    full = full,
    posterior = bind_rows(summaries),
    diagnostics = bind_rows(diagnostics),
    CV_predictions = bind_rows(predictions)
  )
}

all_results <- run_analysis(d, "All 29")
if (isTRUE(cfg$exclude_2801_sensitivity)) {
  no_2801 <- run_analysis(
    d %>% filter(Site_ID != "2801"),
    "Excluding 2801"
  )
  results <- list(all_results, no_2801)
} else {
  results <- list(all_results)
}
posterior <- bind_rows(map(results, "posterior"))
diagnostics <- bind_rows(map(results, "diagnostics"))
cv_predictions <- bind_rows(map(results, "CV_predictions"))

cv_summary <- cv_predictions %>%
  group_by(sample, model) %>%
  summarise(
    n = n(),
    CV_RMSE_mg_L = sqrt(mean(
      (observed_DO - predicted_DO)^2
    )),
    CV_MAE_mg_L = mean(abs(observed_DO - predicted_DO)),
    summed_log_predictive_density = sum(log_predictive_density),
    .groups = "drop"
  )

write_csv(posterior, "outputs/LSR_loss_posterior_summary.csv")
write_csv(diagnostics, "outputs/LSR_loss_convergence.csv")
write_csv(cv_predictions, "outputs/LSR_loss_CV_predictions.csv")
write_csv(cv_summary, "outputs/LSR_loss_CV_summary.csv")

if (any(diagnostics$k_Rhat_basic > 1.05, na.rm = TRUE) ||
    any(diagnostics$LSR_SD_Rhat_basic > 1.05, na.rm = TRUE)) {
  warning("Some chains may not have converged; increase warmup/draws.")
}

# ================================================================
# PLOTS
# ================================================================
theme_set(theme_bw(base_size = cfg$plot_font_size))
full_draws <- map_dfr(
  names(models),
  function(m) all_results$full[[m]]$summary %>%
    mutate(model = models[[m]])
)
p_k <- full_draws %>%
  filter(model != models[["D"]]) %>%
  ggplot(aes(apparent_k, colour = model)) +
  geom_density(linewidth = .9) +
  labs(
    title = "Apparent first-order DO depletion: LSR-loss uncertainty",
    subtitle = "Age-associated coefficient, not a measured aquifer reaction rate",
    x = expression(k~(year^{-1})),
    y = "Posterior density", colour = "Model"
  )
print(p_k)
ggsave("plots/LSR_loss_posterior_k.png", p_k,
       width = 10, height = 6, dpi = cfg$plot_dpi)

p_predictions <- cv_predictions %>%
  filter(sample == "All 29") %>%
  ggplot(aes(observed_DO, predicted_DO)) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", colour = "grey50") +
  geom_point(aes(colour = LSR_fraction), size = 2.8) +
  facet_wrap(~ model, ncol = 2) +
  coord_equal() +
  scale_colour_viridis_c(labels = scales::percent) +
  labs(
    title = "Held-out DO predictions: LSR versus apparent aquifer loss",
    x = "Measured groundwater DO (mg/L)",
    y = "Predicted groundwater DO (mg/L)",
    colour = "LSR fraction"
  )
print(p_predictions)
ggsave("plots/LSR_loss_CV_predictions.png", p_predictions,
       width = 11, height = 9, dpi = cfg$plot_dpi)

p_resid <- cv_predictions %>%
  filter(sample == "All 29") %>%
  mutate(residual = observed_DO - predicted_DO) %>%
  ggplot(aes(LSR_fraction, residual)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             colour = "grey50") +
  geom_point(aes(colour = MRT), size = 2.8) +
  facet_wrap(~ model, ncol = 2) +
  scale_x_continuous(labels = scales::percent) +
  labs(
    title = "Held-out residuals versus LSR contribution",
    x = "LSR fraction", y = "Measured - predicted DO (mg/L)",
    colour = "MRT (years)"
  )
print(p_resid)
ggsave("plots/LSR_loss_CV_residuals.png", p_resid,
       width = 11, height = 9, dpi = cfg$plot_dpi)

cat("\nPosterior summaries:\n")
print(posterior)
cat("\nGrouped CV performance:\n")
print(cv_summary)
cat("\nMCMC convergence (basic Rhat):\n")
print(diagnostics)

cat(
  "\nINTERPRETATION:\n",
  "* LSR loss is defined ONLY for land-surface recharge.\n",
  "* River DO is referenced to monitoring, with age-filtered uncertainty;\n",
  "  no river-recharge DO loss parameter is estimated.\n",
  "* Model C/D LSR entry DO is physically bounded in [0,9.5].\n",
  "* Model D tests zero age-dependent depletion; it is not proof that\n",
  "  recharge pathways can be uniquely distinguished from aquifer transport.\n",
  "* The 1.9 mg/L river SD is DAILY-MEAN variability, not diurnal SD.\n",
  "* The 25% persistent variance and 1 mg/L diurnal SD are assumptions.\n",
  "* Mixing fractions and MRT are treated as known; only five wells are\n",
  "  LSR-dominated, so LSR variability may be prior-sensitive.\n",
  "* The role of aquifer transport cannot be distinguished from\n",
  "  recharge characteristics that systematically covary with MRT.\n",
  sep = ""
)
