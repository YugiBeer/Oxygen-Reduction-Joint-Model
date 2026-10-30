# 14_joint_Bayesian_LSR_DO_grouped_CV.R
#
# Scientific question: does allowing WELL-TO-WELL VARIATION in effective
# LAND-SURFACE-RECHARGE (LSR) DO improve predictions of groundwater DO
# under the FULL joint Bayesian mixing + MRT model?
#
# Compare exactly two models under IDENTICAL priors, training samples,
# group CV folds, river DO prior, mixing tracer and true-MRT likelihoods.
#  A. Common effective LSR DO for every well
#  B. Variable LSR DO with bounded logistic-normal hierarchical prior
#
# Held-out DO is completely EXCLUDED from Stan likelihood (has_DO=0).
# Held-out chloride, delta18O, excess air, reported MRT DO NOT disappear:
# they help infer a NEW WELL's recharge fraction/age and hence predict DO.
# The held-out DO value is replaced by zero in Stan input for clarity;
# it is only used to score predictions AFTER fitting.
#
# xs_N2 is NOT treated as a conservative tracer. use_n2=0 for both
# models to avoid leakage via a reactive observation when predicting DO.
#
# NOTE: This is strictly cross-validated predictive accuracy, not a
# causal identification of where age-associated oxygen loss occurred.
#
# Input: local private data/groundwater.csv. Do NOT commit input or output.
#
# Requires: cmdstanr + CmdStan, posterior, tidyverse, here.
library(tidyverse)
setwd(here::here("Oxygen-Reduction-Joint-Model"))
if (!requireNamespace("cmdstanr", quietly = TRUE) ||
    !requireNamespace("posterior", quietly = TRUE)) {
  stop("Install cmdstanr, CmdStan, posterior; see README.md.")
}

cfg <- list(
  input = "data/groundwater.csv",
  stan_file = "stan/12_joint_bayesian_mixing.stan",
  mean_centered_MRT_log_SD = .25,
  tracer_prior_scale = 1,
  lsr_prior_a = 1,
  lsr_prior_b = 1,
  gamma_age_distribution = FALSE,
  gamma_prior_median = 10,
  seed = 20261030L,
  folds = 5L,
  chains = 4L,
  parallel_chains = 4L,
  warmup = 1000L,
  sampling = 1000L,
  adapt_delta = .97,
  max_treedepth = 13L,
  expected_n = 41L,
  max_rhat = 1.05,
  min_bulk_ESS = 400,
  store_fit_objects = FALSE,
  run_only_fold = NULL  # e.g., 1L for debugging; NULL runs all 5
)

if (!file.exists(cfg$input)) {
  stop("Missing private data/groundwater.csv; see README.")
}
if (!file.exists(cfg$stan_file)) {
  stop("Missing Stan file ", cfg$stan_file)
}
dir.create("outputs", showWarnings = FALSE, recursive = TRUE)
dir.create("plots", showWarnings = FALSE, recursive = TRUE)

raw <- read_csv(
  cfg$input, show_col_types = FALSE,
  col_types = cols(Site_ID = col_character(), .default = col_guess())
)
cols <- c(
  "Site_ID", "DO", "MRT", "O18", "Cl", "xs_Air",
  "xs_Air_uncertainty", "xs_N2", "xs_N2_sig"
)
stopifnot(all(cols %in% colnames(raw)))

d <- raw %>%
  filter(is.finite(DO), DO >= 0,
         is.finite(MRT), MRT >= 0) %>%
  mutate(
    sample_index = row_number(),
    # Hold out wells from related named groups together.
    cv_group = case_when(
      str_detect(Site_ID, "Eastbourne") ~ "Eastbourne",
      str_detect(Site_ID, "Gilbertson") ~ "Gilbertson",
      str_detect(Site_ID, "Flaxmere") ~ "Flaxmere",
      TRUE ~ Site_ID
    )
  )
stopifnot(nrow(d) == cfg$expected_n)

as_obs <- function(x) ifelse(is.finite(x), as.numeric(x), 0)
as_sigma <- function(x, default) {
  ifelse(is.finite(x) & x > 0, as.numeric(x), default)
}

stan_base <- list(
  N = nrow(d),
  has_o18 = as.integer(is.finite(d$O18)),
  has_cl = as.integer(is.finite(d$Cl)),
  has_air = as.integer(is.finite(d$xs_Air)),
  has_n2 = as.integer(is.finite(d$xs_N2)),
  has_DO = rep(1L, nrow(d)),
  o18 = as_obs(d$O18),
  cl = as_obs(d$Cl),
  air = as_obs(d$xs_Air),
  do_obs = as.numeric(d$DO),
  n2_obs = as_obs(d$xs_N2),
  cl_sd = pmax(1, .03 * as_obs(d$Cl)),
  air_sd = as_sigma(d$xs_Air_uncertainty, .15),
  n2_sd = as_sigma(d$xs_N2_sig, .3),
  mrt_reported = as.numeric(d$MRT),
  log_mrt_SD = cfg$mean_centered_MRT_log_SD,
  use_n2 = 0L,
  use_gamma = as.integer(cfg$gamma_age_distribution),
  gamma_prior_median = cfg$gamma_prior_median,
  tracer_prior_scale = cfg$tracer_prior_scale,
  lsr_prior_a = cfg$lsr_prior_a,
  lsr_prior_b = cfg$lsr_prior_b,
  vary_lsr_do = 0L
)

# Create reproducible, near-balanced FOLDS OF GROUPS rather than of rows.
set.seed(cfg$seed)
g <- split(seq_len(nrow(d)), d$cv_group)
g_order <- sample(seq_along(g))
g_order <- g_order[order(-lengths(g)[g_order])]
folds <- vector("list", cfg$folds)
for (i in g_order) {
  j <- which.min(lengths(folds))
  folds[[j]] <- c(folds[[j]], g[[i]])
}
stopifnot(
  all(lengths(folds) > 0L),
  identical(sort(unlist(folds)), seq_len(nrow(d)))
)
fold_map <- map_dfr(seq_along(folds), function(j) {
  tibble(
    fold = j,
    sample_index = folds[[j]],
    Site_ID = d$Site_ID[folds[[j]]],
    Depth = d$Depth[folds[[j]]]
  )
})
write_csv(fold_map, "outputs/joint_LSR_DO_CV_folds.csv")

selected_folds <- if (is.null(cfg$run_only_fold)) {
  seq_along(folds)
} else {
  stopifnot(cfg$run_only_fold %in% seq_along(folds))
  cfg$run_only_fold
}

model <- cmdstanr::cmdstan_model(cfg$stan_file, quiet = TRUE)
out <- list()
diag <- list()
scenario_names <- c(
  "Common effective LSR DO",
  "Variable effective LSR DO"
)

# Five folds x two models. All 41 tracer observations are present in
# BOTH models. The held-out DO observations are completely masked.
for (j in selected_folds) {
  held <- folds[[j]]
  for (varies in 0:1) {
    model_name <- scenario_names[varies + 1L]
    message(
      "\nFold ", j, "/", cfg$folds, " | ", model_name,
      " | held-out DO n=", length(held)
    )
    x <- stan_base
    x$has_DO[held] <- 0L
    x$do_obs[held] <- 0  # held-out DO not supplied to sampler
    x$vary_lsr_do <- as.integer(varies)

    fit <- model$sample(
      data = x,
      seed = cfg$seed + 1000L*j + 10000L*varies,
      chains = cfg$chains,
      parallel_chains = cfg$parallel_chains,
      iter_warmup = cfg$warmup,
      iter_sampling = cfg$sampling,
      adapt_delta = cfg$adapt_delta,
      max_treedepth = cfg$max_treedepth,
      refresh = 250
    )
    vars <- paste0("mean_DO[", held, "]")
    draws <- posterior::as_draws_matrix(
      fit$draws(variables = c(vars, "sigma_DO"))
    )
    # Use the expectation of DO, averaged over latent site-specific
    # recharge fractions, age and source DO, for squared-error scoring.
    point_predictions <- colMeans(draws[, vars, drop = FALSE])

    # Full posterior predictive intervals also include unmodelled DO
    # residual error. These are conditional on the observed tracers/MRT,
    # NOT on the held-out groundwater DO.
    set.seed(cfg$seed + j + 23L*varies)
    s <- draws[, "sigma_DO"]
    pp <- vapply(vars, function(v) {
      draws[, v] + rnorm(length(s), sd = sqrt(s^2 + .2^2))
    }, numeric(length(s)))
    if (is.null(dim(pp))) pp <- matrix(pp, ncol = 1L)

    out[[length(out) + 1L]] <- tibble(
      model = model_name,
      fold = j,
      sample_index = held,
      Site_ID = d$Site_ID[held],
      Depth = d$Depth[held],
      MRT = d$MRT[held],
      observed_DO = d$DO[held],
      predicted_DO = as.numeric(point_predictions),
      predicted_lower95 = apply(pp, 2, quantile, probs = .025),
      predicted_upper95 = apply(pp, 2, quantile, probs = .975),
      # Predictive 95% coverage, not coverage of posterior mean.
      covered95 = observed_DO >= predicted_lower95 &
                  observed_DO <= predicted_upper95
    )

    sumpars <- fit$summary(
      variables = c(
        "k", "central_LSR_DO_loss", "tau_LSR",
        "sigma_DO", "river_DO"
      )
    )
    diag[[length(diag) + 1L]] <- sumpars %>%
      mutate(
        model = model_name,
        fold = j,
        .before = 1
      )
    diagnostic <- fit$diagnostic_summary()
    if (any(sumpars$rhat > cfg$max_rhat, na.rm = TRUE) ||
        any(sumpars$ess_bulk < cfg$min_bulk_ESS, na.rm = TRUE) ||
        any(diagnostic$num_divergent > 0)) {
      warning(
        "MCMC diagnostics require review: ",
        "fold ", j, ", ", model_name
      )
    }
    if (cfg$store_fit_objects) {
      fit$save_object(
        file = file.path(
          "outputs",
          paste0(
            "joint_LSR_CV_fold", j, "_vary", varies, ".rds"
          )
        )
      )
    }
    # Checkpoint after every model.
    write_csv(
      bind_rows(out),
      "outputs/joint_LSR_DO_CV_predictions.csv"
    )
    write_csv(
      bind_rows(diag),
      "outputs/joint_LSR_DO_CV_convergence.csv"
    )
  }
}

preds <- bind_rows(out)
cv <- preds %>%
  group_by(model) %>%
  summarise(
    folds_completed = n_distinct(fold),
    n = n(),
    RMSE = sqrt(mean((observed_DO - predicted_DO)^2)),
    MAE = mean(abs(observed_DO - predicted_DO)),
    bias = mean(predicted_DO - observed_DO),
    coverage95 = mean(covered95),
    .groups = "drop"
  )
write_csv(cv, "outputs/joint_LSR_DO_CV_summary.csv")
write_csv(
  bind_rows(diag),
  "outputs/joint_LSR_DO_CV_convergence.csv"
)

# A paired per-well diagnostic using the IDENTICAL held-out observations.
paired <- preds %>%
  select(sample_index, Site_ID, fold, model, observed_DO, predicted_DO) %>%
  pivot_wider(
    names_from = model,
    values_from = predicted_DO
  ) %>%
  mutate(
    abs_error_common = abs(
      observed_DO - .data[["Common effective LSR DO"]]
    ),
    abs_error_variable = abs(
      observed_DO - .data[["Variable effective LSR DO"]]
    ),
    variable_minus_common_abs_error = (
      abs_error_variable - abs_error_common
    )
  )
write_csv(paired, "outputs/joint_LSR_DO_CV_paired_errors.csv")

p <- preds %>%
  ggplot(aes(observed_DO, predicted_DO)) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", colour = "grey60") +
  geom_point(size = 2.5) +
  facet_wrap(~ model) +
  coord_equal() +
  labs(
    title = "Joint Bayesian mixing: held-out groundwater DO",
    subtitle = "Held-out DO masked; conservative tracers and reported MRT retained",
    x = "Observed DO (mg/L)",
    y = "Bayesian predicted DO (mg/L)"
  ) +
  theme_bw(base_size = 12)
print(p)
ggsave(
  "plots/joint_LSR_DO_grouped_CV_predictions.png",
  p, width = 10, height = 5.5, dpi = 300
)

cat("\nGROUPED FIVE-FOLD CV SUMMARY:\n")
print(cv, n = Inf)
cat("\nPAIRWISE ERRORS (variable minus common; negative = better):\n")
print(summary(paired$variable_minus_common_abs_error))
cat(
  "\nINTERPRETATION:\n",
  "1. Held-out DO values are masked; their tracers and MRT are retained.\n",
  "2. Compare BOTH RMSE and calibration/convergence diagnostics.\n",
  "3. Similar RMSE indicates the data do not justify more complex LSR DO.\n",
  "4. This tests prediction, not exact LSR vs saturated aquifer location.\n",
  "5. Excess N2 is excluded: it is not a conservative mixing tracer.\n",
  "6. The endmember prior centres derive partly from the same dataset;\n",
  "   strict external validation needs independent endmember data.\n",
  "7. Some sample groups share coordinates; this is grouped by well\n",
  "   labels, not by spatially independent recharge domains.\n",
  "8. Saved tables and plot are LOCAL and gitignored.\n",
  sep = ""
)
