# 15_tracer_only_mixing_CV.R
#
# Compare the SOURCE MIXING models using only conservative tracer
# observations. DO, MRT, excess N2 and nitrate are not used.
#
# A: two recharge endmembers and common compositions.
# B: A + variable LSR delta18O signature (hierarchical SD).
# C: B + positive, non-source-specific chloride enrichment.
#
# Bayesian prediction test: each sample has one tracer withheld at
# a time, rotating over Cl, delta18O and excess air in 3 experiments.
# Each fit retains the OTHER TWO measured tracers at that sample.
# Related sample groups get the same withheld tracer in each run.
# Same observation mask for every candidate model. NO groundwater
# DO participates in the fit. Source-fraction posteriors are inferred
# only from observed tracers.
#
# Requires CmdStanR, CmdStan, posterior, tidyverse and here.
# Local private data only; outputs/ and plots/ must be .gitignored.
library(tidyverse)
setwd(here::here("Oxygen-Reduction-Joint-Model"))

cfg <- list(
  input = "data/groundwater.csv",
  stan_file = "stan/15_tracer_only_mixing_models.stan",
  endmember_prior_scale = 1,
  variants = 1:3,                   # 1 common; 2 variable O18; 3 +Cl
  rotations = 1:3,                  # together cover all 3 tracers
  seed = 20261031L,
  chains = 4L,
  parallel_chains = 4L,
  warmup = 1000L,
  sampling = 1000L,
  adapt_delta = .97,
  max_treedepth = 13L,
  fit_full_data = TRUE,             # Save all-data fLSR intervals
  minimum_bulk_ESS = 400,
  maximum_Rhat = 1.05
)

if (!requireNamespace("cmdstanr", quietly = TRUE) ||
    !requireNamespace("posterior", quietly = TRUE)) {
  stop("Requires cmdstanr, CmdStan and posterior; see README.")
}
if (!file.exists(cfg$input) || !file.exists(cfg$stan_file)) {
  stop("Private data/groundwater.csv or Stan model missing.")
}
dir.create("outputs", showWarnings = FALSE, recursive = TRUE)
dir.create("plots", showWarnings = FALSE, recursive = TRUE)

d <- read_csv(
  cfg$input, show_col_types = FALSE,
  col_types = cols(Site_ID = col_character(), .default = col_guess())
)
req <- c(
  "Site_ID", "Depth", "O18", "Cl", "xs_Air",
  "xs_Air_uncertainty"
)
stopifnot(all(req %in% names(d)), nrow(d) >= 40L)
d <- d %>% mutate(
  sample_index = row_number(),
  cv_group = case_when(
    str_detect(Site_ID, "Eastbourne") ~ "Eastbourne",
    str_detect(Site_ID, "Gilbertson") ~ "Gilbertson",
    str_detect(Site_ID, "Flaxmere") ~ "Flaxmere",
    TRUE ~ Site_ID
  )
)

numeric_or_zero <- function(x) {
  ifelse(is.finite(x), as.numeric(x), 0)
}
positive_or_default <- function(x, def) {
  ifelse(is.finite(x) & x > 0, as.numeric(x), def)
}

present <- cbind(
  O18 = is.finite(d$O18),
  Cl = is.finite(d$Cl),
  excess_air = is.finite(d$xs_Air)
)
stopifnot(all(rowSums(present) >= 2L))
stan_base <- list(
  N = nrow(d),
  O18 = numeric_or_zero(d$O18),
  Cl = numeric_or_zero(d$Cl),
  air = numeric_or_zero(d$xs_Air),
  Cl_sd = pmax(1, .03 * numeric_or_zero(d$Cl)),
  air_sd = positive_or_default(d$xs_Air_uncertainty, .15),
  train_O18 = as.integer(present[, 1]),
  train_Cl = as.integer(present[, 2]),
  train_air = as.integer(present[, 3]),
  variant = 1L,
  endmember_prior_scale = cfg$endmember_prior_scale
)

# Three rotated masks, balanced by SAMPLE GROUP. In each rotation, every
# well is asked to predict ONE withheld tracer, from its other two
# measured tracers. Over the 3 rotations, every available observation
# is withheld once. Missing tracers are never used as test observations.
groups <- unique(d$cv_group)
set.seed(cfg$seed)
shuffled <- sample(groups)
group_assignment <- setNames(
  (seq_along(shuffled) - 1L) %% 3L, shuffled
)
base_tracer <- unname(group_assignment[d$cv_group])
holdout_map <- map_dfr(cfg$rotations, function(rot) {
  which_tracer <- ((base_tracer + rot - 1L) %% 3L) + 1L
  tibble(
    rotation = rot,
    sample_index = d$sample_index,
    Site_ID = d$Site_ID,
    Depth = d$Depth,
    withheld_tracer = c(
      "O18", "Cl", "excess_air"
    )[which_tracer],
    observed = present[cbind(
      seq_len(nrow(d)), which_tracer
    )]
  )
})
write_csv(holdout_map, "outputs/tracer_only_CV_masks.csv")
stopifnot(
  all(holdout_map %>%
        filter(observed) %>%
        count(sample_index, withheld_tracer) %>%
        pull(n) == 1L),
  nrow(holdout_map %>% filter(observed)) == sum(present)
)

model <- cmdstanr::cmdstan_model(
  cfg$stan_file, quiet = TRUE
)
labels <- c(
  "A: fixed endmembers",
  "B: variable LSR O18",
  "C: variable O18 + Cl enrichment"
)

CV <- list()
diagnostics <- list()
fit_number <- 0L

for (variant in cfg$variants) {
  for (rot in cfg$rotations) {
    fit_number <- fit_number + 1L
    held <- holdout_map %>%
      filter(rotation == rot, observed)
    X <- stan_base
    X$variant <- as.integer(variant)
    # Mask withheld measurements from the model likelihood AND from
    # the corresponding data values passed to Stan.
    h_o <- held$sample_index[held$withheld_tracer == "O18"]
    h_c <- held$sample_index[held$withheld_tracer == "Cl"]
    h_a <- held$sample_index[held$withheld_tracer == "excess_air"]
    X$train_O18[h_o] <- 0L
    X$train_Cl[h_c] <- 0L
    X$train_air[h_a] <- 0L
    X$O18[h_o] <- 0
    X$Cl[h_c] <- 0
    X$air[h_a] <- 0

    message(
      "Tracer CV: ", labels[variant],
      "; rotation ", rot, " / 3"
    )
    fit <- model$sample(
      data = X,
      seed = cfg$seed + variant * 101L + rot * 1009L,
      chains = cfg$chains,
      parallel_chains = cfg$parallel_chains,
      iter_warmup = cfg$warmup,
      iter_sampling = cfg$sampling,
      adapt_delta = cfg$adapt_delta,
      max_treedepth = cfg$max_treedepth,
      refresh = 250
    )

    # Each masked observation is now a posterior predictive forecast.
    prediction_arrays <- c(
      "predicted_O18", "predicted_Cl", "predicted_air"
    )
    replicate_arrays <- c(
      "replicated_O18", "replicated_Cl", "replicated_air"
    )
    posterior_vars <- as.vector(outer(
      c(prediction_arrays, replicate_arrays),
      paste0("[", seq_len(nrow(d)), "]"),
      paste0
    ))
    draws <- posterior::as_draws_matrix(
      fit$draws(variables = posterior_vars)
    )
    pred_means <- list()
    pred_low <- list()
    pred_high <- list()
    for (j in seq_len(3L)) {
      nm <- c("O18", "Cl", "excess_air")[j]
      pred_cols <- paste0(prediction_arrays[j],
                          "[", seq_len(nrow(d)), "]")
      repl_cols <- paste0(replicate_arrays[j],
                          "[", seq_len(nrow(d)), "]")
      pred_means[[nm]] <- colMeans(
        draws[, pred_cols, drop = FALSE]
      )
      pred_low[[nm]] <- apply(
        draws[, repl_cols, drop = FALSE], 2,
        quantile, probs = .025
      )
      pred_high[[nm]] <- apply(
        draws[, repl_cols, drop = FALSE], 2,
        quantile, probs = .975
      )
    }
    CV[[fit_number]] <- held %>%
      mutate(
        model = labels[variant],
        variant = variant,
        truth = case_when(
          withheld_tracer == "O18" ~ d$O18[sample_index],
          withheld_tracer == "Cl" ~ d$Cl[sample_index],
          withheld_tracer == "excess_air" ~ d$xs_Air[sample_index]
        ),
        prediction = map2_dbl(
          withheld_tracer, sample_index,
          ~ pred_means[[.x]][.y]
        ),
        lower95 = map2_dbl(
          withheld_tracer, sample_index,
          ~ pred_low[[.x]][.y]
        ),
        upper95 = map2_dbl(
          withheld_tracer, sample_index,
          ~ pred_high[[.x]][.y]
        ),
        residual = truth - prediction,
        within95 = truth >= lower95 &
                   truth <= upper95
      )

    selected <- c(
      "O18_R", "O18_LSR", "Cl_R", "Cl_LSR",
      "air_R", "air_LSR", "extra_O18_sd",
      "extra_Cl_sd", "extra_air_sd",
      if (variant >= 2) "tau_O18_LSR",
      if (variant == 3) "tau_Cl_enrichment",
      paste0("f_LSR[", held$sample_index, "]")
    )
    diagnostics[[fit_number]] <- fit$summary(
      variables = selected
    ) %>%
      mutate(
        model = labels[variant],
        rotation = rot,
        .before = 1
      )
    mc <- fit$diagnostic_summary()
    if (any(diagnostics[[fit_number]]$rhat >
            cfg$maximum_Rhat, na.rm = TRUE) ||
        any(diagnostics[[fit_number]]$ess_bulk <
            cfg$minimum_bulk_ESS, na.rm = TRUE) ||
        any(mc$num_divergent > 0)) {
      warning(
        "Review convergence / HMC divergences: ",
        labels[variant], ", rotation ", rot
      )
    }
    # Save checkpoint after every fit.
    write_csv(
      bind_rows(CV),
      "outputs/tracer_only_CV_predictions.csv"
    )
    write_csv(
      bind_rows(diagnostics),
      "outputs/tracer_only_CV_convergence.csv"
    )
  }
}

pred <- bind_rows(CV)
by_tracer <- pred %>%
  group_by(model, withheld_tracer) %>%
  summarise(
    n = n(),
    RMSE = sqrt(mean(residual^2)),
    MAE = mean(abs(residual)),
    bias = mean(prediction - truth),
    coverage95 = mean(within95),
    .groups = "drop"
  ) %>%
  mutate(
    tracer_SD_in_dataset = case_when(
      withheld_tracer == "O18" ~ sd(d$O18, na.rm = TRUE),
      withheld_tracer == "Cl" ~ sd(d$Cl, na.rm = TRUE),
      TRUE ~ sd(d$xs_Air, na.rm = TRUE)
    ),
    normalised_RMSE = RMSE / tracer_SD_in_dataset
  )
overall <- by_tracer %>%
  group_by(model) %>%
  summarise(
    n_held_out = sum(n),
    mean_normalised_RMSE = mean(normalised_RMSE),
    mean_predictive_coverage95 = mean(coverage95),
    .groups = "drop"
  )

write_csv(
  by_tracer, "outputs/tracer_only_CV_scores_by_tracer.csv"
)
write_csv(overall, "outputs/tracer_only_CV_scores_overall.csv")
write_csv(
  bind_rows(diagnostics),
  "outputs/tracer_only_CV_convergence.csv"
)

# Fit every candidate to all available tracers to obtain its
# posterior source fraction distribution; this is NOT CV evidence.
if (isTRUE(cfg$fit_full_data)) {
  site_output <- list()
  global_output <- list()
  for (variant in cfg$variants) {
    X <- stan_base
    X$variant <- as.integer(variant)
    message("Full-data tracer mixing: ", labels[variant])
    fit <- model$sample(
      data = X,
      seed = cfg$seed + 100000L + variant,
      chains = cfg$chains,
      parallel_chains = cfg$parallel_chains,
      iter_warmup = cfg$warmup,
      iter_sampling = cfg$sampling,
      adapt_delta = cfg$adapt_delta,
      max_treedepth = cfg$max_treedepth,
      refresh = 250
    )
    site_output[[variant]] <- fit$summary(
      variables = "f_LSR"
    ) %>%
      mutate(
        sample_index = as.integer(
          str_match(variable, "\\[(\\d+)\\]")[, 2]
        ),
        model = labels[variant],
        .before = 1
      ) %>%
      inner_join(
        d %>% select(sample_index, Site_ID, Depth),
        by = "sample_index"
      )
    global_output[[variant]] <- fit$summary(
      variables = c(
        "O18_R", "O18_LSR",
        "Cl_R", "Cl_LSR",
        "air_R", "air_LSR",
        "extra_O18_sd", "extra_Cl_sd", "extra_air_sd",
        if (variant >= 2) "tau_O18_LSR",
        if (variant == 3) "tau_Cl_enrichment"
      )
    ) %>%
      mutate(model = labels[variant], .before = 1)
  }
  write_csv(
    bind_rows(site_output),
    "outputs/tracer_only_model_source_fractions.csv"
  )
  write_csv(
    bind_rows(global_output),
    "outputs/tracer_only_model_endmembers.csv"
  )
}

theme_set(theme_bw(base_size = 12))
p <- ggplot(pred, aes(truth, prediction)) +
  geom_abline(intercept = 0, slope = 1,
              linetype = "dashed", colour = "grey60") +
  geom_point(alpha = .75, size = 1.8) +
  facet_grid(model ~ withheld_tracer, scales = "free") +
  labs(
    title = "Can the other two tracers predict the withheld tracer?",
    subtitle = "Three rotated, identical masks for all Bayesian mixing models",
    x = "Measured value", y = "Posterior predictive mean"
  )
print(p)
ggsave(
  "plots/tracer_only_heldout_predictions.png",
  p, width = 12, height = 10, dpi = 300
)

p2 <- ggplot(
  by_tracer, aes(withheld_tracer, normalised_RMSE,
                 group = model, colour = model)
) +
  geom_line(linewidth = .85) +
  geom_point(size = 2.5) +
  labs(
    title = "Held-out tracer predictive accuracy",
    subtitle = "RMSE divided by observed tracer SD; smaller is better",
    x = "Withheld tracer",
    y = "Normalised RMSE",
    colour = "Mixing model"
  )
print(p2)
ggsave(
  "plots/tracer_only_normalised_CV_error.png",
  p2, width = 9, height = 6, dpi = 300
)

cat("\nHELD-OUT TRACER PERFORMANCE:\n")
print(by_tracer, n = Inf)
cat("\nEQUAL-TRACER-WEIGHTED COMPARISON:\n")
print(overall, n = Inf)
cat(
  "\nLIMITATIONS:\n",
  "* These are predictions of withheld TRACER observations, not DO.\n",
  "* The source-fraction model is conditional on study-derived endpoint\n",
  "  prior centres and on a TWO-SOURCE mixing hypothesis.\n",
  "* Variable O18 applies to the LSR component; variability in river\n",
  "  source signatures is not separately estimated here.\n",
  "* The chloride enrichment term is nonnegative and empirical;\n",
  "  its modelled presence is NOT proof of anthropogenic contamination.\n",
  "* If CV scores are similar, the extra mixing parameters are not\n",
  "  justified for predicting source fractions at other sites.\n",
  "* Compare posterior Rhat, ESS and predictive 95% coverage.\n",
  "* All outputs contain private per-well data; NEVER commit them.\n",
  sep = ""
)
