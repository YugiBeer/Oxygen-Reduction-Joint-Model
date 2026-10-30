# 10_groundwater_age_distribution_sensitivity.R
#
# Fixed, no-new-prior sensitivity test of the age distribution underlying
# the first-order oxygen depletion term.
#
# Model:
#   DO_i = ((1-fLSR_i)*7.5 + fLSR_i*P) * A(MRT_i, k, alpha)
#
#   Gamma age distribution: shape alpha, mean MRT_i
#     A = (1 + k*MRT_i/alpha)^(-alpha)
#   Piston-flow age distribution (alpha -> Inf)
#     A = exp(-k*MRT_i)
#
# P is effective land-surface recharge (LSR) DO at aquifer entry,
# fitted within 0..9.5 mg/L. No river-recharge loss is fitted.
# Each alpha is FIXED in advance, not estimated from the groundwater
# DO data. The distribution from which reported MRT was originally
# derived MUST be checked before physical interpretation.
#
# The gamma shape controls age variation:
#    coefficient of variation of age = 1/sqrt(alpha)
#    alpha = 1: exponential age distribution
#    alpha -> Inf: piston flow
#
# This is a sensitivity experiment, NOT a Bayesian fit or age-tracer
# inversion. Groundwater age distributions of river and LSR components
# are treated as identical, which is another untested assumption.
#
# Run after scripts/01_three_tracer_mixing.R.
library(tidyverse)
setwd(here::here("Oxygen-Reduction-Joint-Model"))

# ===============================================================
# SETTINGS (all editable)
# ===============================================================
cfg <- list(
  river_DO_reference = 7.5,
  maximum_LSR_DO = 9.5,
  gamma_shapes = c(Inf, 100, 50, 20, 10, 5, 3, 2, 1),
  max_k = 3.0,           # year^-1; larger than usual to permit alpha=1
  k_grid_points = 500L,
  folds = 5L,
  seed = 20261030L,
  compare_without_2801 = TRUE,
  font_size = 12,
  plot_dpi = 300
)

dir.create("outputs", showWarnings = FALSE, recursive = TRUE)
dir.create("plots", showWarnings = FALSE, recursive = TRUE)
f <- "outputs/three_tracer_mixing.csv"
if (!file.exists(f)) stop("Run scripts/01_three_tracer_mixing.R first.")

d <- read_csv(
  f, show_col_types = FALSE,
  col_types = cols(Site_ID = col_character(),
                   .default = col_guess())
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
    fR = river_contribution
  ) %>%
  mutate(group = case_when(
    str_detect(Site_ID, "Eastbourne") ~ "Eastbourne",
    str_detect(Site_ID, "Gilbertson") ~ "Gilbertson",
    str_detect(Site_ID, "Flaxmere") ~ "Flaxmere",
    TRUE ~ Site_ID
  ))
stopifnot(
  nrow(d) == 29L,
  n_distinct(d$Site_ID) == 29L,
  sum(d$fLSR >= 0.5) == 5L,
  all(abs(d$fR + d$fLSR - 1) < 1e-8)
)

# A(MRT, k, alpha) = Laplace transform of the gamma age distribution.
# Use the log form for numerical stability, especially for large alpha.
age_attenuation <- function(MRT, k, alpha) {
  if (is.infinite(alpha)) {
    exp(-k * MRT)
  } else {
    stopifnot(is.finite(alpha), alpha > 0)
    exp(-alpha * log1p(k * MRT / alpha))
  }
}

prediction <- function(dat, P, k, alpha) {
  (dat$fR * cfg$river_DO_reference + dat$fLSR * P) *
    age_attenuation(dat$MRT, k, alpha)
}

# Profile P out analytically at a fixed k; enforce 0 <= P <= 9.5.
fit_age_model <- function(dat, alpha) {
  score <- function(k) {
    att <- age_attenuation(dat$MRT, k, alpha)
    river <- dat$fR * cfg$river_DO_reference * att
    w <- dat$fLSR * att
    denom <- sum(w^2)
    P <- if (denom > 0) {
      sum(w * (dat$DO - river)) / denom
    } else 0
    P <- max(0, min(cfg$maximum_LSR_DO, P))
    predicted <- river + w * P
    c(
      P = P, k = k,
      SSE = sum((dat$DO - predicted)^2)
    )
  }

  # The profiled objective can have local minima. Search a dense grid,
  # refine only around the best grid point, and include boundaries.
  ks <- c(0, exp(seq(
    log(1e-6), log(cfg$max_k),
    length.out = cfg$k_grid_points
  )))
  evaluated <- lapply(ks, score)
  ss <- vapply(evaluated, function(x) x[["SSE"]], numeric(1))
  winner <- which.min(ss)
  low <- ks[max(1L, winner - 1L)]
  high <- ks[min(length(ks), winner + 1L)]
  refined <- if (high > low) {
    optimize(
      function(k) score(k)[["SSE"]],
      interval = c(low, high),
      tol = 1e-9
    )$minimum
  } else {
    ks[winner]
  }
  candidates <- lapply(
    unique(c(0, ks[winner], refined, cfg$max_k)),
    score
  )
  best <- candidates[[which.min(vapply(
    candidates, function(x) x[["SSE"]], numeric(1)
  ))]]
  list(
    P = unname(best[["P"]]),
    k = unname(best[["k"]]),
    SSE = unname(best[["SSE"]]),
    RMSE = sqrt(best[["SSE"]] / nrow(dat))
  )
}

# Hold out entire Eastbourne/Gilbertson/Flaxmere groups together.
# Keep identical folds for every gamma-shape scenario, including
# the 28-well sensitivity omitting 2801.
set.seed(cfg$seed)
groups <- split(seq_len(nrow(d)), d$group)
order_group <- sample(seq_along(groups))
order_group <- order_group[order(
  -lengths(groups)[order_group]
)]
folds <- vector("list", cfg$folds)
for (gi in order_group) {
  j <- which.min(lengths(folds))
  folds[[j]] <- c(folds[[j]], groups[[gi]])
}
stopifnot(
  all(lengths(folds) > 0),
  identical(sort(unlist(folds)), seq_len(nrow(d)))
)

evaluate_one <- function(sub, selection, alpha) {
  fit <- fit_age_model(sub, alpha)
  all_indices <- which(d$Site_ID %in% sub$Site_ID)
  pred <- map_dfr(seq_along(folds), function(j) {
    held_ix <- intersect(folds[[j]], all_indices)
    if (!length(held_ix)) return(tibble())
    held <- d[held_ix, , drop = FALSE]
    training <- sub %>%
      filter(!Site_ID %in% held$Site_ID)
    fit_train <- fit_age_model(training, alpha)
    tibble(
      sample = selection,
      gamma_shape = alpha,
      fold = j,
      Site_ID = held$Site_ID,
      MRT = held$MRT,
      LSR_fraction = held$fLSR,
      observed_DO = held$DO,
      CV_predicted_DO = prediction(
        held, fit_train$P, fit_train$k, alpha
      )
    )
  })

  # Leave-one-well-out, used only as an additional diagnostic.
  loo <- map_dfr(seq_len(nrow(sub)), function(j) {
    held <- sub[j, , drop = FALSE]
    training <- sub[-j, , drop = FALSE]
    fit_train <- fit_age_model(training, alpha)
    tibble(
      Site_ID = held$Site_ID,
      LOOWO_predicted_DO = prediction(
        held, fit_train$P, fit_train$k, alpha
      )
    )
  })
  pred <- left_join(pred, loo, by = "Site_ID") %>%
    mutate(
      CV_residual = observed_DO - CV_predicted_DO
    )

  summary <- tibble(
    sample = selection,
    gamma_shape = alpha,
    n = nrow(sub),
    effective_LSR_DO = fit$P,
    apparent_k_per_year = fit$k,
    fitted_RMSE_mg_L = fit$RMSE,
    grouped_CV_RMSE_mg_L = sqrt(mean(pred$CV_residual^2)),
    grouped_CV_MAE_mg_L = mean(abs(pred$CV_residual)),
    LOOWO_RMSE_mg_L = sqrt(mean(
      (pred$observed_DO - pred$LOOWO_predicted_DO)^2
    ))
  )
  list(summary = summary, predictions = pred)
}

sets <- list("All 29" = d)
if (isTRUE(cfg$compare_without_2801)) {
  sets[["Excluding 2801"]] <- d %>%
    filter(Site_ID != "2801")
}
results <- imap(sets, function(sub, name) {
  map(cfg$gamma_shapes, function(a) {
    evaluate_one(sub, name, a)
  })
})
flat <- purrr::flatten(results)
summary_table <- bind_rows(map(flat, "summary"))
predictions <- bind_rows(map(flat, "predictions"))

write_csv(
  summary_table,
  "outputs/groundwater_age_distribution_sensitivity.csv"
)
write_csv(
  predictions,
  "outputs/groundwater_age_distribution_CV_predictions.csv"
)

# Diagnostic plots use categorical age-model labels; an infinite shape
# does not have a well-defined position on a numeric log axis.
shape_levels <- as.character(cfg$gamma_shapes)
shape_labels <- ifelse(
  is.infinite(cfg$gamma_shapes),
  "Piston",
  paste0("Gamma ", cfg$gamma_shapes)
)
format_shape <- function(a) {
  factor(
    as.character(a),
    levels = shape_levels,
    labels = shape_labels
  )
}

plot_data <- summary_table %>%
  mutate(shape_model = format_shape(gamma_shape))
theme_set(theme_bw(base_size = cfg$font_size))

p_cv <- ggplot(
  plot_data,
  aes(shape_model, grouped_CV_RMSE_mg_L, group = sample)
) +
  geom_line(aes(colour = sample), linewidth = .9) +
  geom_point(aes(colour = sample), size = 2.8) +
  coord_flip() +
  labs(
    title = "Does an age distribution improve DO prediction?",
    subtitle = "Same fixed alpha scenarios and grouped folds",
    x = "Groundwater age model",
    y = "Grouped five-fold CV RMSE (mg/L)",
    colour = "Well selection"
  )
print(p_cv)
ggsave(
  "plots/groundwater_age_distribution_CV.png", p_cv,
  width = 9, height = 6, dpi = cfg$plot_dpi
)

p_k <- ggplot(
  plot_data,
  aes(shape_model, apparent_k_per_year, group = sample)
) +
  geom_line(aes(colour = sample), linewidth = .9) +
  geom_point(aes(colour = sample), size = 2.8) +
  coord_flip() +
  labs(
    title = "Sensitivity of apparent DO depletion rate to age mixing",
    subtitle = "Gamma shapes fixed, not inferred from DO",
    x = "Groundwater age model",
    y = expression(k~(year^{-1})),
    colour = "Well selection"
  )
print(p_k)
ggsave(
  "plots/groundwater_age_distribution_k.png", p_k,
  width = 9, height = 6, dpi = cfg$plot_dpi
)

p_resid <- predictions %>%
  filter(
    sample == "All 29",
    (is.infinite(gamma_shape) |
       gamma_shape %in% c(10, 5, 2, 1))
  ) %>%
  mutate(shape_model = format_shape(gamma_shape)) %>%
  ggplot(aes(MRT, CV_residual)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(aes(colour = LSR_fraction), size = 2.2) +
  facet_wrap(~ shape_model, ncol = 2) +
  labs(
    title = "Held-out DO residuals versus reported MRT",
    x = "Reported mean residence time (years)",
    y = "Observed - predicted DO (mg/L)",
    colour = "LSR fraction"
  )
print(p_resid)
ggsave(
  "plots/groundwater_age_distribution_residuals.png", p_resid,
  width = 10, height = 8, dpi = cfg$plot_dpi
)

cat("\nAGE-DISTRIBUTION SENSITIVITY:\n")
print(summary_table, n = Inf)
cat(
  "\nIMPORTANT SCIENTIFIC LIMITATIONS:\n",
  "* These are fixed-shape SENSITIVITY SCENARIOS, not inferred ages.\n",
  "* Check how the reported MRTs were originally estimated and which\n",
  "  lumped-parameter model was used; otherwise shapes may be invalid.\n",
  "* We assumed the same gamma age distribution applies to river and\n",
  "  LSR contributions; the age of each component may differ.\n",
  "* MRT uncertainty, uncertainty in the source fractions, and possible\n",
  "  source-age correlations have NOT been propagated.\n",
  "* A broader age distribution increases young-water contributions,\n",
  "  generally increasing fitted k without identifying where oxygen\n",
  "  consumption occurred along the recharge-to-well pathway.\n",
  "* No Bayesian priors or well-specific reduction rates were introduced.\n",
  sep = ""
)
