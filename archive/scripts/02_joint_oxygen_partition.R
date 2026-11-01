# Joint fit of source-associated DO deficit and aquifer oxygen reduction.
# This script DOES NOT use the old fixed-6.2-mg/L rainfall DO approach.
# Run scripts/01_three_tracer_mixing.R first.
library(tidyverse)

if (!file.exists("outputs/three_tracer_mixing.csv") &&
    file.exists("Oxygen-Reduction-Joint-Model/outputs/three_tracer_mixing.csv")) {
  setwd("Oxygen-Reduction-Joint-Model")
}
if (!file.exists("outputs/three_tracer_mixing.csv")) {
  stop("Run scripts/01_three_tracer_mixing.R first.")
}
dir.create("outputs", showWarnings = FALSE)
dir.create("plots", showWarnings = FALSE)

river_DO <- 7.5
rainfall_reference_DO <- 9.5  # Assumed PRE-infiltration near-saturation DO
bootstrap_replicates <- 150L
set.seed(20261030)

# Only effective rainfall-source DO (reference - source-associated deficit)
# is identified by DO observations. The location of reduction is not.
input <- read_csv(
  "outputs/three_tracer_mixing.csv", show_col_types = FALSE,
  col_types = cols(Site_ID = col_character(), site = col_character(),
                   .default = col_guess())
)
d <- input %>%
  filter(mixing_reliable, is.finite(DO), DO > 0,
         is.finite(MRT), MRT >= 0,
         is.finite(rainfall_contribution),
         is.finite(river_contribution)) %>%
  mutate(fR = river_contribution, fP = rainfall_contribution)

stopifnot(
  nrow(d) == 29, n_distinct(d$Site_ID) == 29,
  sum(d$fP >= 0.5) == 5,
  all(abs(d$fP + d$fR - 1) < 1e-9)
)

# DO on its measured concentration scale.
predict_DO <- function(dat, eff, k, kinetic) {
  mixed <- dat$fR * river_DO + dat$fP * eff
  if (kinetic == "First order") {
    mixed * exp(-k * dat$MRT)
  } else if (kinetic == "Zero order") {
    pmax(0, mixed - k * dat$MRT)
  } else stop("Unknown kinetic model")
}

# Exact piecewise-quadratic minimisation in k for the CLIPPED zero-order model.
# A single optimize() call across 0..10 could settle on a high-k plateau.
zero_k_at_effective <- function(dat, eff) {
  initial <- dat$fR * river_DO + dat$fP * eff
  t <- dat$MRT
  positive_age <- is.finite(t) & t > 0
  knots <- sort(unique(c(0, 10, pmin(10, initial[positive_age] / t[positive_age]))))
  candidates <- knots
  for (j in seq_len(length(knots) - 1L)) {
    lo <- knots[j]
    hi <- knots[j + 1L]
    active <- positive_age & initial > ((lo + hi) / 2) * t
    denominator <- sum(t[active]^2)
    if (denominator > 0) {
      candidate <- sum(t[active] * (initial[active] - dat$DO[active])) / denominator
      candidates <- c(candidates, max(lo, min(hi, candidate)))
    }
  }
  SSE <- vapply(candidates, function(k) {
    sum((dat$DO - pmax(0, initial - k * t))^2)
  }, numeric(1))
  winner <- which.min(SSE)
  c(eff = eff, k = candidates[winner], SSE = SSE[winner])
}

# For first order, the least-squares optimum effective DO can be solved
# analytically for any fixed rate k (subject to physical bounds 0..9.5).
first_at_k <- function(dat, k, fixed_eff = NA_real_) {
  attenuation <- exp(-k * dat$MRT)
  river_component <- dat$fR * river_DO * attenuation
  rainfall_coefficient <- dat$fP * attenuation
  eff <- if (is.finite(fixed_eff)) {
    fixed_eff
  } else {
    denominator <- sum(rainfall_coefficient^2)
    if (denominator > 0) {
      sum(rainfall_coefficient * (dat$DO - river_component)) / denominator
    } else 0
  }
  eff <- max(0, min(rainfall_reference_DO, eff))
  SSE <- sum((dat$DO - river_component - rainfall_coefficient * eff)^2)
  c(eff = eff, k = k, SSE = SSE)
}

# Scan a grid before refining its minimum; explicitly compare boundaries.
profile_minimum <- function(grid, objective) {
  evaluated <- lapply(grid, objective)
  SSE <- vapply(evaluated, function(x) x[["SSE"]], numeric(1))
  winner <- which.min(SSE)
  lower <- grid[max(1, winner - 1)]
  upper <- grid[min(length(grid), winner + 1)]
  refined <- if (upper > lower) {
    optimize(function(x) objective(x)[["SSE"]],
             c(lower, upper), tol = 1e-8)$minimum
  } else grid[winner]
  choices <- lapply(unique(c(grid[winner], refined, grid[1],
                             tail(grid, 1))), objective)
  choices[[which.min(vapply(choices, function(x) x[["SSE"]], numeric(1)))]]
}

fit_model <- function(dat, kinetic, model) {
  stopifnot(model %in% c("aquifer_only", "infiltration_only", "combined"))
  reference <- rainfall_reference_DO

  if (kinetic == "First order") {
    if (model == "infiltration_only") {
      optimum <- first_at_k(dat, 0)
    } else {
      fixed <- if (model == "aquifer_only") reference else NA_real_
      objective <- function(k) first_at_k(dat, k, fixed)
      optimum <- profile_minimum(seq(0, 10, length.out = 801), objective)
    }
  } else if (kinetic == "Zero order") {
    if (model == "aquifer_only") {
      optimum <- zero_k_at_effective(dat, reference)
    } else if (model == "infiltration_only") {
      eff <- first_at_k(dat, 0)[["eff"]]
      optimum <- c(
        eff = eff, k = 0,
        SSE = sum((dat$DO - predict_DO(dat, eff, 0, kinetic))^2)
      )
    } else {
      objective <- function(eff) zero_k_at_effective(dat, eff)
      optimum <- profile_minimum(
        seq(0, reference, length.out = 101), objective
      )
      # Include exact nested alternatives in the candidate set.
      inf_eff <- first_at_k(dat, 0)[["eff"]]
      choices <- list(optimum, objective(inf_eff), objective(reference))
      optimum <- choices[[which.min(vapply(
        choices, function(x) x[["SSE"]], numeric(1)
      ))]]
    }
  } else stop("Unknown kinetic model")

  tibble(
    kinetic = kinetic,
    model = model,
    reference_rainfall_DO = reference,
    source_associated_loss_mg_L = reference - optimum[["eff"]],
    effective_rainfall_DO_mg_L = optimum[["eff"]],
    aquifer_rate = optimum[["k"]],
    aquifer_rate_units = if_else(
      kinetic == "First order", "1/year", "mg/L/year"
    ),
    n = nrow(dat),
    SSE = optimum[["SSE"]],
    RMSE_mg_L = sqrt(optimum[["SSE"]] / nrow(dat))
  )
}

model_grid <- tidyr::crossing(
  kinetic = c("First order", "Zero order"),
  model = c("aquifer_only", "infiltration_only", "combined")
)
run_fits <- function(dat) {
  purrr::map2_dfr(
    model_grid$kinetic, model_grid$model,
    ~fit_model(dat, .x, .y)
  )
}
fits <- run_fits(d)
write_csv(fits, "outputs/joint_model_fits.csv")

# An optimisation guard: the combined model nests the two simpler models.
check_nesting <- function(table) {
  for (kin in unique(table$kinetic)) {
    sub <- table %>% filter(kinetic == kin)
    errors <- setNames(sub$SSE, sub$model)
    if (errors[["combined"]] >
        min(errors[c("aquifer_only", "infiltration_only")]) + 1e-5) {
      stop("Nested model mismatch: ", kin)
    }
  }
}
check_nesting(fits)

# Leave one WELL out: refit both parameters on other wells, then predict
# the held-out DO. The held-out measured DO is not used in model fitting.
well_ids <- unique(d$Site_ID)
LOO_fits <- purrr::map_dfr(well_ids, function(id) {
  training <- d %>% filter(Site_ID != id)
  run_fits(training) %>% mutate(held_out_Site_ID = id)
})
write_csv(LOO_fits, "outputs/joint_model_LOOWO_fits.csv")

predictions <- purrr::map_dfr(seq_len(nrow(LOO_fits)), function(i) {
  fit <- LOO_fits[i, ]
  held <- d %>% filter(Site_ID == fit$held_out_Site_ID)
  predicted <- predict_DO(
    held, fit$effective_rainfall_DO_mg_L,
    fit$aquifer_rate, fit$kinetic
  )
  initial_DO <- held$fR * river_DO +
    held$fP * fit$effective_rainfall_DO_mg_L
  tibble(
    Site_ID = held$Site_ID,
    site = held$site,
    MRT = held$MRT,
    rainfall_fraction = held$fP,
    DO_initial = initial_DO,
    DO_observed = held$DO,
    kinetic = fit$kinetic,
    model = fit$model,
    DO_predicted = predicted,
    residual = held$DO - predicted
  )
})
write_csv(predictions, "outputs/joint_model_LOOWO_predictions.csv")

performance <- predictions %>%
  group_by(kinetic, model) %>%
  summarise(
    n = n(),
    RMSE_mg_L = sqrt(mean(residual^2)),
    MAE_mg_L = mean(abs(residual)),
    bias_mg_L = mean(residual),
    .groups = "drop"
  )
write_csv(performance, "outputs/joint_model_LOOWO_performance.csv")

# Conditional nonparametric well bootstrap. MRT, source endmembers and
# the tracer mixing fractions are held fixed (not total model uncertainty).
bootstrap_fits <- purrr::map_dfr(seq_len(bootstrap_replicates), function(b) {
  draw <- sample(well_ids, length(well_ids), replace = TRUE)
  boot <- bind_rows(lapply(draw, function(id) d %>% filter(Site_ID == id)))
  purrr::map_dfr(c("First order", "Zero order"), function(kin) {
    fit_model(boot, kin, "combined") %>% mutate(iteration = b)
  })
})
bootstrap_summary <- bootstrap_fits %>%
  group_by(kinetic) %>%
  summarise(
    bootstrap_n = n(),
    loss_median = median(source_associated_loss_mg_L),
    loss_lower_95 = unname(quantile(source_associated_loss_mg_L, 0.025)),
    loss_upper_95 = unname(quantile(source_associated_loss_mg_L, 0.975)),
    effective_DO_median = median(effective_rainfall_DO_mg_L),
    effective_DO_lower_95 =
      unname(quantile(effective_rainfall_DO_mg_L, 0.025)),
    effective_DO_upper_95 =
      unname(quantile(effective_rainfall_DO_mg_L, 0.975)),
    rate_median = median(aquifer_rate),
    rate_lower_95 = unname(quantile(aquifer_rate, 0.025)),
    rate_upper_95 = unname(quantile(aquifer_rate, 0.975)),
    pct_at_zero_loss = 100 * mean(source_associated_loss_mg_L < 1e-4),
    .groups = "drop"
  )
write_csv(bootstrap_summary, "outputs/joint_model_bootstrap_summary.csv")

# Profile the best aquifer rate over EFFECTIVE rainfall DO, to assess
# identifiability of the fitted source-associated oxygen deficit.
profile <- map_dfr(c("First order", "Zero order"), function(kin) {
  map_dfr(seq(0, rainfall_reference_DO, length.out = 101), function(eff) {
    if (kin == "First order") {
      objective <- function(k) first_at_k(d, k, fixed_eff = eff)
      optimum <- profile_minimum(
        seq(0, 10, length.out = 801), objective
      )
    } else {
      optimum <- zero_k_at_effective(d, eff)
    }
    tibble(
      kinetic = kin,
      effective_rainfall_DO_mg_L = eff,
      aquifer_rate = optimum[["k"]],
      RMSE_mg_L = sqrt(optimum[["SSE"]] / nrow(d))
    )
  })
}) %>%
  group_by(kinetic) %>%
  mutate(delta_RMSE_mg_L = RMSE_mg_L - min(RMSE_mg_L)) %>%
  ungroup()
write_csv(profile, "outputs/joint_model_effective_DO_profile.csv")

p1 <- ggplot(
  performance, aes(model, RMSE_mg_L, group = kinetic, colour = kinetic)
) +
  geom_line(linewidth = 0.75) +
  geom_point(size = 3) +
  labs(
    x = "Model", y = "Leave-one-well-out RMSE (mg/L)",
    colour = "Kinetics"
  ) + theme_bw()
ggsave(
  "plots/joint_model_LOOWO_RMSE.png",
  p1, width = 8, height = 5, dpi = 220
)

p2 <- ggplot(
  profile,
  aes(effective_rainfall_DO_mg_L, delta_RMSE_mg_L, colour = kinetic)
) +
  geom_line(linewidth = 1) +
  labs(
    x = "Effective rainfall recharge DO (mg/L)",
    y = "Increase in fitted RMSE above optimum (mg/L)",
    colour = "Kinetics"
  ) + theme_bw()
ggsave(
  "plots/joint_model_effective_DO_profile.png",
  p2, width = 8, height = 5, dpi = 220
)

cat("\nJoint model fits:\n")
print(fits %>% filter(model == "combined"))
cat("\nHeld-out predictive performance:\n")
print(performance)
cat("\nConditional bootstrap uncertainty:\n")
print(bootstrap_summary)
