# 07_infiltration_loss_sensitivity.R (legacy filename retained for compatibility)
# In this analysis, "infiltration" means LAND-SURFACE recharge only.
#
# Purpose
# -------
# Test the trade-off between the river-derived recharge DO deficit (L_R)
# and age-associated oxygen depletion (k). Quantify the minimum
# river-recharge DO deficit required at each well IF k = 0.
# Compare fixed/stratified recharge DO alternatives by grouped 5-fold CV.
#
# Run after scripts/01_three_tracer_mixing.R.
#
# DO_i = [fR_i*(river_DO - L_R) + fP_i*P] * exp(-k*MRT_i)
# P = effective LAND-SURFACE recharge DO entering the aquifer.
#       It ALREADY includes DO loss during soil/unsaturated-zone infiltration.
#
# KEY LIMITATIONS
# ---------------
# - k is an apparent age-associated coefficient. It cannot identify
#   where reduction occurred (recharge pathway versus aquifer).
# - "Confined" means sampled well status, NOT recharge location.
# - No spatial coordinates or direct recharge-entry DO measurements available.
# - An age-correlated recharge-entry-DO model can be mathematically
#   IDENTICAL to a model with age-dependent reduction.
# - River-source DO deficit is NOT independently measured riverbed loss.
# - No Bayesian posterior in this script: bounded least squares and
#   grouped CV are transparent, reproducible sensitivity diagnostics.
#
# All groundwater sample data remain in the user's local data directory.
library(tidyverse)
setwd(here::here("Oxygen-Reduction-Joint-Model"))

# ======================================================================
# USER SETTINGS
# ======================================================================
cfg <- list(
  river_DO = 7.5,              # observed river daily-mean reference, mg/L
  river_DO_scenarios = c(5.6, 7.5, 9.4),
  rain_DO_max = 9.5,           # effective rain DO upper bound, mg/L
  river_DO_deficit_step = 0.25, # mg/L; relative to monitored river DO
  max_k = 1.0,                # year^-1
  n_cv_folds = 5L,
  river_dominated_threshold = 0.80,
  old_groundwater_age = 10,   # years
  seed = 20261030L,
  font_size = 12,
  plot_dpi = 300
)

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
         is.finite(river_contribution),
         !is.na(confined)) %>%
  transmute(
    Site_ID, DO, MRT,
    fR = river_contribution, fP = rainfall_contribution,
    confined = as.integer(confined)
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
  all(d$confined %in% c(0L, 1L)),
  all(abs(d$fR + d$fP - 1) < 1e-8)
)

# Grouped five-fold CV: related samples are withheld together.
# The heuristic grouping is not spatial CV.
set.seed(cfg$seed)
grouped <- split(seq_len(nrow(d)), d$cv_group)
order_groups <- sample(seq_along(grouped))
order_groups <- order_groups[order(
  -lengths(grouped)[order_groups]
)]
folds <- vector("list", cfg$n_cv_folds)
for (i in order_groups) {
  j <- which.min(lengths(folds))
  folds[[j]] <- c(folds[[j]], grouped[[i]])
}
stopifnot(
  identical(sort(unlist(folds)), seq_len(nrow(d))),
  all(lengths(folds) > 0)
)

# ======================================================================
# FIXED RIVER-ENTRY DO DEFICIT: refit land-surface DO and k at each L_R.
# ======================================================================
# Multiple start values prevent a suboptimal local solution at bounds.
fit_fixed_deficit <- function(dat, river_deficit) {
  residual <- function(z) {
    P <- z[1]
    k <- z[2]
    predicted <- (
      dat$fR * (cfg$river_DO - river_deficit) + dat$fP * P
    ) * exp(-k * dat$MRT)
    dat$DO - predicted
  }
  starts <- list(c(8, .07), c(4, .03), c(9, .15))
  candidates <- lapply(starts, function(z0) {
    optim(
      z0,
      fn = function(z) sum(residual(z)^2),
      method = "L-BFGS-B",
      lower = c(0, 0),
      upper = c(cfg$rain_DO_max, cfg$max_k),
      control = list(maxit = 1000)
    )
  })
  best <- candidates[[which.min(vapply(
    candidates, function(x) x$value, numeric(1)
  ))]]
  tibble(
    effective_rain_DO = best$par[1],
    k = best$par[2],
    fit_rmse = sqrt(mean(residual(best$par)^2))
  )
}

profile_one <- function(drop_2801) {
  dat <- if (drop_2801) {
    d %>% filter(Site_ID != "2801")
  } else {
    d
  }
  loss_grid <- seq(0, cfg$river_DO, by = cfg$river_DO_deficit_step)
  purrr::map_dfr(loss_grid, function(L) {
    fit <- fit_fixed_deficit(dat, L)
    # Within-fold refitting: held-out wells do NOT enter the fit.
    predictions <- purrr::map_dfr(seq_along(folds), function(j) {
      held_ids <- intersect(folds[[j]], which(d$Site_ID %in% dat$Site_ID))
      if (!length(held_ids)) return(tibble())
      train <- dat %>% filter(!Site_ID %in% d$Site_ID[held_ids])
      held <- d[held_ids, , drop = FALSE]
      fitted <- fit_fixed_deficit(train, L)
      tibble(
        Site_ID = held$Site_ID,
        observed = held$DO,
        predicted = (
          held$fR * (cfg$river_DO - L) +
            held$fP * fitted$effective_rain_DO
        ) * exp(-fitted$k * held$MRT)
      )
    })
    tibble(
      exclude_2801 = drop_2801,
      n = nrow(dat),
      river_recharge_DO_deficit_mg_L = L,
      river_entry_DO_mg_L = cfg$river_DO - L,
      effective_rain_DO_mg_L = fit$effective_rain_DO,
      k_per_year = fit$k,
      fit_RMSE_mg_L = fit$fit_rmse,
      grouped_CV_RMSE_mg_L = sqrt(
        mean((predictions$observed - predictions$predicted)^2)
      ),
      grouped_CV_MAE_mg_L = mean(
        abs(predictions$observed - predictions$predicted)
      )
    )
  })
}

profile <- bind_rows(profile_one(FALSE), profile_one(TRUE))
write_csv(profile, "outputs/infiltration_loss_profile.csv")

# ======================================================================
# REQUIRED RIVER-SOURCE DO DEFICITS, assuming NO age-dependent depletion
# ======================================================================
# DO_i = fR_i * (R - L_R,i) + fP_i * P_i, 0 <= P_i <= 9.5.
#
# To MINIMISE the river-source DO deficit, set land-surface P_i=0.
#   L_R,min = max(0, R - DO_i/fR_i) .
# This is a necessary lower bound on the river DO deficit given the\n# assumed monitored-river reference. It does NOT locate the depletion.\n# Positive P_i increases the required deficit.
# For a sample with DO_i > fR_i*R + fP_i*9.5, no feasible solution
# exists with these fixed source upper bounds (e.g. site 2801 at R=7.5).
required <- purrr::map_dfr(
  cfg$river_DO_scenarios, function(R) {
    d %>% mutate(
      river_reference_DO_mg_L = R,
      recharge_entry_only_feasible = (
        DO <= fR * R + fP * cfg$rain_DO_max + 1e-10
      ),
      min_river_recharge_DO_deficit_mg_L = if_else(
        fR > 0,
        pmax(0, R - DO / pmax(fR, 1e-12)),
        NA_real_
      ),
      min_river_source_DO_deficit_fraction = min_river_recharge_DO_deficit_mg_L / R,
      river_dominated = fR >= cfg$river_dominated_threshold,
      older_group = MRT >= cfg$old_groundwater_age
    ) %>%
      select(
        river_reference_DO_mg_L, Site_ID, DO, MRT, fR, fP,
        confined, river_dominated, older_group,
        recharge_entry_only_feasible,
        min_river_recharge_DO_deficit_mg_L,
        min_river_source_DO_deficit_fraction
      )
  }
)
write_csv(required, "outputs/infiltration_minimum_loss_per_well.csv")

summary_required <- required %>%
  filter(river_dominated, older_group) %>%
  group_by(river_reference_DO_mg_L) %>%
  summarise(
    n = n(),
    median_min_river_source_DO_deficit_fraction = median(
      min_river_source_DO_deficit_fraction, na.rm = TRUE
    ),
    n_min_loss_ge_75pct = sum(
      min_river_source_DO_deficit_fraction >= .75, na.rm = TRUE
    ),
    n_min_loss_ge_90pct = sum(
      min_river_source_DO_deficit_fraction >= .90, na.rm = TRUE
    ),
    .groups = "drop"
  )
write_csv(summary_required, "outputs/infiltration_required_loss_summary.csv")

# ======================================================================
# GROUPED CV: common versus confinement-specific recharge DO WITHOUT
# AGE EFFECT, compared against first-order model with river R=7.5.
# ======================================================================
# No site-level inlet DO values are estimated; that would give an unfair
# advantage to the variable-loss model and cannot predict new wells.
#
# "shared" fits effective river entry DO R and rain DO P, with k=0.
# "confined" fits R_unconfined, R_confined and common rain DO P, k=0.
# "both sources" fits confined/unconfined R and P independently, k=0.
# "age" fixes river reference R=7.5 and fits P,k, with no explicit
# river-source deficit parameter. P accounts for land-surface pathway loss.
#
# A confined / unconfined contrast is a weak proxy for spatial variance.
# It must NOT be interpreted as the actual location of oxygen loss.
estimate_model <- function(dat, model) {
  if (model == "Shared recharge-entry DO (k=0)") {
    bounds_min <- c(0, 0)
    bounds_max <- c(cfg$river_DO, cfg$rain_DO_max)
    starts <- list(c(4, 5), c(7, 3), c(2, 8))
    forward <- function(z, a) a$fR * z[1] + a$fP * z[2]
  } else if (model == "River-entry DO varies by confinement (k=0)") {
    bounds_min <- c(0, 0, 0)
    bounds_max <- c(cfg$river_DO, cfg$river_DO, cfg$rain_DO_max)
    starts <- list(c(6, 2, 5), c(4, 1, 8), c(7, 5, 3))
    forward <- function(z, a) {
      a$fR * if_else(a$confined == 1L, z[2], z[1]) + a$fP * z[3]
    }
  } else if (model == "Both recharge-entry DO values vary (k=0)") {
    bounds_min <- c(0, 0, 0, 0)
    bounds_max <- c(
      cfg$river_DO, cfg$river_DO,
      cfg$rain_DO_max, cfg$rain_DO_max
    )
    starts <- list(
      c(6, 2, 6, 2), c(4, 1, 8, 1), c(7, 5, 3, 6)
    )
    forward <- function(z, a) {
      a$fR * if_else(a$confined == 1L, z[2], z[1]) +
        a$fP * if_else(a$confined == 1L, z[4], z[3])
    }
  } else if (model == "First-order age-dependent") {
    bounds_min <- c(0, 0)
    bounds_max <- c(cfg$rain_DO_max, cfg$max_k)
    starts <- list(c(8, .07), c(5, .04), c(9, .12))
    forward <- function(z, a) (
      a$fR * cfg$river_DO + a$fP * z[1]
    ) * exp(-z[2] * a$MRT)
  } else {
    stop("Unrecognised model")
  }

  result <- lapply(starts, function(z0) optim(
    z0,
    fn = function(z) sum((dat$DO - forward(z, dat))^2),
    method = "L-BFGS-B",
    lower = bounds_min, upper = bounds_max,
    control = list(maxit = 1000)
  ))
  best <- result[[which.min(vapply(
    result, function(x) x$value, numeric(1)
  ))]]
  list(parameters = best$par, predict = function(a) {
    forward(best$par, a)
  })
}

models <- c(
  "Shared recharge-entry DO (k=0)",
  "River-entry DO varies by confinement (k=0)",
  "Both recharge-entry DO values vary (k=0)",
  "First-order age-dependent"
)

CV_predictions <- purrr::map_dfr(c(FALSE, TRUE), function(drop_2801) {
  ids <- which(!drop_2801 | d$Site_ID != "2801")
  purrr::map_dfr(models, function(model) {
    purrr::map_dfr(seq_along(folds), function(j) {
      held_ids <- intersect(folds[[j]], ids)
      if (!length(held_ids)) return(tibble())
      training_ids <- setdiff(ids, held_ids)
      fit <- estimate_model(d[training_ids, ], model)
      held <- d[held_ids, ]
      tibble(
        exclude_2801 = drop_2801,
        model = model,
        fold = j,
        Site_ID = held$Site_ID,
        DO = held$DO,
        MRT = held$MRT,
        fR = held$fR,
        predicted_DO = fit$predict(held)
      )
    })
  })
})
write_csv(CV_predictions, "outputs/infiltration_model_CV_predictions.csv")

CV_summary <- CV_predictions %>%
  group_by(exclude_2801, model) %>%
  summarise(
    n = n(),
    grouped_CV_RMSE_mg_L = sqrt(mean((DO - predicted_DO)^2)),
    grouped_CV_MAE_mg_L = mean(abs(DO - predicted_DO)),
    .groups = "drop"
  )
write_csv(CV_summary, "outputs/infiltration_model_CV_summary.csv")

# ======================================================================
# FIGURES
# ======================================================================
theme_set(theme_bw(base_size = cfg$font_size))

p1 <- ggplot(profile, aes(
  river_recharge_DO_deficit_mg_L, grouped_CV_RMSE_mg_L,
  colour = factor(exclude_2801)
)) +
  geom_line(linewidth = .9) +
  labs(
    title = "Sensitivity to river-derived recharge DO deficit",
    x = "Assumed river-source DO deficit (mg/L)",
    y = "Grouped five-fold CV RMSE (mg/L)",
    colour = "Exclude 2801",
    subtitle = "Effective rainfall DO and apparent k refitted at each loss"
  )
print(p1)
ggsave("plots/infiltration_profile_CV.png", p1,
       width = 9, height = 6, dpi = cfg$plot_dpi)

p2 <- profile %>%
  filter(!exclude_2801) %>%
  ggplot(aes(river_recharge_DO_deficit_mg_L, k_per_year)) +
  geom_line(linewidth = 1) +
  labs(
    title = "River-source DO deficit versus apparent k",
    x = "Assumed river-source DO deficit (mg/L)",
    y = expression(k~(year^{-1}))
  )
print(p2)
ggsave("plots/infiltration_profile_k.png", p2,
       width = 9, height = 6, dpi = cfg$plot_dpi)

p3 <- required %>%
  filter(abs(river_reference_DO_mg_L - cfg$river_DO) < 1e-10,
         river_dominated, recharge_entry_only_feasible) %>%
  ggplot(aes(
    MRT, 100 * min_river_source_DO_deficit_fraction,
    shape = factor(confined)
  )) +
  geom_point(size = 3) +
  labs(
    title = "Minimum river-source DO deficit if k = 0",
    subtitle = paste(
      "River fraction >= 0.8; rainfall component may be fully anoxic;",
      "site 2801 omitted because its DO exceeds the fixed river reference"
    ),
    x = "Mean residence time (years)",
    y = "Minimum river-source DO deficit (%)",
    shape = "Confined"
  )
print(p3)
ggsave("plots/infiltration_minimum_loss_by_age.png", p3,
       width = 9, height = 6, dpi = cfg$plot_dpi)

cat("\nRiver-recharge DO deficit profile at whole mg/L increments:\n")
print(profile %>%
  filter(abs(river_recharge_DO_deficit_mg_L - round(river_recharge_DO_deficit_mg_L)) < 1e-9)
)
cat("\nModel comparison (grouped 5-fold CV):\n")
print(CV_summary)
cat("\nRequired river-source DO deficits for older river-dominated wells:\n")
print(summary_required)

cat(
  "\nINTERPRETATION:",
  "\n* Lower-bound river-source deficit assumes land-surface recharge is",
  "\n  completely anoxic, minimising the deficit from the river component.",
  "\n* Large older-well losses are mathematically necessary for k=0 under",
  "\n  the stated boundaries, but not proof of aquifer consumption.",
  "\n* Confinement is only a coarse proxy for recharge pathway conditions.",
  "\n* If initial recharge DO is allowed to vary with MRT, the k=0 model",
  "\n  can EXACTLY mimic first-order age-dependent depletion.",
  "\n* The data do not partition DO depletion along recharge pathways versus",
  "\n  during subsequent regional-aquifer transport.\n"
)
