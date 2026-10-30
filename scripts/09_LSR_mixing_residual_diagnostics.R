# 09_LSR_mixing_residual_diagnostics.R
#
# Diagnose whether the scatter around the first-order DO model is
# associated with (a) LSR contribution, (b) mixing-tracer disagreement,
# (c) uncertainty in mixing-tracer choice, or (d) spatial clustering.
#
# This is NOT another fitted oxygen-loss process.
# River DO reference is 7.5 mg/L; effective LSR DO is fitted in [0,9.5].
# "LSR" means LAND-SURFACE recharge; NOT riverbed recharge.
# The first-order k is an apparent age-associated DO depletion coefficient.
#
# Run after scripts/01_three_tracer_mixing.R. The optional spatial test
# requires private coordinates, either in data/groundwater.csv or in a
# local data/Uwe_HP_data_prepared.csv; NO coordinates are committed.
library(tidyverse)
setwd(here::here("Oxygen-Reduction-Joint-Model"))

# =====================================================================
# USER SETTINGS
# =====================================================================
cfg <- list(
  river_DO = 7.5,
  maximum_LSR_entry_DO = 9.5,
  maximum_k = 0.5,
  optional_original_data = "data/Uwe_HP_data_prepared.csv",
  spatial_neighbors = 5L,
  spatial_permutations = 9999L,
  seed = 20261030,
  figure_width = 10,
  figure_dpi = 300,
  font_size = 12
)

dir.create("outputs", showWarnings = FALSE, recursive = TRUE)
dir.create("plots", showWarnings = FALSE, recursive = TRUE)
file <- "outputs/three_tracer_mixing.csv"
if (!file.exists(file)) stop("Run scripts/01_three_tracer_mixing.R first.")

raw <- read_csv(file, show_col_types = FALSE,
                col_types = cols(Site_ID = col_character(),
                                 .default = col_guess()))
d <- raw %>%
  filter(
    mixing_reliable,
    is.finite(DO), DO > 0, is.finite(MRT), MRT >= 0,
    is.finite(rainfall_contribution),
    is.finite(river_contribution)
  ) %>%
  mutate(
    fLSR = rainfall_contribution,
    fRiver = river_contribution,
    tracer_fraction_range = mixing_disagreement
  )
stopifnot(
  nrow(d) == 29L, n_distinct(d$Site_ID) == 29L,
  sum(d$fLSR >= 0.5) == 5L,
  all(abs(d$fRiver + d$fLSR - 1) < 1e-8)
)

# =====================================================================
# BASELINE JOINT FIRST-ORDER FIT, WITH EFFECTIVE LSR DO PROFILED OUT.
# Exactly the model of script 02; no LSR DO spatial variability assumed.
# =====================================================================
predict_DO <- function(dat, P, k, fraction_col) {
  f <- dat[[fraction_col]]
  ((1 - f) * cfg$river_DO + f * P) * exp(-k * dat$MRT)
}

fit_first_order <- function(dat, fraction_col = "fLSR") {
  f <- dat[[fraction_col]]
  t <- dat$MRT
  y <- dat$DO

  profile <- function(k) {
    att <- exp(-k * t)
    river <- (1 - f) * cfg$river_DO * att
    w <- f * att
    P <- if (sum(w^2) > 0) {
      sum(w * (y - river)) / sum(w^2)
    } else {
      0
    }
    P <- pmax(0, pmin(cfg$maximum_LSR_entry_DO, P))
    pred <- river + w * P
    c(P = P, k = k, SSE = sum((y - pred)^2))
  }

  optimised <- optimize(
    function(k) profile(k)[["SSE"]],
    interval = c(0, cfg$maximum_k),
    tol = 1e-9
  )
  candidates <- lapply(
    c(0, optimised$minimum, cfg$maximum_k),
    profile
  )
  best <- candidates[[which.min(vapply(
    candidates, function(v) v[["SSE"]], numeric(1)
  ))]]
  list(
    P = unname(best[["P"]]),
    k = unname(best[["k"]]),
    SSE = unname(best[["SSE"]])
  )
}

LOO_predictions <- function(dat, fraction_col) {
  vapply(seq_len(nrow(dat)), function(i) {
    fit <- fit_first_order(dat[-i, , drop = FALSE], fraction_col)
    predict_DO(dat[i, , drop = FALSE], fit$P, fit$k,
               fraction_col)[[1]]
  }, numeric(1))
}

fit <- fit_first_order(d)
d <- d %>%
  mutate(
    model_prediction = predict_DO(d, fit$P, fit$k, "fLSR"),
    LOOWO_prediction = LOO_predictions(d, "fLSR"),
    LOOWO_residual = DO - LOOWO_prediction,
    absolute_LOOWO_residual = abs(LOOWO_residual)
  )

# =====================================================================
# DIRECT PHYSICAL SENSITIVITY OF DO PREDICTION TO LSR ENTRY DO
# =====================================================================
# Hold river DO and first-order k at the fitted values, then give EACH
# well completely independent effective LSR DO anywhere in [0, 9.5].
# This is deliberately more permissive than hierarchical LSR loss.
#
# For well i:
#   DOmin = (1-fLSR)*R*exp(-k*MRT)
#   DOmax = ((1-fLSR)*R + fLSR*9.5)*exp(-k*MRT)
# Even perfect knowledge of a well-specific LSR loss cannot reproduce
# DO outside these bounds unless another assumption changes.
atten <- exp(-fit$k * d$MRT)
d <- d %>%
  mutate(
    LSR_DO_prediction_min = fRiver * cfg$river_DO * atten,
    LSR_DO_prediction_max = (
      fRiver * cfg$river_DO +
        fLSR * cfg$maximum_LSR_entry_DO
    ) * atten,
    LSR_DO_flexibility_mg_L =
      fLSR * cfg$maximum_LSR_entry_DO * atten,
    LSR_variation_could_match_DO =
      DO >= LSR_DO_prediction_min - 1e-9 &
      DO <= LSR_DO_prediction_max + 1e-9,
    unexplained_after_any_LSR_DO_mg_L = pmax(
      0, LSR_DO_prediction_min - DO,
      DO - LSR_DO_prediction_max
    ),
    # An even more permissive MIXING test: allow fLSR to be ANY value
    # in [0,1], keeping the fitted common-source DO values unchanged.
    unconstrained_mixing_min =
      min(cfg$river_DO, fit$P) * atten,
    unconstrained_mixing_max =
      max(cfg$river_DO, fit$P) * atten,
    any_mixing_fraction_could_match_DO =
      DO >= unconstrained_mixing_min - 1e-9 &
      DO <= unconstrained_mixing_max + 1e-9
  )

# At fixed fitted model parameters, examine the physical impact of
# switching among the three INDIVIDUAL tracer-derived LSR fractions.
tracers <- c(
  "rain_fraction_18O",
  "rain_fraction_Cl",
  "rain_fraction_excess_air"
)
stopifnot(all(tracers %in% names(d)))
individual_predictions <- lapply(tracers, function(col) {
  predict_DO(d, fit$P, fit$k, col)
})
d$tracer_choice_DO_range_mg_L <-
  do.call(pmax, individual_predictions) -
  do.call(pmin, individual_predictions)

# =====================================================================
# RESIDUAL ASSOCIATIONS: NO CORRELATION ASSUMED A PRIORI.
# Screened tracer disagreement is truncated to <=0.35 by script 01.
# Failure to find a correlation DOES NOT rule out systematic mixing error.
# =====================================================================
correlation_test <- function(dat, variable, signed, subset_name) {
  response <- if (signed) {
    dat$LOOWO_residual
  } else {
    dat$absolute_LOOWO_residual
  }
  test <- suppressWarnings(cor.test(
    dat[[variable]], response,
    method = "spearman", exact = FALSE
  ))
  tibble(
    subset = subset_name, n = nrow(dat), variable = variable,
    residual = if (signed) "Signed" else "Absolute",
    Spearman_rho = unname(test$estimate),
    p_value_exploratory = test$p.value
  )
}

correlations <- map_dfr(
  c("All 29", "Excluding 2801"),
  function(selection) {
    sub <- if (selection == "All 29") {
      d
    } else {
      d %>% filter(Site_ID != "2801")
    }
    map_dfr(c("fLSR", "tracer_fraction_range"), function(var) {
      bind_rows(
        correlation_test(sub, var, FALSE, selection),
        correlation_test(sub, var, TRUE, selection)
      )
    })
  }
)

# Same 29 wells, alternative tracer fractions: refit k and effective LSR
# DO, and recompute LOOWO predictions. Do NOT change selection here.
alternate_tracers <- map_dfr(
  c("fLSR", tracers), function(col) {
    alt <- fit_first_order(d, col)
    cv <- LOO_predictions(d, col)
    tibble(
      LSR_fraction_method = col,
      n = nrow(d),
      effective_LSR_DO = alt$P,
      apparent_k_per_year = alt$k,
      fit_RMSE_mg_L = sqrt(alt$SSE / nrow(d)),
      LOOWO_RMSE_mg_L = sqrt(mean((d$DO - cv)^2))
    )
  }
)

# =====================================================================
# OPTIONAL SPATIAL AUTOCORRELATION
# =====================================================================
# Optionally join coordinates from either locally available table.
# The original prepared CSV can include same Site_ID at multiple depths.
# Join using both Site_ID and Depth when available.
spatial_result <- tibble()
loc <- NULL

if (all(c("E", "N") %in% names(d))) {
  loc <- d
} else if (file.exists(cfg$optional_original_data)) {
  original <- read_csv(
    cfg$optional_original_data, show_col_types = FALSE,
    col_types = cols(Site_ID = col_character(),
                     .default = col_guess())
  )
  if (all(c("Site_ID", "E", "N") %in% names(original))) {
    keys <- if ("Depth" %in% names(original) &&
                "Depth" %in% names(d)) {
      c("Site_ID", "Depth")
    } else {
      "Site_ID"
    }
    locations <- original %>%
      select(all_of(keys), E, N) %>%
      distinct(across(all_of(keys)), .keep_all = TRUE)
    loc <- left_join(d, locations, by = keys)
  }
}

if (!is.null(loc) &&
    all(c("E", "N") %in% names(loc))) {
  loc <- loc %>%
    filter(is.finite(E), is.finite(N)) %>%
    mutate(E = as.numeric(E), N = as.numeric(N))
  n <- nrow(loc)
  if (n > cfg$spatial_neighbors + 2L) {
    XY <- as.matrix(loc[, c("E", "N")])
    dist_matrix <- as.matrix(dist(XY))
    diag(dist_matrix) <- Inf
    W <- matrix(0, nrow = n, ncol = n)
    for (i in seq_len(n)) {
      idx <- order(dist_matrix[i, ])[seq_len(cfg$spatial_neighbors)]
      W[i, idx] <- 1
    }
    W <- 1 * ((W + t(W)) > 0)
    W <- W / rowSums(W)
    centered <- loc$LOOWO_residual - mean(loc$LOOWO_residual)
    moran <- function(x) {
      as.numeric(crossprod(x, W %*% x) / crossprod(x, x))
    }
    observed <- moran(centered)
    set.seed(cfg$seed)
    null <- replicate(
      cfg$spatial_permutations,
      moran(sample(centered))
    )
    spatial_result <- tibble(
      n = n,
      neighbors = cfg$spatial_neighbors,
      Moran_I = observed,
      permutation_p_positive = (
        1 + sum(null >= observed)
      ) / (cfg$spatial_permutations + 1),
      permutation_p_two_sided = (
        1 + sum(abs(null - mean(null)) >=
                  abs(observed - mean(null)))
      ) / (cfg$spatial_permutations + 1)
    )
  }
} else {
  message(
    "Spatial Moran test skipped: coordinates absent. To enable,",
    " provide E and N columns in the local groundwater CSV or place ",
    "a local original file at ", cfg$optional_original_data
  )
}

# =====================================================================
# DIAGNOSTICS, TABLES AND FIGURES
# =====================================================================
write_csv(d, "outputs/LSR_mixing_residuals_by_well.csv")
write_csv(correlations, "outputs/LSR_mixing_residual_correlations.csv")
write_csv(alternate_tracers,
          "outputs/LSR_mixing_alternative_tracers.csv")
if (nrow(spatial_result)) {
  write_csv(spatial_result,
            "outputs/LSR_mixing_residual_spatial_test.csv")
}

theme_set(theme_bw(base_size = cfg$font_size))
top_errors <- d %>% slice_max(absolute_LOOWO_residual, n = 5)

p1 <- ggplot(d, aes(fLSR, absolute_LOOWO_residual)) +
  geom_point(size = 2.8) +
  geom_text(
    data = top_errors, aes(label = Site_ID),
    hjust = -0.08, vjust = -0.4, size = 3,
    check_overlap = TRUE
  ) +
  scale_x_continuous(labels = scales::percent,
                     limits = c(0, 1)) +
  labs(
    title = "Does greater LSR contribution mean larger DO residuals?",
    subtitle = "Leave-one-well-out first-order DO residuals",
    x = "Land-surface recharge fraction",
    y = "Absolute DO residual (mg/L)"
  )
print(p1)
ggsave("plots/LSR_mixing_abs_residual_vs_LSR.png", p1,
       width = cfg$figure_width, height = 6, dpi = cfg$figure_dpi)

p2 <- ggplot(
  d, aes(tracer_fraction_range, absolute_LOOWO_residual)
) +
  geom_point(size = 2.8) +
  geom_text(
    data = top_errors, aes(label = Site_ID),
    hjust = -0.08, vjust = -0.4, size = 3,
    check_overlap = TRUE
  ) +
  labs(
    title = "Do uncertain mixing fractions produce larger DO residuals?",
    subtitle = "Range of three UNCLIPPED tracer-derived LSR fractions",
    x = "Intertracer fraction disagreement",
    y = "Absolute leave-one-well-out DO residual (mg/L)"
  )
print(p2)
ggsave("plots/LSR_mixing_abs_residual_vs_tracer_disagreement.png",
       p2, width = cfg$figure_width, height = 6,
       dpi = cfg$figure_dpi)

p3 <- d %>%
  mutate(
    Site_ID_plot = forcats::fct_reorder(Site_ID, fLSR)
  ) %>%
  ggplot(aes(y = Site_ID_plot)) +
  geom_segment(
    aes(x = LSR_DO_prediction_min,
        xend = LSR_DO_prediction_max,
        yend = Site_ID_plot),
    linewidth = 1.5
  ) +
  geom_point(aes(x = DO, shape = "Measured DO"), size = 2.6) +
  geom_point(aes(x = model_prediction,
                 shape = "Common-rate prediction"),
             size = 2.2) +
  labs(
    title = "Can any LSR DO between 0 and 9.5 mg/L fit each well?",
    subtitle = "Horizontal interval: physically available LSR DO variation; river DO and k fixed",
    x = "Groundwater DO (mg/L)",
    y = "Well (ordered by LSR contribution)",
    shape = NULL
  ) +
  theme(legend.position = "bottom")
print(p3)
ggsave("plots/LSR_mixing_possible_DO_by_well.png",
       p3, width = 11, height = 10, dpi = cfg$figure_dpi)

cat("\nBaseline first-order P,k:\n")
print(tibble(effective_LSR_DO = fit$P, apparent_k = fit$k,
             fitted_RMSE = sqrt(fit$SSE / nrow(d)),
             LOOWO_RMSE = sqrt(mean(d$LOOWO_residual^2))))
cat("\nResidual correlations (exploratory, unadjusted p-values):\n")
print(correlations)
cat("\nAlternative tracer predictions (same 29 wells):\n")
print(alternate_tracers)
cat("\nPhysical flexibility of mixing and LSR DO:\n")
print(tibble(
  n = nrow(d),
  n_fit_if_LSR_DO_free = sum(d$LSR_variation_could_match_DO),
  n_fit_if_mixing_free = sum(d$any_mixing_fraction_could_match_DO),
  median_LSR_DO_prediction_span =
    median(d$LSR_DO_flexibility_mg_L),
  median_tracer_DO_prediction_range =
    median(d$tracer_choice_DO_range_mg_L),
  max_tracer_DO_prediction_range =
    max(d$tracer_choice_DO_range_mg_L)
))
if (nrow(spatial_result)) {
  cat("\nSpatial diagnostic (optional coordinates):\n")
  print(spatial_result)
}
cat(
  "\nIMPORTANT:\n",
  "* Lack of correlation does NOT exclude LSR DO spatial variability.\n",
  "* Selection on tracer agreement restricts the range of disagreement.\n",
  "* Tracer disagreement does NOT detect shared systematic bias.\n",
  "* LSR flexibility bounds fix river DO and apparent k; a different\n",
  "  k or groundwater age distribution changes the prediction range.\n",
  "* Spatial autocorrelation is low power with 29, clustered wells.\n",
  "* DO–MRT relationship alone does not locate the oxygen reduction.\n",
  sep = ""
)
