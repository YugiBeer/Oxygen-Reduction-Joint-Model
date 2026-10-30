# 05_recharge_pathway_confounding.R
#
# Sensitivity analysis: Could correlated recharge conditions explain the
# observed DO–groundwater-age relationship, without aquifer reduction?
#
# Four Bayesian models are compared with identical source-DO priors:
#   A) Recharge only: one mean initial DO per source; no age term
#   B) Pathway proxy only: confined/unconfined source DO means; k=0
#   C) Age effect only: common source DO means; fit k
#   D) Both age effect and pathway proxy
#
# A/B do not imply all loss occurs at recharge. The "pathway" model is
# a sensitivity experiment: confinement is measured at the WELL and is
# not necessarily a proxy for where oxygen was consumed.
#
# IDENTIFIABILITY: If source DO varies as R_i=R*exp(-b*t_i) and
# P_i=P*exp(-b*t_i), this predicts exactly the same DO as a
# first-order model with k=b. No DO/MRT model fit can tell these apart.
#
# This script reuses outputs/three_tracer_mixing.csv from script 01.
# It does NOT modify existing outputs from scripts 02 or 04.
library(tidyverse)
setwd(here::here("Oxygen-Reduction-Joint-Model"))

# ----------------------------------------------------------------------
# User settings / priors
# ----------------------------------------------------------------------
cfg <- list(
  river_mean = 7.5,            # mg/L, river-monitoring mean
  river_mean_sd = 0.7,         # assumed uncertainty of the MEAN
  river_between_sd = 1.9,      # half-normal scale, varying inlet DO
  rainfall_upper = 9.5,       # uniform prior on effective rain DO mean
  rainfall_between_sd = 2.5,   # half-normal scale
  unexplained_sd = 1.5,       # half-normal residual SD scale
  measurement_sd = 0.2,       # mg/L, assumed
  group_offset_sd = 3.0,      # normal prior for confined minus unconfined
  k_sd = 0.10,                # half-normal scale, yr^-1
  n_full_chains = 4L,
  full_warmup = 1400L,
  full_draws = 1800L,         # samples retained per chain
  n_cv_chains = 2L,
  cv_warmup = 1000L,
  cv_draws = 900L,
  thin = 2L,
  n_cv_folds = 5L,
  seed = 20261030L,
  plot_font_size = 12,
  plot_dpi = 300
)

dir.create("outputs", showWarnings = FALSE, recursive = TRUE)
dir.create("plots", showWarnings = FALSE, recursive = TRUE)
input <- "outputs/three_tracer_mixing.csv"
if (!file.exists(input)) stop("Run script 01 first.")
d <- read_csv(input, show_col_types = FALSE,
              col_types = cols(Site_ID = col_character(),
                               .default = col_guess())) %>%
  filter(mixing_reliable, is.finite(DO), DO > 0,
         is.finite(MRT), MRT >= 0,
         is.finite(rainfall_contribution),
         is.finite(river_contribution),
         !is.na(confined)) %>%
  transmute(Site_ID, DO, MRT, Depth, confined = as.integer(confined),
            fR = river_contribution, fP = rainfall_contribution) %>%
  mutate(
    cv_group = case_when(
      str_detect(Site_ID, "Eastbourne") ~ "Eastbourne",
      str_detect(Site_ID, "Gilbertson") ~ "Gilbertson",
      str_detect(Site_ID, "Flaxmere") ~ "Flaxmere",
      TRUE ~ Site_ID
    )
  )
stopifnot(nrow(d) == 29L, n_distinct(d$Site_ID) == 29L,
          all(d$confined %in% c(0L, 1L)),
          all(abs(d$fR + d$fP - 1) < 1e-8))

model_specs <- tribble(
  ~model, ~pathway, ~age_term,
  "Recharge only", FALSE, FALSE,
  "Confinement proxy only", TRUE, FALSE,
  "Age effect only", FALSE, TRUE,
  "Confinement + age", TRUE, TRUE
)

# ----------------------------------------------------------------------
# Hierarchical model: well-specific river and rain DO are Gaussian
# and integrated out analytically. Physically bounded source DO could
# be investigated in a later analysis.
# ----------------------------------------------------------------------
unpack <- function(z, pathway, age_term) {
  idx <- 6L
  if (pathway) {
    dR <- z[idx]
    dP <- z[idx + 1L]
    idx <- idx + 2L
  } else {
    dR <- 0
    dP <- 0
  }
  k <- if (age_term) exp(z[idx]) else 0
  c(
    mu_R = z[1], mu_P = z[2],
    sigma_R = exp(z[3]), sigma_P = exp(z[4]),
    sigma_extra = exp(z[5]),
    delta_R = dR, delta_P = dP, k = k
  )
}

log_posterior <- function(z, dat, pathway, age_term) {
  p <- unpack(z, pathway, age_term)
  if (any(!is.finite(p)) ||
      p["mu_R"] <= 0 || p["mu_R"] >= 12 ||
      p["mu_P"] <= 0 || p["mu_P"] >= cfg$rainfall_upper ||
      any(p[c("sigma_R", "sigma_P", "sigma_extra")] >= 30) ||
      p["k"] >= 3) return(-Inf)

  if (pathway && (
    p["mu_R"] + p["delta_R"] <= 0 ||
    p["mu_R"] + p["delta_R"] >= 12 ||
    p["mu_P"] + p["delta_P"] <= 0 ||
    p["mu_P"] + p["delta_P"] >= cfg$rainfall_upper
  )) return(-Inf)

  # Half-normal priors: log(2) constants cancel from within-model MCMC.
  lp <- dnorm(p["mu_R"], cfg$river_mean, cfg$river_mean_sd,
              log = TRUE) +
    dnorm(p["sigma_R"], 0, cfg$river_between_sd, log = TRUE) + z[3] +
    dnorm(p["sigma_P"], 0, cfg$rainfall_between_sd, log = TRUE) + z[4] +
    dnorm(p["sigma_extra"], 0, cfg$unexplained_sd, log = TRUE) + z[5]

  if (pathway) {
    lp <- lp +
      dnorm(p["delta_R"], 0, cfg$group_offset_sd, log = TRUE) +
      dnorm(p["delta_P"], 0, cfg$group_offset_sd, log = TRUE)
  }
  if (age_term) {
    lp <- lp + dnorm(p["k"], 0, cfg$k_sd, log = TRUE) +
      z[length(z)]
  }

  att <- exp(-p["k"] * dat$MRT)
  mu <- (
    dat$fR * (p["mu_R"] + dat$confined * p["delta_R"]) +
    dat$fP * (p["mu_P"] + dat$confined * p["delta_P"])
  ) * att

  variance <- (
    dat$fR^2 * p["sigma_R"]^2 +
    dat$fP^2 * p["sigma_P"]^2
  ) * att^2 + p["sigma_extra"]^2 + cfg$measurement_sd^2

  lp + sum(dnorm(dat$DO, mu, sqrt(variance), log = TRUE))
}

# Componentwise Metropolis-Hastings. Tune ONLY in warmup.
sample_chain <- function(dat, pathway, age_term, chain_seed,
                         warmup, n_draws, thin) {
  set.seed(chain_seed)
  z <- c(cfg$river_mean, 6, log(cfg$river_between_sd),
         log(cfg$rainfall_between_sd), log(1))
  if (pathway) z <- c(z, -3, -2)
  if (age_term) z <- c(z, log(0.07))

  steps <- c(0.26, 0.45, 0.20, 0.20, 0.20)
  if (pathway) steps <- c(steps, 0.38, 0.38)
  if (age_term) steps <- c(steps, 0.20)

  current <- log_posterior(z, dat, pathway, age_term)
  stopifnot(is.finite(current))
  accepted <- integer(length(z))
  keep <- matrix(NA_real_, nrow = n_draws, ncol = 8)
  colnames(keep) <- c("mu_R", "mu_P", "sigma_R", "sigma_P",
                      "sigma_extra", "delta_R", "delta_P", "k")
  counter <- 0L

  for (iteration in seq_len(warmup + n_draws * thin)) {
    for (j in sample.int(length(z))) {
      trial <- z
      trial[j] <- trial[j] + rnorm(1, sd = steps[j])
      new_lp <- log_posterior(trial, dat, pathway, age_term)
      if (is.finite(new_lp) &&
          log(runif(1)) < new_lp - current) {
        z <- trial
        current <- new_lp
        accepted[j] <- accepted[j] + 1L
      }
    }
    if (iteration <= warmup && iteration %% 100L == 0L) {
      rates <- accepted / 100
      steps <- pmax(0.015, pmin(2, steps *
        exp(pmax(-0.25, pmin(0.25, rates - 0.34)))))
      accepted[] <- 0L
    }
    if (iteration > warmup &&
        (iteration - warmup) %% thin == 0L) {
      counter <- counter + 1L
      keep[counter, ] <- unpack(z, pathway, age_term)
    }
  }
  as_tibble(keep)
}

sample_model <- function(dat, pathway, age_term,
                         n_chains, warmup, n_draws, seed) {
  map_dfr(seq_len(n_chains), function(ch) {
    sample_chain(dat, pathway, age_term,
                 chain_seed = seed + 1009L * ch,
                 warmup = warmup, n_draws = n_draws,
                 thin = cfg$thin) %>%
      mutate(chain = ch)
  })
}

predict_draws <- function(draws, dat) {
  n <- nrow(draws)
  m <- nrow(dat)
  att <- exp(-tcrossprod(draws$k, dat$MRT))
  # Each group delta is applied to population source means.
  initial <- tcrossprod(draws$mu_R, dat$fR) +
    tcrossprod(draws$mu_P, dat$fP) +
    tcrossprod(draws$delta_R, dat$fR * dat$confined) +
    tcrossprod(draws$delta_P, dat$fP * dat$confined)
  expected <- initial * att

  variance <- (
    tcrossprod(draws$sigma_R^2, dat$fR^2) +
    tcrossprod(draws$sigma_P^2, dat$fP^2)
  ) * att^2 +
    matrix(draws$sigma_extra^2 + cfg$measurement_sd^2,
           nrow = n, ncol = m)
  list(mean = expected, variance = variance)
}

posterior_summary <- function(draws, model) {
  draws %>%
    select(mu_R, mu_P, sigma_R, sigma_P, sigma_extra,
           delta_R, delta_P, k) %>%
    pivot_longer(cols = everything(), names_to = "parameter") %>%
    group_by(parameter) %>%
    summarise(
      median = median(value),
      q025 = unname(quantile(value, 0.025)),
      q975 = unname(quantile(value, 0.975)),
      .groups = "drop"
    ) %>%
    mutate(model = model, .before = 1)
}

# ----------------------------------------------------------------------
# Full posterior estimates
# ----------------------------------------------------------------------
posterior <- map(model_specs$model, function(model) {
  spec <- model_specs %>% filter(.data$model == model)
  cat("Posterior:", model, "\n")
  sample_model(
    d, spec$pathway, spec$age_term,
    cfg$n_full_chains, cfg$full_warmup, cfg$full_draws,
    cfg$seed + 50000L * match(model, model_specs$model)
  )
})
names(posterior) <- model_specs$model

summaries <- map2_dfr(posterior, names(posterior), posterior_summary)
write_csv(summaries, "outputs/confounding_posterior_summary.csv")

# ----------------------------------------------------------------------
# Grouped five-fold CV: keep obvious multi-depth well groups together.
# The grouping is HEURISTIC and should not be interpreted as spatial CV.
# ----------------------------------------------------------------------
set.seed(cfg$seed)
groups <- split(seq_len(nrow(d)), d$cv_group)
group_order <- sample(seq_along(groups))
group_order <- group_order[order(
  lengths(groups)[group_order], decreasing = TRUE
)]
folds <- vector("list", cfg$n_cv_folds)
for (g in group_order) {
  fold_idx <- which.min(lengths(folds))
  folds[[fold_idx]] <- c(folds[[fold_idx]], groups[[g]])
}
stopifnot(all(lengths(folds) > 0L),
          identical(sort(unlist(folds)), seq_len(nrow(d))))

cv <- map_dfr(seq_len(nrow(model_specs)), function(i) {
  spec <- model_specs[i, ]
  map_dfr(seq_along(folds), function(j) {
    held <- folds[[j]]
    train <- d[-held, , drop = FALSE]
    test <- d[held, , drop = FALSE]
    cat("Grouped CV:", spec$model, "fold", j, "of", length(folds), "\n")

    draws <- sample_model(
      train, spec$pathway, spec$age_term,
      cfg$n_cv_chains, cfg$cv_warmup, cfg$cv_draws,
      cfg$seed + 500000L + i * 10000L + j * 311L
    )
    pp <- predict_draws(draws, test)

    map_dfr(seq_len(nrow(test)), function(k) {
      logdensity <- dnorm(
        test$DO[k], mean = pp$mean[, k],
        sd = sqrt(pp$variance[, k]), log = TRUE
      )
      maximum <- max(logdensity)
      tibble(
        model = spec$model,
        fold = j,
        Site_ID = test$Site_ID[k],
        cv_group = test$cv_group[k],
        confined = test$confined[k],
        MRT = test$MRT[k],
        fP = test$fP[k],
        observed_DO = test$DO[k],
        predicted_DO = mean(pp$mean[, k]),
        log_predictive_density =
          maximum + log(mean(exp(logdensity - maximum)))
      )
    })
  })
})
write_csv(cv, "outputs/confounding_group_CV_predictions.csv")

performance <- cv %>%
  group_by(model) %>%
  summarise(
    n = n(),
    RMSE_mg_L = sqrt(mean((observed_DO - predicted_DO)^2)),
    MAE_mg_L = mean(abs(observed_DO - predicted_DO)),
    sum_log_predictive_density = sum(log_predictive_density),
    .groups = "drop"
  )
write_csv(performance, "outputs/confounding_group_CV_summary.csv")

# ----------------------------------------------------------------------
# Descriptive checks and figures
# ----------------------------------------------------------------------
correlations <- d %>%
  group_by(confined) %>%
  summarise(
    n = n(),
    spearman_DO_MRT = unname(cor(DO, MRT, method = "spearman")),
    mean_MRT = mean(MRT), mean_DO = mean(DO),
    .groups = "drop"
  )
write_csv(correlations, "outputs/confounding_within_group_correlations.csv")

theme_set(theme_bw(base_size = cfg$plot_font_size))
p <- ggplot(d, aes(MRT, DO, colour = factor(confined))) +
  geom_point(size = 3, alpha = 0.9) +
  facet_wrap(~ confined, labeller = as_labeller(
    c("0" = "Unconfined wells", "1" = "Confined wells")
  )) +
  labs(
    title = "DO and MRT within confinement groups",
    subtitle = "Confinement and MRT are correlated; neither defines recharge location",
    x = "Mean residence time (years)", y = "Measured DO (mg/L)",
    colour = "Confined"
  ) +
  guides(colour = "none")
print(p)
ggsave("plots/confounding_DO_MRT_by_confinement.png", p,
       width = 10, height = 5, dpi = cfg$plot_dpi)

p2 <- ggplot(cv, aes(observed_DO, predicted_DO)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_point(aes(colour = MRT), size = 2.7) +
  facet_wrap(~ model, ncol = 2) +
  coord_equal() +
  scale_colour_viridis_c() +
  labs(
    title = "Grouped cross-validation: recharge pathway versus age",
    subtitle = "Eastbourne / Gilbertson / Flaxmere groups held out together",
    x = "Observed DO (mg/L)", y = "Predicted DO (mg/L)",
    colour = "MRT (years)"
  )
print(p2)
ggsave("plots/confounding_group_CV_DO.png", p2,
       width = 11, height = 9, dpi = cfg$plot_dpi)

print(performance)
print(correlations)
cat(
  "\nCAUTION: A nonzero k supports age-associated DO variation only",
  "\nwhen initial recharge DO differences are assumed independent of age.",
  "\nAllowing initial recharge DO to fall with MRT exactly reproduces",
  "\nthe age-term formulation. Physical site of reduction is not identified.",
  "\nDO levels in confined wells are not measurements of their inlet DO.\n"
)
