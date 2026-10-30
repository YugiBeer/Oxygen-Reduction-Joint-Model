# 13_joint_posterior_partition_diagnostics.R
#
# Assess whether the joint Bayesian tracer-mixing model can identify
# LSR oxygen loss separately from apparent groundwater-age depletion.
#
# This reads locally saved posterior objects from script 12.
# It NEVER publishes per-well measurements or posterior outputs.
# The pre-recharge DO of LSR is assumed 9.5 mg/L, not measured.
#
# Important: "age-associated DO depletion" below is an apparent
# reaction term. This statistical decomposition does NOT establish
# that all of it occurs in the saturated aquifer.
#
# Run AFTER scripts/12_joint_bayesian_mixing.R has completed.
library(tidyverse)
setwd(here::here("Oxygen-Reduction-Joint-Model"))

if (!requireNamespace("posterior", quietly = TRUE)) {
  stop("Install posterior before running this script.")
}
source_file <- "data/groundwater.csv"
if (!file.exists(source_file)) {
  stop("Private data/groundwater.csv is required (never commit it).")
}

dir.create("outputs", showWarnings = FALSE, recursive = TRUE)
dir.create("plots", showWarnings = FALSE, recursive = TRUE)

d <- read_csv(
  source_file, show_col_types = FALSE,
  col_types = cols(Site_ID = col_character(), .default = col_guess())
) %>%
  filter(is.finite(DO), DO >= 0, is.finite(MRT), MRT >= 0) %>%
  mutate(sample_index = row_number())
n <- nrow(d)
stopifnot(n >= 20L)

files <- list.files(
  "outputs",
  pattern = "^joint_bayes_mixing_.*\\.rds$",
  full.names = TRUE
)
if (!length(files)) {
  stop("Run scripts/12_joint_bayesian_mixing.R first.")
}

scalar_fields <- c(
  "k", "river_DO", "central_LSR_entry_DO",
  "central_LSR_DO_loss", "sample_mean_LSR_DO_loss",
  "sigma_DO", "cl_R", "cl_LSR",
  "o18_R", "o18_LSR", "air_R", "air_LSR"
)

quantiles <- function(x) c(
  median = median(x),
  lower95 = unname(quantile(x, .025)),
  upper95 = unname(quantile(x, .975))
)
summ <- function(x, metric) {
  v <- quantiles(x)
  tibble(metric = metric, median = v["median"],
         lower95 = v["lower95"], upper95 = v["upper95"])
}

analyse_fit <- function(path) {
  scenario <- sub("^joint_bayes_mixing_", "", basename(path))
  scenario <- sub("\\.rds$", "", scenario)
  message("Posterior diagnostics: ", scenario)
  fit <- readRDS(path)

  vars <- c(
    scalar_fields,
    "f_LSR", "LSR_entry_DO", "mean_DO", "true_MRT"
  )
  M <- as.matrix(posterior::as_draws_matrix(
    fit$draws(variables = vars)
  ))

  f_cols <- paste0("f_LSR[", seq_len(n), "]")
  P_cols <- paste0("LSR_entry_DO[", seq_len(n), "]")
  DO_cols <- paste0("mean_DO[", seq_len(n), "]")
  MRT_cols <- paste0("true_MRT[", seq_len(n), "]")
  required <- c(
    "k", "river_DO", "central_LSR_DO_loss",
    f_cols, P_cols, DO_cols, MRT_cols
  )
  if (!all(required %in% colnames(M))) {
    stop("Missing variables from posterior: ",
         paste(setdiff(required, colnames(M)), collapse = ", "))
  }

  f <- M[, f_cols, drop = FALSE]
  P <- M[, P_cols, drop = FALSE]
  mean_DO <- M[, DO_cols, drop = FALSE]
  latent_age <- M[, MRT_cols, drop = FALSE]
  R <- M[, "river_DO"]

  # Each draw gives a full model-consistent oxygen budget FOR EACH WELL.
  # Before LSR infiltration: river contribution R and LSR reference 9.5.
  # At aquifer entry: river component R and effective LSR component P_i.
  # At sampling: modelled mean_DO after the apparent age-associated term.
  LSR_loss <- f * (9.5 - P)
  entry_DO <- (1-f) * R + f * P
  age_loss <- entry_DO - mean_DO
  total_loss <- LSR_loss + age_loss

  # These are arithmetic summaries across wells, NOT flow-volume-weighted
  # recharge budgets. Wells are not a random sample of the regional aquifer.
  mean_lsr_loss <- rowMeans(LSR_loss)
  mean_age_loss <- rowMeans(age_loss)
  mean_total_loss <- rowMeans(total_loss)
  lsr_fraction_of_loss <- if_else(
    mean_total_loss > 0,
    mean_lsr_loss / mean_total_loss,
    NA_real_
  )

  d_summary <- bind_rows(
    summ(mean_lsr_loss, "Mean modelled LSR pathway DO loss, per well"),
    summ(mean_age_loss, "Mean modelled age-associated DO loss, per well"),
    summ(lsr_fraction_of_loss, "LSR share of modelled per-well DO loss")
  ) %>%
    mutate(scenario = scenario, .before = 1)

  draws <- tibble(
    scenario = scenario,
    k = M[, "k"],
    LSR_loss_central = M[, "central_LSR_DO_loss"],
    river_DO = R,
    mean_fLSR = rowMeans(f),
    LSR_loss_mean_well = mean_lsr_loss,
    age_loss_mean_well = mean_age_loss,
    fraction_LSR_of_loss = lsr_fraction_of_loss
  )

  corr_pairs <- tribble(
    ~x, ~y,
    "k", "LSR_loss_central",
    "k", "river_DO",
    "k", "mean_fLSR",
    "LSR_loss_central", "mean_fLSR"
  )
  corrs <- pmap_dfr(corr_pairs, function(x, y) {
    tibble(
      scenario = scenario,
      parameter_x = x, parameter_y = y,
      posterior_pearson = cor(draws[[x]], draws[[y]]),
      posterior_spearman = cor(
        draws[[x]], draws[[y]], method = "spearman"
      )
    )
  })

  # Detailed per-well summaries remain PRIVATE in outputs/.
  per_well <- map_dfr(seq_len(n), function(i) {
    bind_rows(
      summ(f[, i], "LSR fraction"),
      summ(P[, i], "Effective LSR entry DO"),
      summ(latent_age[, i], "Latent MRT"),
      summ(LSR_loss[, i], "LSR pathway loss contribution"),
      summ(age_loss[, i], "Apparent age-associated loss"),
      summ(mean_DO[, i], "Modelled sample DO")
    ) %>%
      mutate(
        scenario = scenario,
        sample_index = i,
        Site_ID = d$Site_ID[i], Depth = d$Depth[i],
        MRT_reported = d$MRT[i], DO_observed = d$DO[i],
        .before = 1
      )
  })

  diag <- fit$summary(variables = scalar_fields) %>%
    select(variable, rhat, ess_bulk, ess_tail) %>%
    mutate(scenario = scenario, .before = 1)

  list(draws = draws, correlations = corrs,
       budget = d_summary, per_well = per_well,
       convergence = diag)
}

results <- map(files, analyse_fit)
names(results) <- basename(files)

draws <- bind_rows(map(results, "draws"))
correlations <- bind_rows(map(results, "correlations"))
budget <- bind_rows(map(results, "budget"))
wells <- bind_rows(map(results, "per_well"))
convergence <- bind_rows(map(results, "convergence"))

write_csv(correlations,
          "outputs/joint_bayes_partition_correlations.csv")
write_csv(budget,
          "outputs/joint_bayes_partition_summary.csv")
write_csv(wells,
          "outputs/joint_bayes_partition_by_well.csv")
write_csv(convergence,
          "outputs/joint_bayes_partition_convergence.csv")

# Reproducible graphics. The baseline is shown with alternative priors
# so any trade-off is evident without interpreting medians in isolation.
p1 <- ggplot(
  draws,
  aes(LSR_loss_central, k, colour = scenario)
) +
  geom_point(alpha = .12, size = .45) +
  facet_wrap(~ scenario) +
  guides(colour = "none") +
  labs(
    title = "Joint posterior: LSR oxygen deficit versus apparent k",
    subtitle = "Negative correlation indicates partial confounding",
    x = "Central effective LSR DO deficit (mg/L)",
    y = expression(k~(year^{-1}))
  ) +
  theme_bw(base_size = 12)
print(p1)
ggsave(
  "plots/joint_bayes_partition_LSR_vs_k.png",
  p1, width = 11, height = 7, dpi = 250
)

p2 <- ggplot(
  draws,
  aes(river_DO, k, colour = scenario)
) +
  geom_point(alpha = .12, size = .45) +
  facet_wrap(~ scenario) +
  guides(colour = "none") +
  labs(
    title = "Joint posterior: river DO reference versus apparent k",
    subtitle = "Higher initial river oxygen requires more age-associated loss",
    x = "River DO reference (mg/L)",
    y = expression(k~(year^{-1}))
  ) +
  theme_bw(base_size = 12)
print(p2)
ggsave(
  "plots/joint_bayes_partition_riverDO_vs_k.png",
  p2, width = 11, height = 7, dpi = 250
)

p3 <- ggplot(
  budget,
  aes(median, scenario)
) +
  geom_segment(aes(x = lower95, xend = upper95,
                   y = scenario, yend = scenario), linewidth = .8) +
  geom_point(size = 2.5) +
  facet_wrap(~ metric, scales = "free_x") +
  labs(
    title = "Prior sensitivity of model-implied per-well oxygen depletion",
    subtitle = "Unweighted well summaries: NOT an aquifer or catchment flux budget",
    x = "Posterior median and 95% credible interval", y = "Scenario"
  ) +
  theme_bw(base_size = 11)
print(p3)
ggsave(
  "plots/joint_bayes_partition_loss_components.png",
  p3, width = 12, height = 6, dpi = 250
)

cat("\nPOSTERIOR PARAMETER CORRELATIONS:\n")
print(correlations, n = Inf)
cat("\nMODEL-IMPLIED LOSS DECOMPOSITION (UNWEIGHTED WELLS):\n")
print(budget, n = Inf)
cat("\nCONVERGENCE:\n")
print(convergence, n = Inf)
if (any(convergence$rhat > 1.05, na.rm = TRUE)) {
  warning("Some posterior parameters have Rhat > 1.05.")
}
if (any(convergence$ess_bulk < 400, na.rm = TRUE)) {
  warning("Some posterior parameters have bulk ESS < 400.")
}
cat(
  "\nINTERPRETATION:\n",
  "1. The LSR reference 9.5 mg/L is assumed, not measured.\n",
  "2. The LSR pathway component is fLSR*(9.5 - LSR_entry_DO).\n",
  "3. The age-associated component is entry_DO - modelled_sample_DO.\n",
  "4. Age-associated loss is NOT proof of saturated-aquifer reaction.\n",
  "5. Arithmetic averages across wells are NOT volume-weighted budgets.\n",
  "6. Correlated posterior parameters are not independently identified.\n",
  "7. Compare Rhat and ESS before interpreting differences in medians.\n",
  sep = ""
)
