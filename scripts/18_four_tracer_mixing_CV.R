# 18_four_tracer_mixing_CV.R
#
# Matched Bayesian comparison of 2-SOURCE, 3-TRACER versus 4-TRACER
# source mixing on original conservative-tracer observations.
#
# M3: delta18O, chloride, excess air (NO sodium observed by model).
# M4: delta18O, chloride, sodium, excess air.
# - Both use identical endmember priors for shared variables.
# - Both permit uncertain LSR delta18O signature, fLSR*tau_LSR_O18.
# - Sodium and chloride have a JOINT bivariate-normal observation
#   model in M4; rho_Na_Cl models their conditional error correlation.
# - Sodium is NEVER used as Na/Cl ratio in addition to Cl, which
#   would double-use the chloride observation.
# - NO groundwater DO, MRT, excess N2, or categorical source labels
#   are input to the mixing models.
#
# FOUR deterministic rotations of one-tracer-at-a-time masking:
# each measured tracer (O18, Cl, Na, xs_Air) is withheld ONCE per well.
# At each rotation, an entire named well-group receives the same
# withheld tracer. The other available tracers remain observed.
#
# Fair M3 versus M4 scores compare ONLY the same three shared outcomes
# O18, Cl, xs_Air. M4's withheld Na forecasts are assessed separately;
# M3 is not claimed to predict Na without any sodium-data likelihood.
#
# Endpoint priors are informed by this SAME groundwater study and are
# therefore empirical priors. This is CONDITIONAL internal predictive
# validation, not external source-endmember validation.
#
# Outputs are LOCAL and gitignored; never commit private site records.
library(tidyverse)
setwd(here::here("Oxygen-Reduction-Joint-Model"))

cfg <- list(
  input = "data/groundwater.csv",
  model = "stan/18_four_tracer_joint_mixing.stan",
  endmember_prior_scale = 1,  # also rerun with 2 as sensitivity
  seed = 20261031L,
  chains = 4L,
  parallel_chains = 4L,
  warmup = 1100L,
  sampling = 1100L,
  adapt_delta = .97,
  max_treedepth = 13L,
  fit_full_data = TRUE,
  maximum_Rhat = 1.05,
  minimum_ESS_bulk = 400L,
  run_rotations = 1:4   # reduce to a subset for a quick smoke test
)

if (!requireNamespace("cmdstanr", quietly = TRUE) ||
    !requireNamespace("posterior", quietly = TRUE)) {
  stop("Requires cmdstanr + installed CmdStan + posterior. See README.")
}
if (!file.exists(cfg$input) || !file.exists(cfg$model))
  stop("Missing private groundwater CSV or Stan model.")
tryCatch(cmdstanr::cmdstan_path(), error = function(e) {
  stop("CmdStan is missing: ", conditionMessage(e))
})
dir.create("outputs", showWarnings = FALSE, recursive = TRUE)
dir.create("plots", showWarnings = FALSE, recursive = TRUE)

d <- read_csv(
  cfg$input, show_col_types = FALSE,
  col_types = cols(Site_ID = col_character(), .default = col_guess())
) %>%
  mutate(sample_index = row_number(), .before = 1)
required <- c("Site_ID", "Depth", "O18", "Cl", "Na",
              "xs_Air", "xs_Air_uncertainty")
stopifnot(all(required %in% names(d)), nrow(d) >= 35)
names_tracers <- c("O18", "Cl", "Na", "xs_Air")

present <- sapply(names_tracers, function(k) {
  is.finite(d[[k]])
})
stopifnot(is.matrix(present), all(rowSums(present) >= 2L),
          all(d$Cl[is.finite(d$Cl)] > 0),
          all(d$Na[is.finite(d$Na)] > 0))

numeric_or_zero <- function(v) ifelse(is.finite(v), v, 0)
positive_or_default <- function(v, default) {
  ifelse(is.finite(v) & v > 0, v, default)
}

stan_base <- list(
  N = nrow(d),
  train_O18 = as.integer(present[, "O18"]),
  train_Cl = as.integer(present[, "Cl"]),
  train_Na = as.integer(present[, "Na"]),
  train_air = as.integer(present[, "xs_Air"]),
  O18 = numeric_or_zero(d$O18),
  Cl = numeric_or_zero(d$Cl),
  Na = numeric_or_zero(d$Na),
  air = numeric_or_zero(d$xs_Air),
  Cl_sd = pmax(1, .03 * numeric_or_zero(d$Cl)),
  Na_sd = pmax(.5, .03 * numeric_or_zero(d$Na)),
  air_sd = positive_or_default(d$xs_Air_uncertainty, .15),
  use_na = 0L,
  endpoint_prior_scale = cfg$endmember_prior_scale
)

# One held-out tracer per group per rotation. Rotating across 1:4
# with a fixed group offset means EVERY present sample/tracer
# observation is withheld exactly once over the four rotations.
d <- d %>%
  mutate(
    validation_group = case_when(
      str_detect(Site_ID, "Eastbourne") ~ "Eastbourne",
      str_detect(Site_ID, "Gilbertson") ~ "Gilbertson",
      str_detect(Site_ID, "Flaxmere") ~ "Flaxmere",
      TRUE ~ Site_ID
    )
  )
set.seed(cfg$seed)
group_order <- sample(unique(d$validation_group))
group_offset <- setNames(
  (seq_along(group_order)-1L) %% 4L, group_order
)
row_offset <- unname(group_offset[d$validation_group])
masks <- map_dfr(1:4, function(rotation) {
  idx <- ((row_offset + rotation - 1L) %% 4L) + 1L
  tibble(
    rotation = rotation,
    sample_index = d$sample_index,
    Site_ID = d$Site_ID,
    Depth = d$Depth,
    withheld = names_tracers[idx],
    observed = present[cbind(
      seq_len(nrow(d)), idx
    )]
  )
})
stopifnot(
  nrow(masks %>% filter(observed)) == sum(present),
  all(masks %>% filter(observed) %>%
        count(sample_index, withheld) %>% pull(n) == 1L)
)
write_csv(masks, "outputs/four_tracer_CV_masks.csv")

stan_model <- cmdstanr::cmdstan_model(
  cfg$model, quiet = TRUE
)

models <- tribble(
  ~model,                    ~use_na,
  "3 tracers: no Na",         0L,
  "4 tracers: Na and Cl",     1L
)
prediction_name <- c(
  O18 = "predicted_O18",
  Cl = "predicted_Cl",
  Na = "predicted_Na",
  xs_Air = "predicted_air"
)
replicate_name <- c(
  O18 = "replicated_O18",
  Cl = "replicated_Cl",
  Na = "replicated_Na",
  xs_Air = "replicated_air"
)

results <- list()
diagnostics <- list()
draw_count <- 0L

for (rotation in cfg$run_rotations) {
  held_in_rotation <- masks %>%
    filter(.data$rotation == rotation, observed)

  for (j in seq_len(nrow(models))) {
    mod <- models[j, ]
    # M3 has no Na likelihood: when Na is withheld it is a
    # redundant full-data fit, so skip this rotation/model pair.
    held <- held_in_rotation %>%
      filter(mod$use_na == 1L | withheld != "Na")
    if (!nrow(held)) next

    X <- stan_base
    X$use_na <- mod$use_na

    masks_o <- held$sample_index[held$withheld=="O18"]
    masks_c <- held$sample_index[held$withheld=="Cl"]
    masks_n <- held$sample_index[held$withheld=="Na"]
    masks_a <- held$sample_index[held$withheld=="xs_Air"]
    X$train_O18[masks_o] <- 0L
    X$train_Cl[masks_c] <- 0L
    X$train_Na[masks_n] <- 0L
    X$train_air[masks_a] <- 0L
    # Do not even pass the withheld values to Stan.
    X$O18[masks_o] <- 0
    X$Cl[masks_c] <- 0
    X$Na[masks_n] <- 0
    X$air[masks_a] <- 0

    message("CV rotation ",rotation," | ",mod$model,
            " | n withheld = ",nrow(held))

    fit <- stan_model$sample(
      data = X,
      seed = cfg$seed + 1000L*rotation + 10000L*mod$use_na,
      chains = cfg$chains,
      parallel_chains = cfg$parallel_chains,
      iter_warmup = cfg$warmup,
      iter_sampling = cfg$sampling,
      adapt_delta = cfg$adapt_delta,
      max_treedepth = cfg$max_treedepth,
      refresh = 250
    )

    vars_pred <- unname(paste0(
      prediction_name[held$withheld], "[", held$sample_index, "]"
    ))
    vars_rep <- unname(paste0(
      replicate_name[held$withheld], "[", held$sample_index, "]"
    ))
    draws <- posterior::as_draws_matrix(
      fit$draws(variables = unique(c(vars_pred, vars_rep)))
    )

    # Predictions refer to the held-out ion conditional on an
    # observed correlated ion when available, not to an unconditionally
    # independent sodium/chloride observation.
    mean_pred <- vapply(vars_pred, function(v) {
      mean(draws[, v])
    }, numeric(1))
    lower95 <- vapply(vars_rep, function(v) {
      unname(quantile(draws[, v], .025))
    }, numeric(1))
    upper95 <- vapply(vars_rep, function(v) {
      unname(quantile(draws[, v], .975))
    }, numeric(1))

    truth <- mapply(function(i, t) {
      d[[t]][i]
    }, held$sample_index, held$withheld)

    draw_count <- draw_count+1L
    results[[draw_count]] <- held %>%
      mutate(
        model = mod$model,
        use_na = mod$use_na,
        observed_value = as.numeric(truth),
        predicted_mean = as.numeric(mean_pred),
        predictive_lower95 = as.numeric(lower95),
        predictive_upper95 = as.numeric(upper95),
        residual = observed_value - predicted_mean,
        inside95 = observed_value >= predictive_lower95 &
                   observed_value <= predictive_upper95
      )

    param_names <- c(
      "O18_R","O18_LSR","Cl_R","Cl_LSR",
      "air_R","air_LSR","sigma_O18",
      "tau_LSR_O18","sigma_Cl","sigma_air",
      if (mod$use_na == 1L)
        c("Na_R","Na_LSR","sigma_Na","rho_Na_Cl"),
      paste0("f_LSR[",held$sample_index,"]")
    )
    diagnostic <- fit$summary(
      variables = param_names
    )
    diagnostics[[draw_count]] <- diagnostic %>%
      mutate(
        model = mod$model,
        rotation = rotation,
        .before = 1
      )
    hmc <- fit$diagnostic_summary()
    if (any(diagnostic$rhat > cfg$maximum_Rhat,na.rm=TRUE) ||
        any(diagnostic$ess_bulk < cfg$minimum_ESS_bulk,
            na.rm=TRUE) ||
        any(hmc$num_divergent > 0)) {
      warning(
        "MCMC diagnostics need review: ",
        mod$model," rotation ",rotation
      )
    }

    write_csv(
      bind_rows(results),
      "outputs/four_tracer_CV_predictions.csv"
    )
    write_csv(
      bind_rows(diagnostics),
      "outputs/four_tracer_CV_convergence.csv"
    )
  }
}
cv <- bind_rows(results)
if (!nrow(cv)) stop("No CV runs completed.")

by_tracer <- cv %>%
  group_by(model,withheld) %>%
  summarise(
    n=n(),
    RMSE=sqrt(mean(residual^2)),
    MAE=mean(abs(residual)),
    bias=mean(predicted_mean-observed_value),
    coverage95=mean(inside95),
    .groups="drop"
  ) %>%
  mutate(
    observed_sample_sd=map_dbl(withheld,function(name) {
      sd(d[[name]],na.rm=TRUE)
    }),
    normalised_RMSE=RMSE/observed_sample_sd
  )
# Fair paired comparison: ONLY outcomes measured/predicted by both
# models (O18, Cl, excess air). Na is additional 4-model diagnostic.
common <- by_tracer %>%
  filter(withheld!="Na") %>%
  group_by(model) %>%
  summarise(
    n=sum(n),
    mean_normalised_RMSE=mean(normalised_RMSE),
    mean_coverage95=mean(coverage95),
    .groups="drop"
  )
# Compare paired absolute errors on EXACT same observations.
paired <- cv %>%
  filter(withheld!="Na") %>%
  select(sample_index,Site_ID,rotation,withheld,
         model,observed_value,predicted_mean) %>%
  pivot_wider(
    names_from=model, values_from=predicted_mean
  ) %>%
  mutate(
    absolute_error_3=abs(
      observed_value-.data[["3 tracers: no Na"]]
    ),
    absolute_error_4=abs(
      observed_value-.data[["4 tracers: Na and Cl"]]
    ),
    difference_absolute_error=absolute_error_4-
      absolute_error_3
  )
write_csv(
  by_tracer,"outputs/four_tracer_CV_scores_by_tracer.csv"
)
write_csv(
  common,"outputs/four_tracer_CV_scores_shared_outcomes.csv"
)
write_csv(
  paired,"outputs/four_tracer_CV_paired_predictions.csv"
)

if (isTRUE(cfg$fit_full_data)) {
  f_summary <- list()
  end_summary <- list()
  for (j in seq_len(nrow(models))) {
    mod <- models[j,]
    X <- stan_base
    X$use_na <- mod$use_na
    message("Full-data posterior source fractions: ",mod$model)
    fit <- stan_model$sample(
      data=X,
      seed=cfg$seed + 60000L + mod$use_na,
      chains=cfg$chains,
      parallel_chains=cfg$parallel_chains,
      iter_warmup=cfg$warmup,
      iter_sampling=cfg$sampling,
      adapt_delta=cfg$adapt_delta,
      max_treedepth=cfg$max_treedepth,
      refresh=250
    )
    f_summary[[j]] <- fit$summary(variables="f_LSR") %>%
      mutate(
        sample_index=as.integer(
          str_match(variable,"\\[(\\d+)\\]")[,2]
        ),
        model=mod$model, .before=1
      ) %>%
      inner_join(
        d %>% select(sample_index,Site_ID,Depth),
        by="sample_index"
      )
    end_summary[[j]] <- fit$summary(variables=c(
      "O18_R","O18_LSR","Cl_R","Cl_LSR",
      "air_R","air_LSR","tau_LSR_O18",
      if (mod$use_na==1L)
        c("Na_R","Na_LSR","rho_Na_Cl","sigma_Na")
    )) %>%
      mutate(model=mod$model,.before=1)
    if (any(fit$diagnostic_summary()$num_divergent>0)) {
      warning("Divergences in full-data ",mod$model)
    }
  }
  fs <- bind_rows(f_summary)
  write_csv(
    fs,"outputs/four_tracer_full_data_LSR_fractions.csv"
  )
  write_csv(
    bind_rows(end_summary),
    "outputs/four_tracer_full_data_endmembers.csv"
  )
  shift <- fs %>%
    select(sample_index,Site_ID,Depth,model,
           median,q5,q95) %>%
    pivot_wider(
      names_from=model,
      values_from=c(median,q5,q95)
    ) %>%
    mutate(
      median_LSR_difference=(
        .data[["median_4 tracers: Na and Cl"]] -
        .data[["median_3 tracers: no Na"]]
      )
    )
  write_csv(shift,"outputs/four_tracer_LSR_fraction_change.csv")
  cat("\nCHANGE IN INFERRED LSR FRACTIONS:\n")
  print(
    shift %>%
      summarise(
        n=n(),
        median_abs_shift=median(abs(median_LSR_difference)),
        maximum_abs_shift=max(abs(median_LSR_difference))
      )
  )
}

theme_set(theme_bw(base_size=12))
p <- cv %>%
  filter(withheld!="Na") %>%
  ggplot(aes(observed_value,predicted_mean)) +
  geom_abline(slope=1,intercept=0,linetype="dashed",
              colour="grey55") +
  geom_point(alpha=.75,size=1.6) +
  facet_grid(model~withheld,scales="free") +
  labs(
    title="Does sodium improve prediction of other mixing tracers?",
    subtitle="Matched held-out observations; joint Na-Cl errors",
    x="Observed tracer",
    y="Posterior predictive mean"
  )
print(p)
ggsave(
  "plots/four_tracer_CV_predictions.png",
  p,width=12,height=8,dpi=300
)

p2 <- by_tracer %>%
  filter(withheld!="Na") %>%
  ggplot(aes(withheld,normalised_RMSE,
             group=model,colour=model)) +
  geom_line(linewidth=.9) +
  geom_point(size=2.6) +
  labs(
    title="Predictive information added by sodium",
    subtitle="Common outcomes only; lower normalised RMSE is better",
    x="Held-out tracer",y="RMSE / observed tracer SD",
    colour="Bayesian model"
  )
print(p2)
ggsave(
  "plots/four_tracer_CV_normalised_error.png",
  p2,width=9,height=6,dpi=300
)
cat("\nHELD-OUT TRACER SCORES:\n")
print(by_tracer,n=Inf)
cat("\nSHARED-OUTCOME MODEL COMPARISON:\n")
print(common,n=Inf)
cat("\nPAIRWISE CHANGE IN ABSOLUTE ERROR:\n")
print(summary(paired$difference_absolute_error))
cat(
  "\nLIMITATIONS:\n",
  "* No DO, MRT, N2 or predefined source class was fitted.\n",
  "* Na and Cl measurement errors are approximated at 3%, subject\n",
  "  to minimum assumed analytical SD; verify laboratory precision.\n",
  "* Their residual correlation is statistical, NOT proof of an\n",
  "  ion-exchange mechanism or independent analytical calibration.\n",
  "* Rainfall delta18O variability is marginalised, not inferred\n",
  "  from source-area/time-series samples.\n",
  "* Endmember priors use study-derived well endpoints (empirical Bayes).\n",
  "  Repeat with cfg$endmember_prior_scale=2, review divergences,\n",
  "  Rhat and bulk ESS before interpreting predictive improvement.\n",
  "* Source fractions are latent estimates, not observed ground truth.\n",
  "* This is not spatially independent validation; endmember priors\n",
  "  and source chemistry may differ between recharge environments.\n",
  sep=""
)
