# 21_excess_air_ablation_paired_DO_N2.R
#
# Run IDENTICAL paired DO/excess-N2 models with vs without the
# excess-air observation in the recharge-mixing likelihood.
#
# WITH: source fractions inferred from delta18O, Cl, Na, excess air
# WITHOUT: source fractions inferred from delta18O, Cl, Na.
# Excess N2 still uses its ORIGINAL gas correction; xs_Air is removed
# only as an EXPLANATORY MIXING TRACER. The same DO/N2/MRT observations
# and priors are retained in both models.
#
# Each held-out group loses BOTH DO and N2 before Stan sampling.
# Source tracers and MRT remain; the five group folds are identical.
# Results must pass chain Rhat, ESS and divergence checks.
#
# Full comparison: 6 full fits + 30 grouped CV fits. It can be
# computationally intensive. Run_only_fold=1L for smoke testing.
# All data and outputs are private, ignored by Git.
library(tidyverse)
if (!file.exists("scripts/20_joint_paired_DO_excessN2.R")) {
  if (file.exists(
    "Oxygen-Reduction-Joint-Model/scripts/20_joint_paired_DO_excessN2.R"
  )) setwd("Oxygen-Reduction-Joint-Model")
  else stop("Open the Oxygen-Reduction-Joint-Model repository.")
}

comparison <- list(
  input = "data/Uwe_HP_data_prepared.csv",
  exclude_16758 = FALSE,  # rerun TRUE: influential high-N2 well
  low_DO_threshold = 2,  # sensitivity: 1, 2, 3 mg/L
  run_DOloss = FALSE,    # DO_only, N2_age and N2_lowDO
  run_CV = TRUE,
  run_only_fold = NULL,  # 1L smoke test
  chains = 4L, parallel_chains = 4L,
  warmup = 1400L, sampling = 1400L,
  adapt_delta = .98, max_treedepth = 13L,
  seed = 20261031L,
  endmember_prior_scale = 1,
  log_mrt_SD = .25,
  lsr_do_prior_a = 1, lsr_do_prior_b = 1,
  save_fits = TRUE
)
if (!file.exists(comparison$input))
  stop("Missing ORIGINAL PRIVATE data/Uwe_HP_data_prepared.csv")

run_comparison <- function() {
  previous <- getOption("oxygen.paired.cfg.override")
  on.exit(options(oxygen.paired.cfg.override=previous),add=TRUE)
  for (include_air in c(TRUE,FALSE)) {
    message("\nPAIRWISE ABLATION: excess air in mixing = ",include_air)
    options(oxygen.paired.cfg.override=c(
      comparison,list(include_excess_air=include_air)
    ))
    sys.source(
      "scripts/20_joint_paired_DO_excessN2.R",
      envir=new.env(parent=globalenv())
    )
  }
}
run_comparison()

# Exact same suffix conventions as script 20. Explicit file names
# avoid silently using stale model outputs from a previous analysis.
d <- read_csv(comparison$input,show_col_types=FALSE) %>%
  filter(is.finite(DO),DO>=0,
         is.finite(xs_N2),xs_N2>=0,
         is.finite(MRT),MRT>=0)
if (comparison$exclude_16758)
  d <- d %>% filter(Site_ID!="16758")
suffix <- paste0(
  "_N",nrow(d),
  if (comparison$exclude_16758) "_without16758" else "",
  "_threshold",comparison$low_DO_threshold
)
filename <- function(prefix,include_air) paste0(
  "outputs/",prefix,suffix,
  if (include_air) "_withair" else "_noair",".csv"
)
load_both <- function(prefix) map_dfr(
  c(TRUE,FALSE),
  function(air) read_csv(
    filename(prefix,air),show_col_types=FALSE
  ) %>% mutate(include_excess_air=air,.before=1)
)
out <- function(name) paste0(
  "outputs/joint_DO_N2_21_",name,suffix,".csv"
)

post <- load_both("joint_DO_N2_20_full_posteriors")
diag <- load_both("joint_DO_N2_20_convergence")
wells <- load_both("joint_DO_N2_20_per_well")
budget <- post %>%
  filter(variable %in% c(
    "k","LSR_DO_deficit","river_DO",
    "mean_well_LSR_pathway_loss",
    "mean_well_age_associated_loss",
    "sampled_well_fraction_LSR_of_loss",
    "beta_n2_link","sigma_N2"
  )) %>%
  select(include_excess_air,model,variable,
         median,lower95,upper95) %>%
  arrange(model,variable,desc(include_excess_air))
write_csv(budget,out("air_ablation_budget"))
write_csv(diag,out("air_ablation_convergence"))

# cmdstanr summary q5/q95 covers 90%, NOT 95%.
f <- wells %>%
  filter(str_detect(variable,"^f_LSR\\[")) %>%
  select(include_excess_air,model,sample_index,Site_ID,
         Depth,median,q5,q95) %>%
  mutate(width90=q95-q5) %>%
  pivot_wider(
    names_from=include_excess_air,
    values_from=c(median,width90),
    names_glue="{.value}_{include_excess_air}"
  ) %>%
  mutate(
    LSR_change_noair_minus_air=median_FALSE-median_TRUE
  )
write_csv(f,out("air_ablation_fractions"))
cat("\nLSR FRACTION CHANGES, FULL JOINT POSTERIOR:\n")
print(f %>% group_by(model) %>% summarise(
  n=n(),
  median_absolute_shift=median(abs(LSR_change_noair_minus_air)),
  max_absolute_shift=max(abs(LSR_change_noair_minus_air)),
  median_width90_air=median(width90_TRUE),
  median_width90_noair=median(width90_FALSE),
  .groups="drop"
),n=Inf)

conv <- diag %>% group_by(model,include_excess_air) %>%
  summarise(
    max_rhat=max(rhat,na.rm=TRUE),
    min_bulk_ESS=min(ess_bulk,na.rm=TRUE),
    .groups="drop"
  )
print(conv,n=Inf)
if (any(conv$max_rhat>1.05 |
        conv$min_bulk_ESS<400,na.rm=TRUE))
  warning("At least one full posterior fails convergence checks.")

if (comparison$run_CV) {
  scores <- load_both("joint_DO_N2_20_CV_summary")
  pred <- load_both("joint_DO_N2_20_CV_predictions")
  write_csv(scores,out("air_ablation_CV_scores"))
  paired <- pred %>%
    select(include_excess_air,model,fold,sample_index,
           Site_ID,Depth,observed_DO,predicted_DO,
           observed_N2,predicted_N2) %>%
    pivot_wider(
      names_from=include_excess_air,
      values_from=c(predicted_DO,predicted_N2),
      names_glue="{.value}_{include_excess_air}"
    ) %>%
    mutate(
      DO_delta_abs_error=
        abs(predicted_DO_FALSE-observed_DO)-
        abs(predicted_DO_TRUE-observed_DO),
      N2_delta_abs_error=
        abs(predicted_N2_FALSE-observed_N2)-
        abs(predicted_N2_TRUE-observed_N2)
    )
  write_csv(paired,out("air_ablation_paired_errors"))
  cat("\nMATCHED GROUPED CV (air versus no air):\n")
  print(
    scores %>% arrange(model,desc(include_excess_air)),n=Inf
  )
}
p <- budget %>% filter(
  variable %in% c(
    "k","LSR_DO_deficit",
    "mean_well_LSR_pathway_loss",
    "mean_well_age_associated_loss"
  )
) %>%
  mutate(label=paste(
    model,ifelse(include_excess_air,"With air","No air")
  )) %>%
  ggplot(aes(median,label)) +
  geom_segment(aes(
    x=lower95,xend=upper95,
    y=label,yend=label
  )) +
  geom_point(size=2.1) +
  facet_wrap(~variable,scales="free_x") +
  labs(
    title="Joint DO/excess N2: with and without excess-air mixing",
    subtitle="Posterior median and central 95% credible interval",
    x=NULL,y=NULL
  ) +
  theme_bw(base_size=11)
print(p)
ggsave(
  paste0("plots/joint_DO_N2_21_air_ablation",suffix,".png"),
  p,width=11,height=7,dpi=300
)
cat("\nFULL BAYESIAN OXYGEN POSTERIOR COMPARISON:\n")
print(budget,n=Inf)
cat(
  "\nCAUTIONS:\n",
  "* Air removed from MIXING likelihood only, not N2 gas correction.\n",
  "* DO and N2 jointly held out at the same wells.\n",
  "* Source priors are study-derived; rainfall isotope varies widely.\n",
  "* Age-associated k does not prove saturated aquifer oxygen loss.\n",
  "* LSR DO deficit depends on assumed 9.5 mg/L reference.\n",
  "* Site means are not flow-weighted oxygen budgets.\n",
  "* Compare Rhat, ESS, divergences before interpreting posteriors.\n",
  sep=""
)
