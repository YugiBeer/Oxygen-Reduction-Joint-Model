# 22_joint_DO_denitrification_progress.R
#
# Fully joint Bayesian DO / denitrification reaction-progress analysis.
#
# REMOVE WELL 16758 by study decision. Use the 40 complete records:
# DO, MRT, and observed reaction_progress in 0..1. The missing
# Pelomino Rd reaction progress and Frimley MRT are not imputed.
#
# Reaction progress is a DERIVED response (not independent tracer):
# xs_N2 / (xs_N2 + NO3N + NO2).
# The two exact zeros are modelled via a zero-hurdle. Positive
# values are modelled with a robust logit-scale likelihood.
# NO numerator/denominator variables, initial_NO3N, or derived
# reaction_progress are input as predictors. xs_N2_sig enters only
# as APPROXIMATE response observation uncertainty.
#
# Three mixing tracers by default: delta18O, Na, Cl.
# Set include_excess_air=TRUE for optional four-tracer sensitivity.
# Sodium/chloride residual correlation modelled jointly.
#
# Four models: DO-only; uncoupled progress age/source; progress
# linked to cumulative DO depletion; progress linked to time
# below a fixed DO threshold. All include DO, MRT, mixing data.
#
# Matched grouped five-fold CV: DO and progress BOTH masked/zeroed
# at the same held-out wells. Response prediction cannot use
# measured DO or progress at held-out sites.
#
# Run in Positron with CmdStanR + CmdStan and original PRIVATE
# data/Uwe_HP_data_prepared.csv. Outputs local/gitignored.
library(tidyverse)
if (file.exists("stan/22_joint_DO_denitrification_progress.stan")) {
  # at repo root
} else if (file.exists(
  "Oxygen-Reduction-Joint-Model/stan/22_joint_DO_denitrification_progress.stan"
)) {
  setwd("Oxygen-Reduction-Joint-Model")
} else stop("Set working directory to Oxygen-Reduction-Joint-Model.")

cfg <- list(
  input = "data/Uwe_HP_data_prepared.csv",
  stan = "stan/22_joint_DO_denitrification_progress.stan",
  exclude_site = "16758",    # obligatory exclusion per study decision
  include_excess_air = FALSE, # prior ablation indicates little benefit
  log_mrt_SD = .25,
  endmember_prior_scale = 1,
  lsr_do_prior_a = 1, lsr_do_prior_b = 1,
  low_DO_threshold = 2,       # test sensitivity at 1 and 3 mg/L
  nfolds = 5L,
  seed = 20261031L,
  chains = 4L, parallel_chains = 4L,
  warmup = 1400L, sampling = 1400L,
  adapt_delta = .98, max_treedepth = 13L,
  max_rhat = 1.05, min_bulk_ESS = 400,
  run_CV = TRUE,              # four full + 20 CV fits
  run_only_fold = NULL,       # use 1L for a CV smoke test
  run_cumulative_loss = TRUE,# FALSE omits progress_DOloss
  save_full_fit = TRUE
)
if (!requireNamespace("cmdstanr",quietly=TRUE) ||
    !requireNamespace("posterior",quietly=TRUE))
  stop("Install cmdstanr + CmdStan and posterior.")
if (!file.exists(cfg$input))
  stop("Copy the ORIGINAL private prepared CSV to ",cfg$input,
       " locally; the reduced groundwater.csv lacks progress/Na.")
tryCatch(cmdstanr::cmdstan_path(),error=function(e)
  stop("CmdStan is missing: ",conditionMessage(e))
)
dir.create("outputs",showWarnings=FALSE,recursive=TRUE)
dir.create("plots",showWarnings=FALSE,recursive=TRUE)

raw <- read_csv(
  cfg$input,show_col_types=FALSE,
  col_types=cols(Site_ID=col_character(),.default=col_guess())
)
needed <- c(
  "Site_ID","Depth","DO","MRT","O18","Na","Cl","xs_Air",
  "xs_Air_uncertainty","xs_N2","xs_N2_sig",
  "NO3N","NO2","reaction_progress"
)
if (!all(needed %in% names(raw)))
  stop("Missing fields: ",paste(setdiff(needed,names(raw)),
                                 collapse=", "))
if (!any(raw$Site_ID==cfg$exclude_site))
  warning("Excluded site not found: ",cfg$exclude_site)
raw <- raw %>% filter(Site_ID!=cfg$exclude_site)
d <- raw %>%
  filter(
    is.finite(DO),DO>=0,
    is.finite(MRT),MRT>=0,
    is.finite(reaction_progress),
    reaction_progress>=0,reaction_progress<1,
    is.finite(xs_N2),xs_N2>=0,
    is.finite(NO3N),is.finite(NO2),
    is.finite(xs_N2_sig),xs_N2_sig>=0
  ) %>%
  mutate(
    sample_index=row_number(),
    group=case_when(
      str_detect(Site_ID,"Eastbourne") ~ "Eastbourne",
      str_detect(Site_ID,"Gilbertson") ~ "Gilbertson",
      str_detect(Site_ID,"Flaxmere") ~ "Flaxmere",
      TRUE ~ Site_ID
    )
  )
stopifnot(nrow(d)>=30, !any(d$Site_ID==cfg$exclude_site))
# Guard against changes to the reaction-progress construction
# or accidental mixing of N2 gas-concentration and N equivalents.
recalc <- with(d,xs_N2/(xs_N2+NO3N+NO2))
if (any(!is.finite(recalc) |
        abs(d$reaction_progress-recalc)>1e-9))
  stop("reaction_progress no longer matches N2/(N2+NO3N+NO2).")
message("Included ",nrow(d)," paired DO/reaction-progress samples; ",
        "excluded site ",cfg$exclude_site,".")
message("Two-valued response: exact zeros = ",
        sum(d$reaction_progress==0),", positive = ",
        sum(d$reaction_progress>0),".")
if (nrow(d)!=40L)
  warning("Expected 40 complete rows from current study CSV, found ",
          nrow(d),". Inspect missing data/version.")

# Approximate excess-N2 measurement error in progress logit.
# logit(P) = log(N2)-log(NO3N+NO2) (for positive N2).
# First-order SD from N2 measurement alone: sigma_N2/N2.
# NO3/NO2 analytical uncertainties unavailable; the fitted
# sigma_progress_logit captures extra scatter.
d <- d %>% mutate(
  analytical_progress_logit_sd=if_else(
    xs_N2>0,pmin(2,pmax(.1,xs_N2_sig/pmax(xs_N2,.2))),
    .5
  )
)
as_obs <- function(x) ifelse(is.finite(x),as.numeric(x),0)
positive <- function(x,value) ifelse(
  is.finite(x)&x>0,as.numeric(x),value
)
stan_base <- list(
  N=nrow(d),
  has_O18=as.integer(is.finite(d$O18)),
  has_Cl=as.integer(is.finite(d$Cl)),
  has_Na=as.integer(is.finite(d$Na)),
  has_air=if (isTRUE(cfg$include_excess_air))
    as.integer(is.finite(d$xs_Air)) else rep(0L,nrow(d)),
  has_DO=rep(1L,nrow(d)),
  has_progress=rep(1L,nrow(d)),
  O18=as_obs(d$O18),
  Cl=as_obs(d$Cl),
  Na=as_obs(d$Na),
  air=if (isTRUE(cfg$include_excess_air))
    as_obs(d$xs_Air) else rep(0,nrow(d)),
  DO_obs=as.numeric(d$DO),
  progress_obs=as.numeric(d$reaction_progress),
  Cl_sd=pmax(1,.03*as_obs(d$Cl)),
  Na_sd=pmax(.5,.03*as_obs(d$Na)),
  air_sd=positive(d$xs_Air_uncertainty,.15),
  progress_logit_sd=d$analytical_progress_logit_sd,
  mrt_reported=as.numeric(d$MRT),
  use_na=1L,
  use_progress=0L,
  progress_link=0L,
  endmember_prior_scale=cfg$endmember_prior_scale,
  log_mrt_SD=cfg$log_mrt_SD,
  lsr_do_prior_a=cfg$lsr_do_prior_a,
  lsr_do_prior_b=cfg$lsr_do_prior_b,
  low_DO_threshold=cfg$low_DO_threshold
)
models <- tribble(
  ~model,               ~use_progress,~progress_link,
  "DO_only",                0L,          0L,
  "progress_age",           1L,          0L,
  "progress_DOloss",        1L,          1L,
  "progress_lowDO",         1L,          2L
)
if (!cfg$run_cumulative_loss)
  models <- models %>% filter(model!="progress_DOloss")
suffix <- paste0(
  "_N",nrow(d),
  if(cfg$include_excess_air) "_withair" else "_noair",
  "_threshold",cfg$low_DO_threshold
)
file_csv <- function(name)
  file.path("outputs",paste0("joint_DO_progress_22_",name,suffix,".csv"))
fit_stan <- cmdstanr::cmdstan_model(cfg$stan,quiet=TRUE)
fit_spec <- function(s,fold=0L,x=stan_base) {
  x$use_progress <- as.integer(s$use_progress)
  x$progress_link <- as.integer(s$progress_link)
  # DO-only does not use any progress observations.
  if (s$use_progress==0L) {
    x$has_progress[] <- 0L
    x$progress_obs[] <- 0
  }
  fit_stan$sample(
    data=x,
    seed=cfg$seed + 10000L*fold +
      101L*match(s$model,models$model),
    chains=cfg$chains,
    parallel_chains=cfg$parallel_chains,
    iter_warmup=cfg$warmup,iter_sampling=cfg$sampling,
    adapt_delta=cfg$adapt_delta,
    max_treedepth=cfg$max_treedepth,
    refresh=250
  )
}
global_names <- function(s) c(
  "k","LSR_entry_DO","LSR_DO_deficit","river_DO",
  "mean_well_LSR_pathway_loss",
  "mean_well_age_associated_loss",
  "sampled_well_fraction_LSR_of_loss",
  "mean_well_LSR_fraction",
  "sigma_DO","tau_LSR_O18","rho_Na_Cl",
  if(s$use_progress==1L) c(
    "progress_zero_probability","progress_intercept",
    "progress_source_beta","progress_age_beta",
    "sigma_progress_logit",
    if(s$progress_link!=0L) "progress_link_beta"
  )
)
check_convergence <- function(fit,s,fold=0L,held=NULL) {
  vars <- c(
    global_names(s),
    if(is.null(held)) c("f_LSR","true_MRT")
    else c(paste0("f_LSR[",held,"]"),
           paste0("true_MRT[",held,"]"))
  )
  out <- fit$summary(variables=vars)
  hmc <- fit$diagnostic_summary()
  bad <- any(out$rhat>cfg$max_rhat,na.rm=TRUE) ||
    any(out$ess_bulk<cfg$min_bulk_ESS,na.rm=TRUE) ||
    any(hmc$num_divergent>0,na.rm=TRUE)
  if (bad) warning(
    "Convergence problems: ",s$model," fold ",fold,
    "; do NOT rely on posterior intervals or model rankings.")
  out %>% mutate(model=s$model,fold=fold,
                 any_HMC_divergences=any(
                   hmc$num_divergent>0,na.rm=TRUE),
                 .before=1)
}
posterior95 <- function(fit,vars) {
  x <- posterior::as_draws_matrix(
    fit$draws(variables=vars)
  )
  tibble(
    variable=colnames(x),
    median=apply(x,2,median),
    lower95=apply(x,2,quantile,.025),
    upper95=apply(x,2,quantile,.975)
  )
}
full <- list();diag <- list();well <- list()
for (j in seq_len(nrow(models))) {
  s <- models[j,]
  message("Full joint DO + progress posterior: ",s$model)
  fit <- fit_spec(s)
  full[[j]] <- posterior95(
    fit,global_names(s)
  ) %>% mutate(model=s$model,.before=1)
  diag[[j]] <- check_convergence(fit,s)
  well[[j]] <- fit$summary(variables=c(
    "f_LSR","mean_DO","true_MRT",
    "low_DO_duration_years","LSR_pathway_DO_loss",
    "age_associated_DO_loss","progress_logit_center"
  )) %>%
    mutate(
      sample_index=as.integer(str_match(
        variable,"\\[(\\d+)\\]")[,2]),
      model=s$model,.before=1
    ) %>%
    inner_join(
      d %>% select(
        sample_index,Site_ID,Depth,DO,reaction_progress,MRT
      ),by="sample_index"
    )
  if (cfg$save_full_fit) fit$save_object(
    file=file.path("outputs",
      paste0("joint_DO_progress_22_",s$model,suffix,".rds"))
  )
  write_csv(bind_rows(full),file_csv("full_posteriors"))
  write_csv(bind_rows(diag),file_csv("full_convergence"))
  write_csv(bind_rows(well),file_csv("per_well"))
}

# Group all correlated labels together in stable five-fold splits.
set.seed(cfg$seed)
gg <- split(seq_len(nrow(d)),d$group)
g_order <- sample(seq_along(gg))
g_order <- g_order[base::order(-lengths(gg)[g_order])]
folds <- vector("list",cfg$nfolds)
for (i in g_order) {
  j <- which.min(lengths(folds))
  folds[[j]] <- c(folds[[j]],gg[[i]])
}
stopifnot(
  all(lengths(folds)>0),
  identical(sort(unlist(folds)),seq_len(nrow(d)))
)
fold_map <- map_dfr(seq_along(folds),function(j)
  tibble(fold=j,sample_index=folds[[j]],
         Site_ID=d$Site_ID[folds[[j]]],Depth=d$Depth[folds[[j]]])
)
write_csv(fold_map,file_csv("fold_map"))
if (cfg$run_CV) {
  selected <- if (is.null(cfg$run_only_fold))
    seq_along(folds) else cfg$run_only_fold
  stopifnot(all(selected %in% seq_along(folds)))
  cvpred <- list();cvdiag <- list()
  for (fold in selected) {
    held <- folds[[fold]]
    for (j in seq_len(nrow(models))) {
      s <- models[j,]
      x <- stan_base
      x$has_DO[held] <- 0L
      x$DO_obs[held] <- 0
      x$has_progress[held] <- 0L
      x$progress_obs[held] <- 0
      # The original progress error was estimated from measured
      # xs_N2/xs_N2_sig. At held-out wells it MUST NOT be passed
      # to the posterior predictive generator: use a training-only
      # representative error instead to prevent response leakage.
      x$progress_logit_sd[held] <- median(
        stan_base$progress_logit_sd[-held]
      )
      message("Paired CV ",fold,"/",cfg$nfolds,": ",s$model)
      fit <- fit_spec(s,fold,x)
      var_DO <- paste0("mean_DO[",held,"]")
      rep_DO <- paste0("replicated_DO[",held,"]")
      rep_P <- paste0("replicated_progress[",held,"]")
      vars <- c(var_DO,rep_DO)
      if(s$use_progress==1L) vars <- c(vars,rep_P)
      draws <- posterior::as_draws_matrix(
        fit$draws(variables=vars)
      )
      series <- function(v,fun) vapply(v,function(z)
        fun(draws[,z]),numeric(1))
      q <- function(pr) function(v)
        unname(quantile(v,probs=pr))
      n <- length(held)
      pred_p <- if(s$use_progress==1L)
        series(rep_P,mean) else rep(NA_real_,n)
      lo_p <- if(s$use_progress==1L)
        series(rep_P,q(.025)) else rep(NA_real_,n)
      hi_p <- if(s$use_progress==1L)
        series(rep_P,q(.975)) else rep(NA_real_,n)
      pp <- tibble(
        model=s$model,fold=fold,sample_index=held,
        Site_ID=d$Site_ID[held],Depth=d$Depth[held],
        observed_DO=d$DO[held],
        predicted_DO=series(var_DO,mean),
        DO_lower95=series(rep_DO,q(.025)),
        DO_upper95=series(rep_DO,q(.975)),
        observed_progress=d$reaction_progress[held],
        predicted_progress=pred_p,
        progress_lower95=lo_p,
        progress_upper95=hi_p
      ) %>%
        mutate(
          DO_covered95=observed_DO>=DO_lower95 &
                       observed_DO<=DO_upper95,
          progress_covered95=observed_progress>=progress_lower95 &
                             observed_progress<=progress_upper95
        )
      cvpred[[length(cvpred)+1L]] <- pp
      cvdiag[[length(cvdiag)+1L]] <-
        check_convergence(fit,s,fold,held)
      write_csv(bind_rows(cvpred),file_csv("CV_predictions"))
      write_csv(bind_rows(cvdiag),file_csv("CV_convergence"))
    }
  }
  pp <- bind_rows(cvpred)
  scores <- pp %>% group_by(model) %>% summarise(
    n=n(),folds=n_distinct(fold),
    DO_RMSE=sqrt(mean((observed_DO-predicted_DO)^2)),
    DO_MAE=mean(abs(observed_DO-predicted_DO)),
    DO_bias=mean(predicted_DO-observed_DO),
    DO_coverage95=mean(DO_covered95),
    progress_RMSE=if(all(is.na(predicted_progress)))
      NA_real_ else sqrt(mean(
        (observed_progress-predicted_progress)^2)),
    progress_MAE=if(all(is.na(predicted_progress)))
      NA_real_ else mean(
        abs(observed_progress-predicted_progress)),
    progress_coverage95=if(all(is.na(progress_covered95)))
      NA_real_ else mean(progress_covered95),
    .groups="drop"
  )
  write_csv(scores,file_csv("CV_summary"))
  pair <- pp %>% select(
    sample_index,Site_ID,Depth,fold,model,
    observed_DO,observed_progress,
    predicted_DO,predicted_progress
  ) %>%
    pivot_wider(
      names_from=model,
      values_from=c(predicted_DO,predicted_progress)
    )
  write_csv(pair,file_csv("CV_paired_by_well"))
  cat("\nMATCHED PAIRED DO & PROGRESS GROUPED CV:\n")
  print(scores,n=Inf)

  p <- pp %>%
    filter(model!="DO_only") %>%
    select(model,observed_DO,predicted_DO,
           observed_progress,predicted_progress) %>%
    pivot_longer(
      cols=c(observed_DO,predicted_DO,
             observed_progress,predicted_progress),
      names_to=c(".value","outcome"),
      names_pattern="(observed|predicted)_(DO|progress)"
    ) %>%
    ggplot(aes(observed,predicted))+
    geom_abline(slope=1,intercept=0,
                linetype="dashed",colour="grey60")+
    geom_point(alpha=.8)+
    facet_grid(model~outcome,scales="free")+
    labs(
      title="Paired Bayesian prediction: DO and denitrification progress",
      subtitle="Both outcomes withheld at same wells; site 16758 excluded",
      x="Measured",y="Posterior predictive mean"
    )+theme_bw(base_size=11)
  print(p)
  ggsave(
    file.path("plots",paste0(
      "joint_DO_progress_22_CV",suffix,".png")),
    p,width=12,height=8,dpi=300
  )
}
global <- bind_rows(full)
p <- global %>% filter(variable %in% c(
  "k","LSR_DO_deficit","LSR_entry_DO",
  "mean_well_LSR_pathway_loss",
  "mean_well_age_associated_loss",
  "progress_link_beta"
)) %>%
  ggplot(aes(median,model))+
  geom_segment(aes(x=lower95,xend=upper95,
                   y=model,yend=model))+
  geom_point(size=2.4)+
  facet_wrap(~variable,scales="free_x")+
  labs(
    title="DO and denitrification progress: posterior comparison",
    subtitle="Medians and central 95% credible intervals",
    x=NULL,y=NULL
  )+theme_bw(base_size=11)
print(p)
ggsave(
  file.path("plots",paste0(
    "joint_DO_progress_22_posteriors",suffix,".png")),
  p,width=12,height=7,dpi=300
)
cat("\nPOSTERIOR OXYGEN PARAMETERS:\n")
print(global %>% filter(variable %in% c(
  "k","LSR_entry_DO","LSR_DO_deficit",
  "mean_well_LSR_pathway_loss",
  "mean_well_age_associated_loss",
  "progress_link_beta","sigma_progress_logit"
)),n=Inf)
cat(
  "\nLIMITATIONS:\n",
  "* Reaction progress is derived from excess N2 and remaining\n",
  "  nitrate/nitrite; the reconstruction assumes measured N2 is\n",
  "  denitrification-derived and initial oxidised N is preserved.\n",
  "* Derived initial_NO3N and reaction_progress are NOT predictors.\n",
  "* The positive response has logit-t residuals; exact zeros use a\n",
  "  hurdle. Analytical progress uncertainty is only approximated.\n",
  "* Low-DO history is a piston-flow conceptual approximation.\n",
  "* Source-endmember priors are partly derived from study wells.\n",
  "* Rhat, bulk ESS and divergences must pass before interpretation.\n",
  "* The LSR DO reference (9.5 mg/L) is assumed, not measured.\n",
  "* Arithmetic well-average DO losses are NOT flow-weighted.\n",
  "* No automatic link demonstrates denitrification consumes O2.\n",
  "* Local groundwater data and outputs remain private.\n",
  sep=""
)
