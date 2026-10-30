# 19_joint_three_vs_four_tracer_DO.R
#
# Matched, fully joint Bayesian dissolved-oxygen model comparison.
# 3-tracer control = Cl, delta18O, excess air
# 4-tracer model  = control + Na, using a JOINT Na/Cl likelihood
#
# All priors (except prior-only Na terms in the control), MRT uncertainty,
# sample selection, age-associated k process, LSR-entry and river-entry DO
# are identical. LSR delta18O composition varies by recharge via the
# marginalized extra fLSR*tau_LSR_O18 variance.
#
# Runs two FULL JOINT FITS, then matched grouped 5-fold held-out DO CV.
# Held-out DO inputs are zeroed and their likelihood masked, but their
# measured tracers and reported MRT remain available.
#
# Critical: Na and Cl are not counted independently, and no ratio is
# used as a pseudo-additional tracer. Excess N2 is reactive, omitted.
# Mean well DO losses are NOT flow-weighted and not proof of reaction
# location in the saturated aquifer.
#
# Input from ORIGINAL PRIVATE study: data/Uwe_HP_data_prepared.csv.
# All local outputs are gitignored. Never commit groundwater records.
library(tidyverse)

if (file.exists("stan/19_joint_three_vs_four_tracer_DO.stan")) {
  # project root
} else if (file.exists(
  "Oxygen-Reduction-Joint-Model/stan/19_joint_three_vs_four_tracer_DO.stan"
)) {
  setwd("Oxygen-Reduction-Joint-Model")
} else stop("Open Oxygen-Reduction-Joint-Model as working directory.")

cfg <- list(
  input = "data/Uwe_HP_data_prepared.csv",
  model = "stan/19_joint_three_vs_four_tracer_DO.stan",
  log_mrt_SD = .25,
  endmember_prior_scale = 1,
  lsr_do_prior_a = 1, lsr_do_prior_b = 1, # beta(1,1) = uniform
  nfolds = 5L,
  seed = 20261031L,
  chains = 4L, parallel_chains = 4L,
  warmup = 1250L, sampling = 1250L,
  adapt_delta = .98, max_treedepth = 13L,
  max_rhat = 1.05, min_bulk_ESS = 400L,
  run_CV = TRUE,         # 10 expensive additional CmdStan fits
  run_only_fold = NULL,  # NULL for all; 1L for a smoke test
  run_prior_sensitivity = FALSE, # optional doubled endmember prior widths
  exclude_river_endpoint_from_DO = FALSE, # sensitivity: site 2801 DO=10.2
  save_fits = TRUE
)

if (!file.exists(cfg$input)) {
  stop(
    "Missing original private data/Uwe_HP_data_prepared.csv. ",
    "Copy it from Lincoln-Environmental/Oxygen-Reduction/data/ ",
    "without committing it. The old groundwater.csv lacks sodium."
  )
}
if (!requireNamespace("cmdstanr",quietly=TRUE) ||
    !requireNamespace("posterior",quietly=TRUE)) {
  stop("Requires cmdstanr, CmdStan and posterior; see README.")
}
tryCatch(cmdstanr::cmdstan_path(),error=function(e)
  stop("Install/configure CmdStan: ",conditionMessage(e))
)
dir.create("outputs",showWarnings=FALSE,recursive=TRUE)
dir.create("plots",showWarnings=FALSE,recursive=TRUE)

raw <- read_csv(
  cfg$input,show_col_types=FALSE,
  col_types=cols(Site_ID=col_character(),.default=col_guess())
)
required <- c("Site_ID","Depth","O18","Cl","Na","xs_Air",
              "xs_Air_uncertainty","DO","MRT")
if (!all(required %in% names(raw)))
  stop("Missing data fields: ",
       paste(setdiff(required,names(raw)),collapse=", "))
d <- raw %>%
  filter(is.finite(DO),DO>=0,is.finite(MRT),MRT>=0) %>%
  mutate(
    sample_index=row_number(),
    cv_group=case_when(
      str_detect(Site_ID,"Eastbourne") ~ "Eastbourne",
      str_detect(Site_ID,"Gilbertson") ~ "Gilbertson",
      str_detect(Site_ID,"Flaxmere") ~ "Flaxmere",
      TRUE ~ Site_ID
    )
  )
stopifnot(nrow(d)>=35L)
message("Using ",nrow(d)," wells with observed DO and MRT.")
if (any(d$Site_ID=="2801")) {
  message(
    "River endpoint well 2801: DO=",
    d$DO[d$Site_ID=="2801"][1],
    " mg/L, MRT=",d$MRT[d$Site_ID=="2801"][1],
    " years. This is high relative to river_DO prior 7.5."
  )
}

obs <- function(x) ifelse(is.finite(x),as.numeric(x),0)
positive <- function(x,default) ifelse(
  is.finite(x)&x>0,as.numeric(x),default
)
base <- list(
  N=nrow(d),
  has_O18=as.integer(is.finite(d$O18)),
  has_Cl=as.integer(is.finite(d$Cl)),
  has_Na=as.integer(is.finite(d$Na)),
  has_air=as.integer(is.finite(d$xs_Air)),
  has_DO=rep(1L,nrow(d)),
  O18=obs(d$O18),
  Cl=obs(d$Cl),
  Na=obs(d$Na),
  air=obs(d$xs_Air),
  DO_obs=as.numeric(d$DO),
  Cl_sd=pmax(1,.03*obs(d$Cl)),
  Na_sd=pmax(.5,.03*obs(d$Na)),
  air_sd=positive(d$xs_Air_uncertainty,.15),
  mrt_reported=as.numeric(d$MRT),
  use_na=0L,
  log_mrt_SD=cfg$log_mrt_SD,
  endmember_prior_scale=cfg$endmember_prior_scale,
  lsr_do_prior_a=cfg$lsr_do_prior_a,
  lsr_do_prior_b=cfg$lsr_do_prior_b
)
if (cfg$exclude_river_endpoint_from_DO) {
  i <- which(d$Site_ID=="2801")
  base$has_DO[i] <- 0L
  base$DO_obs[i] <- 0
}

model <- cmdstanr::cmdstan_model(cfg$model,quiet=TRUE)
names_model <- c("Three tracers","Four tracers")
active_global <- function(use_na) c(
  "O18_R","O18_LSR","Cl_R","Cl_LSR","air_R","air_LSR",
  "sigma_O18","tau_LSR_O18","sigma_Cl","sigma_air",
  "river_DO","LSR_entry_DO","LSR_DO_deficit",
  "k","sigma_DO","mean_well_LSR_pathway_loss",
  "mean_well_age_associated_loss",
  "sampled_well_fraction_LSR_of_loss","mean_well_LSR_fraction",
  if (use_na==1L)
    c("Na_R","Na_LSR","sigma_Na","rho_Na_Cl")
)
active_sites <- function() c(
  "f_LSR","true_MRT","mean_DO",
  "age_associated_DO_loss","LSR_pathway_DO_loss"
)
posterior_tbl <- function(fit,variables) {
  values <- posterior::as_draws_matrix(fit$draws(variables=variables))
  tibble(
    variable=colnames(values),
    median=apply(values,2,median),
    lower95=apply(values,2,quantile,.025),
    upper95=apply(values,2,quantile,.975)
  )
}
fit_args <- function(use_na,fold_id=0L,data=base) {
  data$use_na <- as.integer(use_na)
  model$sample(
    data=data,
    seed=cfg$seed+1000L*fold_id+10000L*use_na,
    chains=cfg$chains,parallel_chains=cfg$parallel_chains,
    iter_warmup=cfg$warmup,iter_sampling=cfg$sampling,
    adapt_delta=cfg$adapt_delta,max_treedepth=cfg$max_treedepth,
    refresh=250
  )
}
check_mcmc <- function(fit,vars,label) {
  s <- fit$summary(variables=vars)
  hmc <- fit$diagnostic_summary()
  problems <- (
    any(s$rhat>cfg$max_rhat,na.rm=TRUE) ||
    any(s$ess_bulk<cfg$min_bulk_ESS,na.rm=TRUE) ||
    any(hmc$num_divergent>0,na.rm=TRUE)
  )
  if (problems) warning("Review MCMC diagnostics: ",label)
  s %>% mutate(model_check=label,.before=1)
}

# FULL DATA FITS: save all private posterior draws, source fractions,
# and oxygen-loss decomposition for a paired parametric comparison.
globals <- list()
sites <- list()
full_diag <- list()
for (use_na in 0:1) {
  label <- names_model[use_na+1L]
  message("Full joint fit: ",label)
  fit <- fit_args(use_na)
  gl <- active_global(use_na)
  globals[[use_na+1L]] <- posterior_tbl(fit,gl) %>%
    mutate(model=label,.before=1)
  varsites <- active_sites()
  ss <- fit$summary(variables=varsites)
  sites[[use_na+1L]] <- ss %>%
    mutate(
      site_parameter=str_remove(variable,"\\[.*$"),
      sample_index=as.integer(str_match(
        variable,"\\[(\\d+)\\]")[,2]),
      model=label,
      .before=1
    ) %>%
    inner_join(
      d %>% select(sample_index,Site_ID,Depth,DO,MRT),
      by="sample_index"
    )
  full_diag[[use_na+1L]] <- check_mcmc(
    fit,
    c(gl,"f_LSR","true_MRT"),
    paste0("full_",label)
  )
  if (cfg$save_fits) fit$save_object(
    file=file.path("outputs",
      paste0("joint_DO_19_full_Na",use_na,".rds"))
  )
  if (use_na==1L && cfg$run_prior_sensitivity) {
    wider <- base
    wider$endmember_prior_scale <- 2
    message("Full four-tracer sensitivity: doubled endmember prior SD")
    fit2 <- fit_args(1L,fold_id=99L,data=wider)
    globals[[3]] <- posterior_tbl(fit2,active_global(1L)) %>%
      mutate(model="Four tracers, doubled prior SD",.before=1)
    full_diag[[3]] <- check_mcmc(
      fit2,c(active_global(1L),"f_LSR","true_MRT"),
      "full_wide_prior"
    )
    if (cfg$save_fits)
      fit2$save_object(file="outputs/joint_DO_19_four_wide_prior.rds")
  }
}
globs <- bind_rows(globals)
well <- bind_rows(sites)
write_csv(globs,"outputs/joint_DO_19_posterior_parameters.csv")
write_csv(well,"outputs/joint_DO_19_per_well_posterior.csv")
write_csv(bind_rows(full_diag),
          "outputs/joint_DO_19_full_diagnostics.csv")
mix_shift <- well %>%
  filter(site_parameter=="f_LSR") %>%
  select(sample_index,Site_ID,Depth,model,median,q5,q95) %>%
  pivot_wider(names_from=model,values_from=c(median,q5,q95)) %>%
  mutate(
    fraction_shift_Four_minus_Three=
      .data[["median_Four tracers"]]-
      .data[["median_Three tracers"]]
  )
write_csv(mix_shift,
          "outputs/joint_DO_19_source_fraction_change.csv")

# FIXED IDENTICAL CV GROUPS AND OBSERVATIONS IN BOTH MODELS.
set.seed(cfg$seed)
grouped <- split(seq_len(nrow(d)),d$cv_group)
group_order <- sample(seq_along(grouped))
group_order <- group_order[
  order(-lengths(grouped)[group_order])
]
folds <- vector("list",cfg$nfolds)
for (i in group_order) {
  j <- which.min(lengths(folds))
  folds[[j]] <- c(folds[[j]],grouped[[i]])
}
stopifnot(
  all(lengths(folds)>0L),
  identical(sort(unlist(folds)),seq_len(nrow(d)))
)
fold_labels <- map_dfr(seq_along(folds),function(j)
  tibble(fold=j,sample_index=folds[[j]],
         Site_ID=d$Site_ID[folds[[j]]],
         Depth=d$Depth[folds[[j]]])
)
write_csv(fold_labels,"outputs/joint_DO_19_CV_folds.csv")

cv_preds <- list()
cv_diag <- list()
if (cfg$run_CV) {
  selected_folds <- if (is.null(cfg$run_only_fold))
    seq_along(folds) else cfg$run_only_fold
  stopifnot(all(selected_folds %in% seq_along(folds)))
  for (fold in selected_folds) {
    held <- folds[[fold]]
    for (use_na in 0:1) {
      label <- names_model[use_na+1L]
      message(
        "DO grouped CV fold ",fold,"/",cfg$nfolds,
        " model ",label
      )
      x <- base
      x$has_DO[held] <- 0L
      x$DO_obs[held] <- 0  # DO absent even as data
      fit <- fit_args(use_na,fold,data=x)
      vars <- paste0("mean_DO[",held,"]")
      draws <- posterior::as_draws_matrix(
        fit$draws(variables=c(vars,"sigma_DO"))
      )
      predicted <- colMeans(draws[,vars,drop=FALSE])
      set.seed(cfg$seed+fold+100*use_na)
      pp <- vapply(vars,function(v) {
        draws[,v]+rnorm(nrow(draws),
          sd=sqrt(draws[,"sigma_DO"]^2+.2^2))
      },numeric(nrow(draws)))
      if (is.null(dim(pp))) pp <- matrix(pp,ncol=1)
      cv_preds[[length(cv_preds)+1L]] <- tibble(
        fold=fold,model=label,
        sample_index=held,Site_ID=d$Site_ID[held],
        Depth=d$Depth[held],observed_DO=d$DO[held],
        predicted_DO=as.numeric(predicted),
        predicted_lower95=apply(pp,2,quantile,.025),
        predicted_upper95=apply(pp,2,quantile,.975),
        covered95=observed_DO>=predicted_lower95 &
          observed_DO<=predicted_upper95
      )
      cv_diag[[length(cv_diag)+1L]] <- check_mcmc(
        fit,c(active_global(use_na),
              paste0("f_LSR[",held,"]"),
              paste0("true_MRT[",held,"]")),
        paste0("fold_",fold,"_",label)
      )
      # Output checkpoints survive interrupted runs.
      write_csv(bind_rows(cv_preds),
                "outputs/joint_DO_19_CV_predictions.csv")
      write_csv(bind_rows(cv_diag),
                "outputs/joint_DO_19_CV_diagnostics.csv")
    }
  }
  preds <- bind_rows(cv_preds)
  score <- preds %>%
    group_by(model) %>%
    summarise(
      n=n(),completed_folds=n_distinct(fold),
      RMSE=sqrt(mean((predicted_DO-observed_DO)^2)),
      MAE=mean(abs(predicted_DO-observed_DO)),
      bias=mean(predicted_DO-observed_DO),
      predictive_coverage95=mean(covered95),
      .groups="drop"
    )
  write_csv(score,"outputs/joint_DO_19_CV_scores.csv")
  paired <- preds %>%
    select(sample_index,fold,Site_ID,observed_DO,
           model,predicted_DO) %>%
    pivot_wider(names_from=model,values_from=predicted_DO) %>%
    mutate(
      abs_error_3=abs(
        observed_DO-.data[["Three tracers"]]),
      abs_error_4=abs(
        observed_DO-.data[["Four tracers"]]),
      difference_abs_error=abs_error_4-abs_error_3
    )
  write_csv(paired,"outputs/joint_DO_19_CV_paired_errors.csv")
  print(score,n=Inf)

  g <- ggplot(preds,aes(observed_DO,predicted_DO)) +
    geom_abline(slope=1,intercept=0,
                linetype="dashed",colour="grey60") +
    geom_point(size=2,alpha=.8) +
    coord_equal() +
    facet_wrap(~model) +
    labs(title="Fully joint Bayesian DO: held-out wells",
         subtitle="Both models use identical well groups, DO/MRT priors",
         x="Observed DO (mg/L)",
         y="Posterior mean prediction (mg/L)") +
    theme_bw(base_size=12)
  print(g)
  ggsave("plots/joint_DO_19_grouped_CV_predictions.png",
         g,width=10,height=6,dpi=300)
}
# Posterior 95% intervals, not the default 90% q5/q95.
p <- globs %>%
  filter(variable %in% c("k","LSR_DO_deficit",
                          "LSR_entry_DO",
                          "mean_well_LSR_pathway_loss",
                          "mean_well_age_associated_loss")) %>%
  ggplot(aes(median,model))+
  geom_segment(aes(x=lower95,xend=upper95,
                   y=model,yend=model))+
  geom_point(size=2.5)+
  facet_wrap(~variable,scales="free_x")+
  labs(title="Joint oxygen depletion: does Na change the posterior?",
       subtitle="Median and central 95% credible interval",
       x="Parameter value",y=NULL)+
  theme_bw(base_size=11)
print(p)
ggsave("plots/joint_DO_19_posterior_comparison.png",
       p,width=12,height=7,dpi=300)

cat("\nGLOBAL POSTERIOR COMPARISON:\n")
print(globs %>% filter(variable %in%
  c("k","LSR_DO_deficit","LSR_entry_DO","river_DO",
    "mean_well_LSR_pathway_loss",
    "mean_well_age_associated_loss",
    "sampled_well_fraction_LSR_of_loss")),n=Inf)
cat("\nFRACTION CHANGE:\n")
print(mix_shift %>% summarise(
  n=n(),
  median_absolute_change=median(
    abs(fraction_shift_Four_minus_Three)),
  max_absolute_change=max(
    abs(fraction_shift_Four_minus_Three))
))
cat(
  "\nINTERPRETATION:\n",
  "* The only added observed tracer is Na, with Na-Cl dependence.\n",
  "* Temporal LSR isotope variation is present in BOTH models.\n",
  "* DO and MRT are modelled jointly, not inserted posterior means.\n",
  "* Validation masks held-out DO completely; tracers remain.\n",
  "* 9.5 mg/L pre-LSR source DO is assumed, not measured.\n",
  "* The age-dependent loss is not proven saturated-aquifer loss.\n",
  "* The arithmetic across-wells loss fraction is NOT volume weighted.\n",
  "* Endmember centres are partly from study wells (empirical Bayes).\n",
  "* Report only after reviewing Rhat, ESS, HMC divergences and\n",
  "  held-out predictive coverage. Data remain LOCAL AND PRIVATE.\n",
  sep=""
)
