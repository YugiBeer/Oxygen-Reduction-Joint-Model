# 20_joint_paired_DO_excessN2.R
#
# Paired dissolved oxygen (DO) and excess-N2 Bayesian prediction.
#
# Four mixing tracers: delta18O, Cl, Na, excess air.
# Two recharge sources: river recharge and land surface recharge (LSR).
# Na-Cl correlated observations, LSR isotope-endmember variability,
# mean-centred uncertain true MRT, uncertain source entry DO, common k.
#
# Compare:
# - DO_only: NO excess N2 in fit; shared DO model reference.
# - N2_age: additionally model excess N2 using source + ln(1+MRT).
# - N2_DOloss: add beta * apparent cumulative oxygen depletion.
# - N2_lowDO: add beta * time exposed below low_DO_threshold (mg/L).
#
# The N2 source-associated intercepts are empirical statistical terms,
# NOT conservative N2 mixing endmembers. Beta parameters are NOT
# stoichiometric yields. Nitrate availability is unobserved at entry.
# Excess N2 is NEVER used as a conservative source-mixing tracer.
# Derived initial_NO3N and reaction_progress MUST NEVER be predictors
# because they already contain the measured excess-N2 quantity.
#
# Both DO and N2 observations are held out together in grouped CV;
# conservative tracers and reported MRT remain available. Matched
# groups and same data for every model; compare paired DO and N2 RMSE.
#
# Default 4 full + 20 group CV fits; expensive. Start run_CV=FALSE
# or run_only_fold=1L if needed. Refit with exclude_16758=TRUE to
# see whether the high excess-N2 well drives all inferred reaction.
#
# Local original input: data/Uwe_HP_data_prepared.csv.
# Outputs local in gitignored outputs/ and plots/ directories.
library(tidyverse)
if (file.exists("stan/20_joint_paired_DO_excessN2.stan")) {
  # repository root
} else if (file.exists(
  "Oxygen-Reduction-Joint-Model/stan/20_joint_paired_DO_excessN2.stan"
)) {
  setwd("Oxygen-Reduction-Joint-Model")
} else stop("Open Oxygen-Reduction-Joint-Model as the working directory.")

cfg <- list(
  input = "data/Uwe_HP_data_prepared.csv",
  stan_model = "stan/20_joint_paired_DO_excessN2.stan",
  seed = 20261031L,
  log_mrt_SD = 0.25,           # assumed, not measured
  endmember_prior_scale = 1,
  lsr_do_prior_a = 1, lsr_do_prior_b = 1,
  low_DO_threshold = 2,        # mg/L, sensitivity: 1, 2, 3
  exclude_16758 = FALSE,      # recommended sensitivity rerun
  nfolds = 5L,
  run_CV = TRUE,
  run_only_fold = NULL,
  run_DOloss = TRUE,           # omit to save 6 Stan fits
  chains = 4L, parallel_chains = 4L,
  warmup = 1400L, sampling = 1400L,
  adapt_delta = .98, max_treedepth = 13L,
  max_rhat = 1.05, min_ess_bulk = 400,
  save_fits = TRUE
)
if (!requireNamespace("cmdstanr",quietly=TRUE) ||
    !requireNamespace("posterior",quietly=TRUE))
  stop("Install cmdstanr, CmdStan, and posterior first.")
if (!file.exists(cfg$input)) stop(
  "Missing data/Uwe_HP_data_prepared.csv. Copy the ORIGINAL private ",
  "prepared dataset from Lincoln-Environmental/Oxygen-Reduction. ",
  "The reduced groundwater.csv lacks sodium and site 16758."
)
tryCatch(cmdstanr::cmdstan_path(),error=function(e)
  stop("CmdStan not configured: ",conditionMessage(e))
)
dir.create("outputs",recursive=TRUE,showWarnings=FALSE)
dir.create("plots",recursive=TRUE,showWarnings=FALSE)

d <- read_csv(
  cfg$input, show_col_types=FALSE,
  col_types=cols(Site_ID=col_character(),.default=col_guess())
)
required <- c(
  "Site_ID","Depth","DO","MRT","O18","Cl","Na","xs_Air",
  "xs_Air_uncertainty","xs_N2","xs_N2_sig"
)
if (!all(required %in% names(d))) stop(
  "Missing fields: ",
  paste(setdiff(required,names(d)),collapse=", ")
)
d <- d %>%
  filter(
    is.finite(DO),DO>=0,
    is.finite(xs_N2),xs_N2>=0,
    is.finite(MRT),MRT>=0
  ) %>%
  { if (isTRUE(cfg$exclude_16758))
      filter(.,Site_ID!="16758") else . } %>%
  mutate(
    sample_index=row_number(),
    group=case_when(
      str_detect(Site_ID,"Eastbourne")~"Eastbourne",
      str_detect(Site_ID,"Gilbertson")~"Gilbertson",
      str_detect(Site_ID,"Flaxmere")~"Flaxmere",
      TRUE~Site_ID
    )
  )
stopifnot(nrow(d)>=30)
if (any(d$Site_ID=="16758")) message(
  "High-N2 well 16758 included: excess N2 = ",
  d$xs_N2[d$Site_ID=="16758"][1],
  ", MRT = ",d$MRT[d$Site_ID=="16758"][1],
  " years. Refit with cfg$exclude_16758=TRUE as sensitivity."
)
message("Joint paired DO + excess-N2 analysis: n=",nrow(d))

obs <- function(x) ifelse(is.finite(x),as.numeric(x),0)
sd_default <- function(x,v) ifelse(is.finite(x)&x>0,x,v)
base <- list(
  N=nrow(d),
  has_O18=as.integer(is.finite(d$O18)),
  has_Cl=as.integer(is.finite(d$Cl)),
  has_Na=as.integer(is.finite(d$Na)),
  has_air=as.integer(is.finite(d$xs_Air)),
  has_DO=rep(1L,nrow(d)),
  has_N2=rep(1L,nrow(d)),
  O18=obs(d$O18), Cl=obs(d$Cl),
  Na=obs(d$Na), air=obs(d$xs_Air),
  DO_obs=as.numeric(d$DO),
  xsN2_obs=as.numeric(d$xs_N2),
  Cl_sd=pmax(1,.03*obs(d$Cl)),
  Na_sd=pmax(.5,.03*obs(d$Na)),
  air_sd=sd_default(d$xs_Air_uncertainty,.15),
  xsN2_sd=sd_default(d$xs_N2_sig,.3),
  mrt_reported=as.numeric(d$MRT),
  use_na=1L,
  use_n2=0L,n2_link=0L,
  endmember_prior_scale=cfg$endmember_prior_scale,
  log_mrt_SD=cfg$log_mrt_SD,
  lsr_do_prior_a=cfg$lsr_do_prior_a,
  lsr_do_prior_b=cfg$lsr_do_prior_b,
  low_DO_threshold=cfg$low_DO_threshold
)

models <- tribble(
  ~label,             ~use_n2, ~n2_link,
  "DO_only",            0L,       0L,
  "N2_age",             1L,       0L,
  "N2_DOloss",          1L,       1L,
  "N2_lowDO",           1L,       2L
)
if (!cfg$run_DOloss)
  models <- models %>% filter(label!="N2_DOloss")

active_vars <- function(use_n2,link) c(
  "k","LSR_entry_DO","LSR_DO_deficit","river_DO","sigma_DO",
  "mean_well_LSR_pathway_loss",
  "mean_well_age_associated_loss",
  "sampled_well_fraction_LSR_of_loss",
  "mean_well_LSR_fraction",
  "tau_LSR_O18","rho_Na_Cl",
  if (use_n2==1L) c(
    "n2_R","n2_LSR","beta_n2_age","sigma_N2",
    if (link>0) "beta_n2_link"
  )
)

model <- cmdstanr::cmdstan_model(cfg$stan_model,quiet=TRUE)
fit_model <- function(spec,fold=0L,stan_data=base) {
  stan_data$use_n2 <- as.integer(spec$use_n2)
  stan_data$n2_link <- as.integer(spec$n2_link)
  if (spec$use_n2==0L) {
    stan_data$has_N2[] <- 0L
    stan_data$xsN2_obs[] <- 0
  }
  model$sample(
    data=stan_data,
    seed=cfg$seed+1000L*fold+101L*match(spec$label,models$label),
    chains=cfg$chains,parallel_chains=cfg$parallel_chains,
    iter_warmup=cfg$warmup,iter_sampling=cfg$sampling,
    adapt_delta=cfg$adapt_delta,max_treedepth=cfg$max_treedepth,
    refresh=300
  )
}
summary95 <- function(fit,vars) {
  x <- posterior::as_draws_matrix(fit$draws(variables=vars))
  tibble(
    variable=colnames(x),
    median=apply(x,2,median),
    lower95=apply(x,2,quantile,probs=.025),
    upper95=apply(x,2,quantile,probs=.975)
  )
}
diagnostic <- function(fit,spec,fold=0L,held=NULL) {
  vars <- active_vars(spec$use_n2,spec$n2_link)
  vars <- c(
    vars,
    if (is.null(held)) "f_LSR" else
      paste0("f_LSR[",held,"]"),
    if (is.null(held)) "true_MRT" else
      paste0("true_MRT[",held,"]")
  )
  x <- fit$summary(variables=vars)
  hmc <- fit$diagnostic_summary()
  if (any(x$rhat>cfg$max_rhat,na.rm=TRUE) ||
      any(x$ess_bulk<cfg$min_ess_bulk,na.rm=TRUE) ||
      any(hmc$num_divergent>0,na.rm=TRUE))
    warning(
      "MCMC diagnostics need review: ",spec$label,
      ", fold ",fold
    )
  x %>% mutate(model=spec$label,fold=fold,.before=1)
}
suffix <- paste0(
  "_N",nrow(d),
  if (cfg$exclude_16758) "_without16758" else "",
  "_threshold",cfg$low_DO_threshold
)
to_file <- function(prefix) paste0("outputs/",prefix,suffix,".csv")

full <- list();conv <- list();by_well <- list()
for (j in seq_len(nrow(models))) {
  s <- models[j,]
  message("Full paired Bayesian fit ",s$label)
  fit <- fit_model(s)
  full[[j]] <- summary95(
    fit,active_vars(s$use_n2,s$n2_link)
  ) %>% mutate(model=s$label,.before=1)
  conv[[j]] <- diagnostic(fit,s)
  by_well[[j]] <- fit$summary(
    variables=c("f_LSR","mean_DO",
                "mean_N2","true_MRT",
                "LSR_pathway_DO_loss",
                "age_associated_DO_loss",
                "low_DO_duration_years")
  ) %>%
    mutate(
      sample_index=as.integer(str_match(
        variable,"\\[(\\d+)\\]")[,2]),
      model=s$label,.before=1
    ) %>%
    inner_join(
      d %>% select(
        sample_index,Site_ID,Depth,DO,xs_N2,MRT
      ),by="sample_index"
    )
  if (cfg$save_fits) fit$save_object(file=paste0(
    "outputs/joint_DO_N2_20_",s$label,suffix,".rds"
  ))
  write_csv(bind_rows(full),
            to_file("joint_DO_N2_20_full_posteriors"))
  write_csv(bind_rows(conv),
            to_file("joint_DO_N2_20_convergence"))
  write_csv(bind_rows(by_well),
            to_file("joint_DO_N2_20_per_well"))
}

# Paired, grouped five-fold CV. DO and xsN2 are BOTH zeroed for
# held-out wells (regardless of the link). Independent tracers and
# MRT remain; this is intentionally not one-outcome-at-a-time CV.
set.seed(cfg$seed)
g <- split(seq_len(nrow(d)),d$group)
order <- sample(seq_along(g))
order <- order[base::order(-lengths(g)[order])]
folds <- vector("list",cfg$nfolds)
for (i in order) {
  j <- which.min(lengths(folds))
  folds[[j]] <- c(folds[[j]],g[[i]])
}
stopifnot(identical(sort(unlist(folds)),seq_len(nrow(d))))
fold_map <- map_dfr(seq_along(folds),function(i) tibble(
  fold=i,sample_index=folds[[i]],
  Site_ID=d$Site_ID[folds[[i]]],Depth=d$Depth[folds[[i]]]
))
write_csv(fold_map,to_file("joint_DO_N2_20_CV_folds"))

if (cfg$run_CV) {
  selected <- if (is.null(cfg$run_only_fold))
    seq_along(folds) else cfg$run_only_fold
  stopifnot(all(selected %in% seq_along(folds)))
  preds <- list();diags <- list()
  for (fold in selected) {
    held <- folds[[fold]]
    for (j in seq_len(nrow(models))) {
      s <- models[j,]
      x <- base
      x$has_DO[held] <- 0L
      x$DO_obs[held] <- 0
      x$has_N2[held] <- 0L
      x$xsN2_obs[held] <- 0
      message(
        "Paired DO+N2 holdout fold ",fold,
        " model ",s$label
      )
      fit <- fit_model(s,fold,stan_data=x)
      var_do <- paste0("mean_DO[",held,"]")
      var_n2 <- paste0("mean_N2[",held,"]")
      var_rep_do <- paste0("replicated_DO[",held,"]")
      var_rep_n2 <- paste0("replicated_N2[",held,"]")
      m <- posterior::as_draws_matrix(fit$draws(
        variables=c(var_do,var_n2,
                    var_rep_do,var_rep_n2)
      ))
      summarise_var <- function(nms,fn) vapply(nms,function(v)
        fn(m[,v]),numeric(1))
      q <- function(p) function(v) unname(quantile(v,p))
      pp <- tibble(
        model=s$label,fold=fold,
        sample_index=held,
        Site_ID=d$Site_ID[held],Depth=d$Depth[held],
        observed_DO=d$DO[held],
        predicted_DO=summarise_var(var_do,mean),
        DO_lower95=summarise_var(var_rep_do,q(.025)),
        DO_upper95=summarise_var(var_rep_do,q(.975)),
        observed_N2=d$xs_N2[held],
        # Do NOT score N2 from the DO-only model; N2 parameters
        # are prior-only there, hence no learnt N2 prediction.
        predicted_N2=if (s$use_n2==1L)
          summarise_var(var_n2,mean) else NA_real_,
        N2_lower95=if (s$use_n2==1L)
          summarise_var(var_rep_n2,q(.025)) else NA_real_,
        N2_upper95=if (s$use_n2==1L)
          summarise_var(var_rep_n2,q(.975)) else NA_real_
      ) %>%
      mutate(
        covered_DO_95=observed_DO>=DO_lower95 &
                      observed_DO<=DO_upper95,
        covered_N2_95=observed_N2>=N2_lower95 &
                      observed_N2<=N2_upper95
      )
      preds[[length(preds)+1L]] <- pp
      diags[[length(diags)+1L]] <- diagnostic(
        fit,s,fold,held
      )
      write_csv(bind_rows(preds),
                to_file("joint_DO_N2_20_CV_predictions"))
      write_csv(bind_rows(diags),
                to_file("joint_DO_N2_20_CV_diagnostics"))
    }
  }
  pred <- bind_rows(preds)
  cv_score <- pred %>%
    group_by(model) %>%
    summarise(
      n=n(),nfolds=n_distinct(fold),
      DO_RMSE=sqrt(mean(
        (predicted_DO-observed_DO)^2
      )),
      DO_MAE=mean(abs(predicted_DO-observed_DO)),
      DO_coverage95=mean(covered_DO_95),
      N2_RMSE=if (all(is.na(predicted_N2)))
        NA_real_ else sqrt(mean(
          (predicted_N2-observed_N2)^2
        )),
      N2_MAE=if (all(is.na(predicted_N2)))
        NA_real_ else mean(abs(predicted_N2-observed_N2)),
      N2_coverage95=if (all(is.na(covered_N2_95)))
        NA_real_ else mean(covered_N2_95),
      .groups="drop"
    )
  write_csv(cv_score,to_file("joint_DO_N2_20_CV_summary"))
  print(cv_score,n=Inf)

  g1 <- pred %>%
    filter(model!="DO_only") %>%
    select(model,observed_DO,predicted_DO,
           observed_N2,predicted_N2) %>%
    pivot_longer(
      cols=c(observed_DO,predicted_DO,
             observed_N2,predicted_N2),
      names_to=c(".value","analyte"),
      names_pattern="(observed|predicted)_(DO|N2)"
    ) %>%
    ggplot(aes(observed,predicted)) +
    geom_abline(slope=1,intercept=0,
                linetype="dashed",colour="grey60") +
    geom_point(alpha=.8,size=1.8) +
    facet_grid(model~analyte,scales="free") +
    labs(title="Joint Bayesian prediction of DO and excess N2",
         subtitle="Both outcomes withheld at the same wells",
         x="Measured",y="Posterior mean prediction")+
    theme_bw(base_size=11)
  print(g1)
  ggsave(paste0("plots/joint_DO_N2_20_CV",suffix,".png"),
         g1,width=12,height=8,dpi=300)
}

global <- bind_rows(full)
p <- global %>%
  filter(variable %in% c(
    "k","LSR_DO_deficit","river_DO",
    "mean_well_LSR_pathway_loss",
    "mean_well_age_associated_loss",
    "beta_n2_link","sigma_N2"
  )) %>%
  ggplot(aes(median,model)) +
  geom_segment(aes(
    x=lower95,xend=upper95,y=model,yend=model
  )) +
  geom_point(size=2.2)+
  facet_wrap(~variable,scales="free_x")+
  labs(
    title="Joint DO + excess-N2: model and prior sensitivity",
    subtitle="Medians and central 95% credible intervals",
    x=NULL,y=NULL
  ) + theme_bw(base_size=11)
print(p)
ggsave(paste0("plots/joint_DO_N2_20_posterior",suffix,".png"),
       p,width=12,height=7,dpi=300)

cat("\nFULL POSTERIORS:\n")
print(global %>% filter(variable %in% c(
  "k","LSR_DO_deficit","river_DO",
  "mean_well_LSR_pathway_loss",
  "mean_well_age_associated_loss",
  "beta_n2_link","sigma_N2"
)),n=Inf)
cat(
  "\nSCIENCE AND VALIDATION NOTES:\n",
  "* Excess N2 is a reactive outcome, NOT an O2-loss measurement.\n",
  "* Denitrification consumes nitrate as electron acceptor; it\n",
  "  generally follows O2 depletion rather than consuming O2.\n",
  "* Source-associated N2 terms are empirical backgrounds.\n",
  "* Initial_NO3N and reaction_progress contain excess N2 and\n",
  "  MUST NOT be used to predict held-out excess N2.\n",
  "* DO threshold exposure is a simplified piston-flow history.\n",
  "* Interpret N2 link coefficients statistically, not chemically.\n",
  "* All four scenarios use the same conservative tracers and MRT.\n",
  "* Compare with and without site 16758 and threshold=1,2,3.\n",
  "* The 9.5 mg/L LSR oxygen reference is assumed, not observed.\n",
  "* Age-associated DO loss is not proof of saturated-aquifer loss.\n",
  "* The mean across wells is NOT a recharge-volume-weighted budget.\n",
  "* Results require acceptable Rhat, ESS, and HMC divergences.\n",
  "* Data, predictions, and posterior draws are PRIVATE and local.\n",
  sep=""
)
