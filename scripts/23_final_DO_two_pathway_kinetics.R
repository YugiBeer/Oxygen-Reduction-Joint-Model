# 23_final_DO_two_pathway_kinetics.R
#
# FINAL PRIMARY BAYESIAN MODEL: aquifer zero/first-order DO kinetics
# and pre-water-table LSR oxygen deficit from assumed saturated rainfall.
#
# PHYSICAL BOUNDARY CONDITIONS (explicit assumptions):
#   Rain at land surface: DO = Csat(T_recharge, P_atm), saturated.
#   LSR at water table: DO = Csat * survival_LSR; 0 < survival < 1.
#   River recharge at aquifer entry: DO=7.5 mg/L FIXED.
#   River-dominated well 2801 is set to 7.5 in model-only data,
#   never in the original source file or the fitting likelihood.
#   Aquifer: both endmembers attenuate for reported MRT independently,
#     then mix. Same aquifer reaction coefficient for both sources.
#     Zero and first order are COMPETING MODELS, not fitted together.
#
# CRUCIAL: Without vadose-zone travel time we estimate LSR DO deficit,
# NOT a separate measured LSR zero- or first-order reaction rate.
# Conditional LSR rates for ASSUMED t_vadose are exported solely
# to show nonidentifiability; these are reparameterisations of one
# inferred entry concentration.
#
# DO-only reactive response (progress, xs_N2 not fitted), mixing from
# O18, Na and Cl (joint Na-Cl residuals), optional excess air ablation.
# Variable LSR isotope source uncertainty and true groundwater MRT.
# Exclude anomalous far-south site 16758 from all analyses.
# The model-only 2801 DO is 7.5 but is ALWAYS masked in the
# likelihood and NEVER CV-scored (no pseudo-observation).
#
# Needs CmdStanR, CmdStan, tidyverse and posterior. Data PRIVATE.
library(tidyverse)

if (file.exists("stan/23_final_DO_two_pathway_kinetics.stan")) {
  # Running from repository root.
} else if (file.exists(
  "Oxygen-Reduction-Joint-Model/stan/23_final_DO_two_pathway_kinetics.stan"
)) {
  setwd("Oxygen-Reduction-Joint-Model")
} else {
  stop("Set WD to Oxygen-Reduction-Joint-Model repository root.")
}

cfg <- list(
  input="data/Uwe_HP_data_prepared.csv",
  stan="stan/23_final_DO_two_pathway_kinetics.stan",
  excluded_site="16758",
  river_proxy_site="2801",     # river-dominated well; NOT a river sample
  river_DO_boundary=7.5,       # fixed river-entry DO, mg/L
  set_2801_model_DO=7.5,       # local modelling value only
  mask_river_proxy_DO=TRUE,      # mandatory to avoid double-use
  recharge_temperature_C=15,   # assumed; vary 10/20 C in sensitivities
  atmospheric_pressure_ratio=1,# P/P_sea_level, approximate
  C_ref=7.5,                    # mg/L fractional-kinetics comparison
  survival_prior_a=1,          # Beta on LSR DO surviving to water table
  survival_prior_b=1,
  endmember_prior_scale=1,     # empirical endpoints - test 2
  log_mrt_SD=.25,              # assume reported MRT mostly aquifer age
  vadose_years_assumptions=c(.25,1,5),
  include_excess_air=FALSE,    # previous held-out tests show little gain
  folds=5L,
  seed=20261031L,
  chains=4L,
  parallel_chains=4L,
  warmup=1400L,
  sampling=1400L,
  adapt_delta=.98,
  max_treedepth=13L,
  max_rhat=1.05,
  min_bulk_ESS=400,
  run_CV=TRUE,                # 2 full and 10 paired grouped CV fits
  run_only_fold=NULL,        # 1L for smoke test, else NULL
  run_prior_sensitivity=FALSE, # full fits for selected BC/priors
  save_fits=TRUE
)

if (!cfg$mask_river_proxy_DO)
  stop("Proxy site DO must stay masked when used as river prior.")
if (!requireNamespace("cmdstanr",quietly=TRUE) ||
    !requireNamespace("posterior",quietly=TRUE))
  stop("Requires cmdstanr, installed CmdStan and posterior.")
tryCatch(cmdstanr::cmdstan_path(),
         error=function(e) stop("CmdStan unavailable: ",conditionMessage(e)))
if (!file.exists(cfg$input))
  stop("Copy the ORIGINAL PRIVATE prepared CSV locally to ",cfg$input)
stopifnot(
  cfg$recharge_temperature_C>=0,cfg$recharge_temperature_C<=30,
  cfg$atmospheric_pressure_ratio>0,
  cfg$C_ref>0,cfg$river_DO_boundary>0,cfg$log_mrt_SD>0,
  all(cfg$vadose_years_assumptions>0),
  abs(cfg$set_2801_model_DO-cfg$river_DO_boundary)<1e-12
)

# Freshwater equilibrium oxygen saturation in mg/L; polynomial
# approximation for 0-30 degrees C at ~1 atm. Pressure-ratio
# correction is APPROXIMATE (ignores water vapour correction);
# use external measured T, atmospheric pressure if available.
do_sat_mgL <- function(T,pressure_ratio=1) {
  (14.652-0.41022*T+0.007991*T^2-
     0.000077774*T^3)*pressure_ratio
}
Csat <- do_sat_mgL(
  cfg$recharge_temperature_C,cfg$atmospheric_pressure_ratio
)
message(
  "FIXED SATURATED RAINFALL BC: T=",
  cfg$recharge_temperature_C," C, pressure ratio=",
  cfg$atmospheric_pressure_ratio,", Csat=",
  signif(Csat,4)," mg/L."
)
message(
  "FIXED RIVER-RECHARGE BC: ",cfg$river_DO_boundary,
  " mg/L; model-only site 2801 DO assigned ",
  cfg$set_2801_model_DO,
  " mg/L, not used as an independent observation."
)
dir.create("outputs",recursive=TRUE,showWarnings=FALSE)
dir.create("plots",recursive=TRUE,showWarnings=FALSE)
# Explicit BC audit: these are assumptions or posterior variables,
# never independent measurement evidence.
write_csv(
  tibble(
    boundary=c("rainfall_surface","river_aquifer_entry",
               "LSR_water_table","aquifer"),
    DO_mgL=c(Csat,cfg$river_DO_boundary,NA_real_,NA_real_),
    source=c("fixed_saturation_assumption",
             "fixed_7p5_mgL",
             "posterior_uncertain",
             "posterior_reaction_coefficient"),
    note=c(
      paste0("T_recharge=",cfg$recharge_temperature_C,
             " deg C, pressure ratio=",
             cfg$atmospheric_pressure_ratio),
      "2801 model-only DO is 7.5 and excluded from fitting",
      "LSR at water table is Csat multiplied by posterior survival",
      paste0("Zero/first-order k_ref prior common at ",
             cfg$C_ref," mg/L")
    )
  ),
  "outputs/final_DO_23_boundary_conditions.csv"
)

raw <- readr::read_csv(
  cfg$input,show_col_types=FALSE,
  col_types=cols(Site_ID=col_character(),.default=col_guess())
)
required <- c(
  "Site_ID","Depth","DO","MRT","O18","Cl","Na",
  "xs_Air","xs_Air_uncertainty"
)
if (!all(required %in% names(raw)))
  stop("Missing required columns: ",
       paste(setdiff(required,names(raw)),collapse=", "))
stopifnot(
  nrow(raw)>=40,
  any(raw$Site_ID==cfg$excluded_site),
  sum(raw$Site_ID==cfg$river_proxy_site)==1L
)
d <- raw %>%
  filter(
    Site_ID!=cfg$excluded_site,
    is.finite(DO),DO>=0,
    is.finite(MRT),MRT>=0
  ) %>%
  mutate(
    sample_index=row_number(),
    cv_group=case_when(
      str_detect(Site_ID,"Eastbourne")~"Eastbourne",
      str_detect(Site_ID,"Gilbertson")~"Gilbertson",
      str_detect(Site_ID,"Flaxmere")~"Flaxmere",
      TRUE~Site_ID
    )
  )
stopifnot(nrow(d)>=35L)
proxy_index <- which(d$Site_ID==cfg$river_proxy_site)
stopifnot(length(proxy_index)==1L)
# Never alter the private input CSV or substitute values in-place.
d$DO_original <- d$DO
message(
  "River proxy 2801 original DO=",d$DO_original[proxy_index],
  " mg/L, model-only assigned value=",cfg$set_2801_model_DO,
  " mg/L (DO observation masked)."
)
d$DO[proxy_index] <- cfg$set_2801_model_DO

as_obs <- function(x) ifelse(is.finite(x),as.numeric(x),0)
as_sigma <- function(x,default) ifelse(
  is.finite(x)&x>0,as.numeric(x),default
)
stan_data <- list(
  N=nrow(d),
  has_O18=as.integer(is.finite(d$O18)),
  has_Cl=as.integer(is.finite(d$Cl)),
  has_Na=as.integer(is.finite(d$Na)),
  has_air=if (cfg$include_excess_air)
    as.integer(is.finite(d$xs_Air)) else rep(0L,nrow(d)),
  has_DO=rep(1L,nrow(d)),
  O18=as_obs(d$O18),
  Cl=as_obs(d$Cl),
  Na=as_obs(d$Na),
  air=if (cfg$include_excess_air)
    as_obs(d$xs_Air) else rep(0,nrow(d)),
  DO_obs=as.numeric(d$DO),
  Cl_sd=pmax(1,.03*as_obs(d$Cl)),
  Na_sd=pmax(.5,.03*as_obs(d$Na)),
  air_sd=as_sigma(d$xs_Air_uncertainty,.15),
  mrt_reported=as.numeric(d$MRT),
  aq_order=0L,
  Csat_rain=Csat,
  C_ref=cfg$C_ref,
  endmember_prior_scale=cfg$endmember_prior_scale,
  log_mrt_SD=cfg$log_mrt_SD,
  river_DO_boundary=cfg$river_DO_boundary,
  beta_survival_a=cfg$survival_prior_a,
  beta_survival_b=cfg$survival_prior_b
)
# Site 2801 is a model-only source boundary value, NOT evidence.
stan_data$has_DO[proxy_index] <- 0L
stan_data$DO_obs[proxy_index] <- cfg$set_2801_model_DO
eligible <- which(stan_data$has_DO==1L)
message(
  "After removing ",cfg$excluded_site,": ",
  nrow(d)," wells with DO+MRT and chemistry; ",
  length(eligible)," independent DO responses after reserving ",
  cfg$river_proxy_site," as boundary proxy."
)

# Full posterior reports are central 95% intervals, NOT q5/q95.
posterior95 <- function(fit,vars) {
  sm <- fit$summary(variables=vars)
  mat <- posterior::as_draws_matrix(fit$draws(variables=vars))
  sm %>%
    select(variable,mean,median,sd,rhat,ess_bulk,ess_tail) %>%
    mutate(
      lower95=vapply(variable,function(v)
        unname(quantile(mat[,v],.025)),numeric(1)),
      upper95=vapply(variable,function(v)
        unname(quantile(mat[,v],.975)),numeric(1))
    )
}
names_global <- c(
  "river_DO","survival_LSR","LSR_DO_at_water_table",
  "LSR_DO_deficit","k_ref","sigma_DO",
  "aquifer_rate_actual",
  "LSR_zero_order_exposure","LSR_first_order_exposure",
  "mean_vadose_LSR_DO_loss","mean_aquifer_DO_loss",
  "sampled_well_vadose_loss_fraction",
  "mean_f_LSR","O18_R","O18_LSR","Cl_R","Cl_LSR",
  "Na_R","Na_LSR","tau_LSR_O18","rho_Na_Cl"
)
select_diag <- function(fit,order,fold=0L,held=NULL) {
  extra <- if (is.null(held)) c("f_LSR","true_MRT")
    else c(paste0("f_LSR[",held,"]"),
           paste0("true_MRT[",held,"]"))
  z <- fit$summary(variables=c(names_global,extra))
  hm <- fit$diagnostic_summary()
  bad <- any(z$rhat>cfg$max_rhat,na.rm=TRUE) ||
    any(z$ess_bulk<cfg$min_bulk_ESS,na.rm=TRUE) ||
    any(hm$num_divergent>0,na.rm=TRUE)
  if (bad) warning(
    "MCMC diagnostics not acceptable: order=",order,
    " fold=",fold,
    " (review Rhat, ESS, divergences before inference)."
  )
  z %>%
    mutate(
      aquifer_order=order,fold=fold,
      any_divergence=any(hm$num_divergent>0),
      .before=1
    )
}
mod <- cmdstanr::cmdstan_model(cfg$stan,quiet=TRUE)
sample_model <- function(order,x=stan_data,fold=0L,scenario=0L) {
  x$aq_order <- as.integer(order)
  mod$sample(
    data=x,
    seed=cfg$seed+10000L*fold+101L*order+23001L*scenario,
    chains=cfg$chains,parallel_chains=cfg$parallel_chains,
    iter_warmup=cfg$warmup,iter_sampling=cfg$sampling,
    adapt_delta=cfg$adapt_delta,
    max_treedepth=cfg$max_treedepth,
    refresh=250
  )
}
post <- list();diag <- list();wells <- list();conditional <- list()
for (order in 0:1) {
  message("FULL JOINT FIT aquifer reaction ORDER ",order,
          if (order==0) " (zero-order)" else " (first-order)")
  fit <- sample_model(order)
  post[[order+1L]] <- posterior95(fit,names_global) %>%
    mutate(
      aquifer_order=order,
      kinetic_rate_units=if_else(
        variable=="aquifer_rate_actual",
        if (order==0) "mg/L/year" else "year^-1",
        NA_character_
      ),
      .before=1
    )
  diag[[order+1L]] <- select_diag(fit,order)
  w <- fit$summary(variables=c(
    "f_LSR","true_MRT","mean_DO",
    "vadose_LSR_DO_loss","aquifer_DO_loss",
    "river_DO_after_aquifer","LSR_DO_after_aquifer"
  ))
  wells[[order+1L]] <- w %>%
    mutate(
      sample_index=as.integer(str_match(
        variable,"\\[(\\d+)\\]")[,2]),
      aquifer_order=order,.before=1
    ) %>%
    inner_join(
      d %>% select(sample_index,Site_ID,Depth,DO,MRT),
      by="sample_index"
    )
  # Given T_vadose, convert POSTERIOR survival into one of two
  # conditional kinetic rates. These are NOT separately observed
  # LSR reaction rates and do not distinguish LSR reaction orders.
  sub <- posterior::as_draws_matrix(fit$draws(
    variables=c("LSR_zero_order_exposure",
                "LSR_first_order_exposure")
  ))
  rates <- map_dfr(cfg$vadose_years_assumptions,function(t) {
    k0 <- sub[,"LSR_zero_order_exposure"]/t
    k1 <- sub[,"LSR_first_order_exposure"]/t
    bind_rows(
      tibble(
        assumed_vadose_years=t,
        LSR_kinetic_order=0L,
        rate_units="mg/L/year",
        posterior_draw=k0
      ),
      tibble(
        assumed_vadose_years=t,
        LSR_kinetic_order=1L,
        rate_units="year^-1",
        posterior_draw=k1
      )
    )
  }) %>%
    group_by(assumed_vadose_years,LSR_kinetic_order,rate_units) %>%
    summarise(
      median=median(posterior_draw),
      lower95=unname(quantile(posterior_draw,.025)),
      upper95=unname(quantile(posterior_draw,.975)),
      .groups="drop"
    ) %>%
    mutate(aquifer_order=order,.before=1)
  conditional[[order+1L]] <- rates
  if (cfg$save_fits) fit$save_object(file=paste0(
    "outputs/final_DO_23_full_order",order,".rds"
  ))
  write_csv(bind_rows(post),"outputs/final_DO_23_global_posteriors.csv")
  write_csv(bind_rows(diag),"outputs/final_DO_23_full_diagnostics.csv")
  write_csv(bind_rows(wells),"outputs/final_DO_23_well_posteriors.csv")
  write_csv(bind_rows(conditional),
            "outputs/final_DO_23_conditional_LSR_rates.csv")
}
main <- bind_rows(post)
write_csv(main,"outputs/final_DO_23_global_posteriors.csv")

# Five GROUPED folds; 2801 is only a source/boundary proxy, never
# a held-out DO target. Group labels and fold membership identical
# between zero- and first-order aquifer models.
set.seed(cfg$seed)
g <- split(eligible,d$cv_group[eligible])
order_g <- sample(seq_along(g))
order_g <- order_g[base::order(-lengths(g)[order_g])]
folds <- vector("list",cfg$folds)
for (i in order_g) {
  j <- which.min(lengths(folds))
  folds[[j]] <- c(folds[[j]],g[[i]])
}
stopifnot(
  all(lengths(folds)>0),
  identical(sort(unlist(folds)),eligible)
)
foldmap <- map_dfr(seq_along(folds),function(j)
  tibble(
    fold=j,sample_index=folds[[j]],
    Site_ID=d$Site_ID[folds[[j]]],
    Depth=d$Depth[folds[[j]]]
  )
)
write_csv(foldmap,"outputs/final_DO_23_cv_folds.csv")

cvpred <- list();cvdiag <- list()
if (cfg$run_CV) {
  chosen <- if (is.null(cfg$run_only_fold))
    seq_along(folds) else cfg$run_only_fold
  stopifnot(all(chosen %in% seq_along(folds)))
  for (fold in chosen) {
    held <- folds[[fold]]
    for (order in 0:1) {
      message(
        "Grouped held-out DO: fold ",fold,
        "/",cfg$folds," aquifer order ",order
      )
      x <- stan_data
      x$has_DO[held] <- 0L
      x$DO_obs[held] <- 0
      fit <- sample_model(order,x=x,fold=fold)
      v <- paste0("mean_DO[",held,"]")
      pr <- paste0("replicated_DO[",held,"]")
      draws <- posterior::as_draws_matrix(
        fit$draws(variables=c(v,pr))
      )
      predmean <- vapply(v,function(z)
        mean(draws[,z]),numeric(1))
      qvar <- function(names,p) vapply(names,function(z)
        unname(quantile(draws[,z],p)),numeric(1))
      pp <- tibble(
        aquifer_order=order,
        fold=fold,
        sample_index=held,
        Site_ID=d$Site_ID[held],
        Depth=d$Depth[held],
        observed_DO=d$DO[held],
        predicted_DO=as.numeric(predmean),
        lower95=qvar(pr,.025),
        upper95=qvar(pr,.975)
      ) %>%
        mutate(
          residual=predicted_DO-observed_DO,
          covered95=(observed_DO>=lower95&
                     observed_DO<=upper95)
        )
      cvpred[[length(cvpred)+1L]] <- pp
      cvdiag[[length(cvdiag)+1L]] <-
        select_diag(fit,order,fold,held)
      write_csv(bind_rows(cvpred),
                "outputs/final_DO_23_cv_predictions.csv")
      write_csv(bind_rows(cvdiag),
                "outputs/final_DO_23_cv_diagnostics.csv")
    }
  }
  preds <- bind_rows(cvpred)
  scores <- preds %>% group_by(aquifer_order) %>%
    summarise(
      n=n(),
      folds_completed=n_distinct(fold),
      DO_RMSE=sqrt(mean(residual^2)),
      DO_MAE=mean(abs(residual)),
      bias=mean(residual),
      coverage95=mean(covered95),
      .groups="drop"
    )
  write_csv(scores,"outputs/final_DO_23_cv_scores.csv")
  paired <- preds %>%
    select(
      sample_index,Site_ID,Depth,fold,
      observed_DO,aquifer_order,predicted_DO
    ) %>%
    pivot_wider(
      names_from=aquifer_order,
      values_from=predicted_DO,
      names_prefix="order_"
    ) %>%
    mutate(
      abs_error_zero=abs(observed_DO-order_0),
      abs_error_first=abs(observed_DO-order_1),
      first_minus_zero_abs_error=
        abs_error_first-abs_error_zero
    )
  write_csv(paired,"outputs/final_DO_23_cv_paired_errors.csv")
  cat("\nPAIRED GROUPED HELD-OUT DO SCORES:\n")
  print(scores,n=Inf)
  p <- preds %>%
    ggplot(aes(observed_DO,predicted_DO))+
    geom_abline(slope=1,intercept=0,
                linetype="dashed",colour="grey60")+
    geom_point(alpha=.8,size=2)+
    facet_wrap(~aquifer_order)+
    coord_equal()+
    labs(
      title="Final DO model: aquifer zero versus first order",
      subtitle="Grouped held-out wells; source fractions jointly fitted",
      x="Observed DO (mg/L)",
      y="Posterior predictive mean DO (mg/L)"
    )+theme_bw(base_size=12)
  print(p)
  ggsave("plots/final_DO_23_cv.png",
         p,width=10,height=5.5,dpi=300)
}
p <- main %>%
  filter(variable %in% c(
    "LSR_DO_at_water_table","LSR_DO_deficit",
    "k_ref","river_DO","survival_LSR",
    "mean_vadose_LSR_DO_loss","mean_aquifer_DO_loss"
  )) %>%
  ggplot(aes(median,factor(aquifer_order)))+
  geom_segment(aes(
    x=lower95,xend=upper95,
    y=factor(aquifer_order),yend=factor(aquifer_order)
  ))+
  geom_point(size=2.2)+
  facet_wrap(~variable,scales="free_x")+
  labs(
    title="Aquifer kinetics and LSR DO attenuation",
    subtitle="Joint posterior medians and central 95% credible intervals",
    x="Parameter",y="Aquifer kinetic order"
  )+theme_bw(base_size=11)
print(p)
ggsave("plots/final_DO_23_posteriors.png",
       p,width=12,height=7,dpi=300)

# Optional full-only boundary/age/LSR-prior sensitivity, same
# aquifer-order comparison on each modified physically explicit input.
if (cfg$run_prior_sensitivity) {
  cases <- list(
    list(name="rain10C",T=10,river=cfg$river_DO_boundary,
         sd=cfg$log_mrt_SD,ba=1,bb=1),
    list(name="rain20C",T=20,river=cfg$river_DO_boundary,
         sd=cfg$log_mrt_SD,ba=1,bb=1),
    list(name="river6p5",T=cfg$recharge_temperature_C,
         river=6.5,sd=cfg$log_mrt_SD,ba=1,bb=1),
    list(name="river8p5",T=cfg$recharge_temperature_C,
         river=8.5,sd=cfg$log_mrt_SD,ba=1,bb=1),
    list(name="wideMRT",T=cfg$recharge_temperature_C,
         river=cfg$river_DO_boundary,sd=.5,ba=1,bb=1),
    list(name="LSRstrongloss",T=cfg$recharge_temperature_C,
         river=cfg$river_DO_boundary,sd=cfg$log_mrt_SD,ba=2,bb=5),
    list(name="LSRweakloss",T=cfg$recharge_temperature_C,
         river=cfg$river_DO_boundary,sd=cfg$log_mrt_SD,ba=5,bb=2)
  )
  sensitivity <- list()
  for (caseid in seq_along(cases)) {
    cas <- cases[[caseid]]
    x <- stan_data
    x$Csat_rain <- do_sat_mgL(
      cas$T,cfg$atmospheric_pressure_ratio
    )
    x$river_DO_boundary <- cas$river
    x$log_mrt_SD <- cas$sd
    x$beta_survival_a <- cas$ba
    x$beta_survival_b <- cas$bb
    for (order in 0:1) {
      message("Full posterior sensitivity ",cas$name,
              " order ",order)
      fit <- sample_model(order,x=x,scenario=caseid)
      sensitivity[[length(sensitivity)+1L]] <-
        posterior95(fit,c(
          "river_DO","survival_LSR",
          "LSR_DO_at_water_table","LSR_DO_deficit",
          "k_ref","mean_vadose_LSR_DO_loss",
          "mean_aquifer_DO_loss"
        )) %>% mutate(
          scenario=cas$name,order=order,.before=1
        )
      check <- select_diag(fit,order)
      if (any(check$rhat>cfg$max_rhat,na.rm=TRUE))
        warning("Sensitivity ",cas$name,
                " has unresolved convergence.")
      write_csv(
        bind_rows(sensitivity),
        "outputs/final_DO_23_sensitivity_posteriors.csv"
      )
    }
  }
}

cat("\nPOSTERIOR PARAMETERS (full data):\n")
print(main %>% filter(
  variable %in% c(
    "river_DO","survival_LSR",
    "LSR_DO_at_water_table","LSR_DO_deficit","k_ref",
    "aquifer_rate_actual",
    "mean_vadose_LSR_DO_loss","mean_aquifer_DO_loss",
    "sampled_well_vadose_loss_fraction"
  )
),n=Inf)
cat(
  "\nPHYSICAL INTERPRETATION AND LIMITATIONS:\n",
  "* Rain DO SATURATION concentration fixed by assumed recharge T, P.\n",
  "* LSR LOSS estimates are surface-to-water-table DO deficits,\n",
  "  not independently identifiable LSR reaction RATE constants.\n",
  "* River aquifer-entry DO is FIXED at 7.5 mg/L.\n",
  "  Site 2801's model-only DO=7.5 is excluded from fitting.\n",
  "* ZERO- and FIRST-order AQUFER kinetics are alternative models.\n",
  "* Zero-order k0=k_ref*C_ref [mg/L/yr], DO clipped to zero.\n",
  "* First-order k1=k_ref [1/yr], exponential DO loss.\n",
  "* LSR rates are conditional zero/first order for assumed vadose\n",
  "  travel times: no unique LSR rate from these groundwater data.\n",
  "* Chemical source fractions are latent and correlated Na/Cl is\n",
  "  handled jointly; excess air excluded by default.\n",
  "* The aquifer travel time is represented by uncertain MRT, but\n",
  "  assumes the reported MRT primarily represents aquifer residence.\n",
  "* Different source-specific MRTs and recharge-age distributions\n",
  "  could alter these apparent reaction coefficients substantially.\n",
  "* Arithmetic well-mean losses are not a catchment flow budget.\n",
  "* No DO/NO3 or denitrification data are used to double-count loss.\n",
  "* Model order chosen by held-out DO prediction AND convergence.\n",
  "* Run boundary, MRT and LSR prior sensitivities before publication.\n",
  sep=""
)
