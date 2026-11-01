# 01_DO_dual_pathway_kinetics.R — Bayesian fitting and diagnostics
#
# The ONLY active MCMC runner. Source 00_setup.R first to prepare a fresh,
# coherent configuration and input dataset, then fit competing zero- and
# first-order aquifer oxygen-reduction models. No source or Stan model changes.
#
# Run from the repository root:
#   source("scripts/00_setup.R")   # optional: inspect setup without sampling
#   source("scripts/01_DO_dual_pathway_kinetics.R") # fits both models
#
# Script 01 ALWAYS re-sources 00_setup.R to prevent stale configurations.
# To configure run_CV or sample settings, edit cfg in 00_setup.R.
# The sole private input is data/Uwe_HP_data_prepared.csv.

if (!file.exists("scripts/00_setup.R"))
  stop("Set working directory to repository root (missing scripts/00_setup.R).")
source("scripts/00_setup.R", local=TRUE)

# Setup is side-effect free apart from reading data. Write model audit and
# create output directories ONLY when the fit runner is called.
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
dir.create(plotdir, recursive=TRUE, showWarnings=FALSE)
write_csv(
  tibble(
    quantity=c(
      "rain_DO_surface","river_DO_entry",
      "LSR_DO_water_table","aquifer_reaction",
      "n_wells","n_fit_DO","n_O18","n_Cl"
    ),
    assumed_value=c(
      Csat,cfg$river_DO_boundary,NA_real_,NA_real_,
      nrow(d),length(eligible),
      sum(stan_data$has_O18),sum(stan_data$has_Cl)
    ),
    description=c(
      paste0(
        "Fixed oxygen saturation at ",cfg$recharge_temperature_C,
        " C; atmospheric pressure ratio ",cfg$atmospheric_pressure_ratio
      ),
      paste0("Fixed ",cfg$river_DO_boundary,
             " mg/L; proxy well ",cfg$river_proxy_site,
             " DO not fitted"),
      "Posterior: saturated rainfall DO times survival_LSR",
      "Alternative aquifer zero- or first-order kinetics",
      paste0("After excluding ",cfg$excluded_site,
             " and filtering DO/MRT"),
      paste0("Proxy ",cfg$river_proxy_site," DO masked"),
      "O18 present in original data",
      "Cl present in original data"
    )
  ),
  output("boundary_conditions.csv")
)
# Posterior summaries must use numeric vectors: posterior draws_matrix
# columns otherwise carry attributes that can corrupt tibble operations.
posterior95 <- function(fit, vars) {
  sm <- fit$summary(variables=vars)
  x <- posterior::as_draws_matrix(fit$draws(variables=vars))
  sm %>%
    select(variable,mean,median,sd,rhat,ess_bulk,ess_tail) %>%
    mutate(
      lower95=vapply(variable,function(v)
        unname(quantile(as.numeric(x[,v]),0.025)),numeric(1)),
      upper95=vapply(variable,function(v)
        unname(quantile(as.numeric(x[,v]),0.975)),numeric(1))
    )
}

global_vars <- c(
  "river_DO","survival_LSR","LSR_DO_at_water_table",
  "LSR_DO_deficit","k_ref","aquifer_rate_actual","sigma_DO",
  "mean_f_LSR","mean_vadose_LSR_DO_loss",
  "mean_aquifer_DO_loss","sampled_well_vadose_loss_fraction",
  "O18_R","O18_LSR","Cl_R","Cl_LSR",
  "sigma_O18","tau_LSR_O18","sigma_Cl"
)

# Omit the deliberately constant proxy MRT and the fixed river boundary
# from convergence calculations; all estimated latent and mixing
# parameters are included (plus relevant transformed quantities).
check_fit <- function(fit, order, fold=0L) {
  sm <- fit$summary(variables=c(
    global_vars,"f_LSR","true_MRT","z_mrt"
  )) %>%
    filter(
      variable!="river_DO",
      variable!=paste0("true_MRT[",proxy_index,"]")
    )
  hm <- fit$diagnostic_summary()
  sampler <- tibble(
    aquifer_order=order,
    fold=fold,
    chain=seq_len(cfg$chains),
    n_divergences=as.integer(hm$num_divergent),
    n_max_treedepth=as.integer(hm$num_max_treedepth),
    ebfmi=as.numeric(hm$ebfmi)
  )
  bad_number <- any(
    !is.finite(sm$rhat) | !is.finite(sm$ess_bulk)
  ) || any(!is.finite(sampler$ebfmi))
  worst_rhat <- max(sm$rhat,na.rm=TRUE)
  min_ess <- min(sm$ess_bulk,na.rm=TRUE)
  n_div <- sum(sampler$n_divergences)
  n_tree <- sum(sampler$n_max_treedepth)
  min_ebfmi <- min(sampler$ebfmi)
  passed <- !bad_number &&
    worst_rhat<=cfg$max_rhat &&
    min_ess>=cfg$min_bulk_ESS &&
    n_div==0 && n_tree==0 && min_ebfmi>=0.3
  summary <- tibble(
    aquifer_order=order,fold=fold,
    max_rhat=worst_rhat,min_bulk_ess=min_ess,
    n_divergences=n_div,n_max_treedepth=n_tree,
    min_ebfmi=min_ebfmi,all_diagnostics_pass=passed
  )
  list(summary=summary,sampler=sampler,
       parameter_diagnostics=sm)
}

mod <- cmdstanr::cmdstan_model(cfg$stan,quiet=TRUE)
sample_model <- function(order, x=stan_data, fold=0L) {
  x$aq_order <- as.integer(order)
  mod$sample(
    data=x,
    seed=cfg$seed+101L*order+10000L*fold,
    chains=cfg$chains,
    parallel_chains=cfg$parallel_chains,
    iter_warmup=cfg$warmup,
    iter_sampling=cfg$sampling,
    adapt_delta=cfg$adapt_delta,
    max_treedepth=cfg$max_treedepth,
    refresh=250
  )
}

post <- list();diag <- list();sampler_tables <- list()
well_tables <- list();idents <- list();draw_tables <- list()
for (order in 0:1) {
  message("Full fit: aquifer order ",order)
  fit <- sample_model(order)
  # Save FIRST, before generating any reports.
  fit$save_object(file=output(paste0("fit_order",order,".rds")))
  chk <- check_fit(fit,order)
  diag[[order+1L]] <- chk$summary
  sampler_tables[[order+1L]] <- chk$sampler
  write_csv(bind_rows(diag),output("fit_status.csv"))
  write_csv(bind_rows(sampler_tables),output("sampler_chains.csv"))
  write_csv(chk$parameter_diagnostics,
            output(paste0("parameter_diagnostics_order",order,".csv")))

  post[[order+1L]] <- posterior95(fit,global_vars) %>%
    mutate(
      aquifer_order=order,
      rate_units=if_else(
        variable=="aquifer_rate_actual",
        if(order==0) "mg/L/year" else "year^-1",
        NA_character_
      ),
      .before=1
    )
  write_csv(bind_rows(post),output("global_posteriors.csv"))

  wells <- fit$summary(variables=c(
    "f_LSR","true_MRT","mean_DO",
    "vadose_LSR_DO_loss","aquifer_DO_loss"
  ))
  well_tables[[order+1L]] <- wells %>%
    mutate(
      sample_index=as.integer(
        stringr::str_match(variable,"\\[(\\d+)\\]")[,2]
      ),
      aquifer_order=order,.before=1
    ) %>%
    inner_join(
      d %>% select(sample_index,Site_ID,Depth,
                   DO_original,MRT),
      by="sample_index"
    ) %>%
    mutate(DO_used_in_likelihood=
             as.integer(sample_index!=proxy_index))
  write_csv(bind_rows(well_tables),output("well_posteriors.csv"))

  z <- posterior::as_draws_matrix(fit$draws(variables=c(
    "LSR_DO_deficit","aquifer_rate_actual",
    "mean_vadose_LSR_DO_loss","mean_aquifer_DO_loss"
  )))
  loss <- as.numeric(z[,"LSR_DO_deficit"])
  rate <- as.numeric(z[,"aquifer_rate_actual"])
  vadose <- as.numeric(z[,"mean_vadose_LSR_DO_loss"])
  aquifer <- as.numeric(z[,"mean_aquifer_DO_loss"])
  a <- cfg$survival_prior_a
  b <- cfg$survival_prior_b
  prior_sd <- Csat*sqrt(a*b/((a+b)^2*(a+b+1)))
  idents[[order+1L]] <- tibble(
    aquifer_order=order,
    pearson_deficit_rate=cor(loss,rate),
    spearman_deficit_rate=cor(loss,rate,method="spearman"),
    correlation_vadose_aquifer_loss=cor(vadose,aquifer),
    LSR_deficit_posterior_sd=sd(loss),
    LSR_deficit_prior_sd=prior_sd,
    posterior_to_prior_SD=sd(loss)/prior_sd,
    P_LSR_deficit_gt_5=mean(loss>5),
    P_vadose_loss_gt_aquifer_loss=mean(vadose>aquifer)
  )
  set.seed(cfg$seed+order)
  keep <- sample(seq_along(loss),min(length(loss),3500L))
  draw_tables[[order+1L]] <- tibble(
    aquifer_order=order,
    LSR_DO_deficit_mgL=loss[keep],
    aquifer_rate=rate[keep],
    nominal_vadose_loss_mgL=vadose[keep],
    nominal_aquifer_loss_mgL=aquifer[keep],
    rate_units=if(order==0) "mg/L/year" else "year^-1"
  )
  write_csv(bind_rows(idents),output("identifiability.csv"))
  write_csv(bind_rows(draw_tables),output("joint_draws.csv"))
}

full_status <- bind_rows(diag)
cat("\nFULL-FIT MCMC STATUS\n")
print(full_status,n=Inf,width=Inf)
if (!all(full_status$all_diagnostics_pass))
  warning("At least one full model failed MCMC checks. ",
          "Do not interpret its parameter intervals or compare orders yet.")

cat("\nCORE PARAMETER POSTERIORS\n")
print(bind_rows(post) %>%
        filter(variable %in% c(
          "LSR_DO_deficit","LSR_DO_at_water_table",
          "aquifer_rate_actual","mean_f_LSR","sigma_DO"
        )) %>%
        select(aquifer_order,variable,median,lower95,upper95,
               rhat,ess_bulk,rate_units),
      n=Inf,width=Inf)

p <- bind_rows(draw_tables) %>%
  ggplot(aes(LSR_DO_deficit_mgL,aquifer_rate))+
  geom_point(alpha=0.15,size=0.55)+
  facet_wrap(~aquifer_order,scales="free_y",
             labeller=as_labeller(c(
               "0"="Zero order (mg/L/year)",
               "1"="First order (year^-1)"
             )))+
  labs(
    title="Two-endmember DO model: LSR deficit versus aquifer rate",
    subtitle="Delta18O + chloride; posterior correlation is not causal",
    x="LSR DO deficit to water table (mg/L)",
    y="Aquifer reaction rate"
  )+theme_bw(base_size=12)
ggsave(file.path(plotdir,"identifiability.png"),
       p,width=10,height=5,dpi=300)

# OPTIONAL grouped held-out DO prediction; retains chemical measurements
# for all wells but removes DO for the held-out group. Same folds/seed
# in zero and first order. Site 2801 is never a scored DO observation.
if (cfg$run_CV) {
  if (!all(full_status$all_diagnostics_pass))
    stop("CV disabled: full fits must first pass MCMC diagnostics.")
  stopifnot(cfg$folds>=2L, is.null(cfg$run_only_fold) ||
              all(cfg$run_only_fold %in% seq_len(cfg$folds)))

  set.seed(cfg$seed)
  groups <- split(eligible,d$cv_group[eligible])
  g_order <- sample(seq_along(groups))
  g_order <- g_order[base::order(-lengths(groups)[g_order])]
  folds <- vector("list",cfg$folds)
  for (i in g_order) {
    j <- which.min(lengths(folds))
    folds[[j]] <- c(folds[[j]],groups[[i]])
  }
  stopifnot(
    all(lengths(folds)>0),
    identical(sort(unlist(folds)),eligible)
  )
  foldmap <- purrr::map_dfr(seq_along(folds),function(i)
    tibble(
      fold=i, sample_index=folds[[i]],
      Site_ID=d$Site_ID[folds[[i]]],
      cv_group=d$cv_group[folds[[i]]]
    )
  )
  write_csv(foldmap,output("cv_folds.csv"))

  cv_preds <- list();cv_checks <- list()
  chosen <- if (is.null(cfg$run_only_fold))
    seq_along(folds) else cfg$run_only_fold
  for (fold in chosen) {
    held <- folds[[fold]]
    for (order in 0:1) {
      x <- stan_data
      x$has_DO[held] <- 0L
      x$DO_obs[held] <- 0
      message("Held-out fold ",fold," order ",order)
      fit <- sample_model(order,x,fold=fold)
      fit$save_object(file=output(paste0(
        "cv_fit_fold",fold,"_order",order,".rds"
      )))
      chk <- check_fit(fit,order,fold=fold)
      cv_checks[[length(cv_checks)+1L]] <- chk$summary
      write_csv(bind_rows(cv_checks),output("cv_fit_status.csv"))
      if (!chk$summary$all_diagnostics_pass)
        warning("Non-converged CV fit: fold ",fold,
                ", order ",order)

      mvars <- paste0("mean_DO[",held,"]")
      rvars <- paste0("replicated_DO[",held,"]")
      draws <- posterior::as_draws_matrix(
        fit$draws(variables=c(mvars,rvars))
      )
      predmean <- vapply(mvars,function(v)
        mean(as.numeric(draws[,v])),numeric(1))
      qt <- function(vars,p) vapply(vars,function(v)
        unname(quantile(as.numeric(draws[,v]),p)),numeric(1))
      cv_preds[[length(cv_preds)+1L]] <- tibble(
        aquifer_order=order,fold=fold,
        sample_index=held,
        Site_ID=d$Site_ID[held],
        cv_group=d$cv_group[held],
        observed_DO=d$DO_original[held],
        predicted_DO=predmean,
        predictive_lower95=qt(rvars,0.025),
        predictive_upper95=qt(rvars,0.975)
      ) %>%
        mutate(
          residual=predicted_DO-observed_DO,
          covered95=(observed_DO>=predictive_lower95 &
                       observed_DO<=predictive_upper95)
        )
      write_csv(bind_rows(cv_preds),output("cv_predictions.csv"))
    }
  }

  predictions <- bind_rows(cv_preds)
  scores <- predictions %>%
    group_by(aquifer_order) %>%
    summarise(
      n=n(),
      folds_completed=n_distinct(fold),
      DO_RMSE=sqrt(mean(residual^2)),
      DO_MAE=mean(abs(residual)),
      bias=mean(residual),
      predictive_coverage95=mean(covered95),
      .groups="drop"
    )
  write_csv(scores,output("cv_scores.csv"))
  cat("\nGROUPED HELD-OUT DO SCORES\n")
  print(scores,n=Inf,width=Inf)

  if (all(bind_rows(cv_checks)$all_diagnostics_pass) &&
      is.null(cfg$run_only_fold)) {
    paired <- predictions %>%
      select(
        sample_index,Site_ID,cv_group,fold,
        observed_DO,aquifer_order,predicted_DO
      ) %>%
      pivot_wider(
        names_from=aquifer_order,
        values_from=predicted_DO,
        names_prefix="order_"
      )
    write_csv(paired,output("cv_paired_predictions.csv"))
    grouped <- split(paired,paired$cv_group)
    set.seed(cfg$seed+202L)
    boot <- purrr::map_dfr(seq_len(2000L),function(i) {
      pick <- sample(seq_along(grouped),length(grouped),
                     replace=TRUE)
      x <- bind_rows(grouped[pick])
      tibble(
        draw=i,
        delta_RMSE=
          sqrt(mean((x$order_1-x$observed_DO)^2))-
          sqrt(mean((x$order_0-x$observed_DO)^2)),
        delta_MAE=
          mean(abs(x$order_1-x$observed_DO))-
          mean(abs(x$order_0-x$observed_DO))
      )
    })
    write_csv(boot,output("cv_group_bootstrap.csv"))
    boot_summary <- tibble(
      measure=c("first_minus_zero_RMSE","first_minus_zero_MAE"),
      median=c(median(boot$delta_RMSE),median(boot$delta_MAE)),
      lower95=c(unname(quantile(boot$delta_RMSE,0.025)),
                unname(quantile(boot$delta_MAE,0.025))),
      upper95=c(unname(quantile(boot$delta_RMSE,0.975)),
                unname(quantile(boot$delta_MAE,0.975)))
    )
    write_csv(boot_summary,output("cv_group_bootstrap_summary.csv"))
    cat("\nGROUP-BOOTSTRAP PREDICTIVE DIFFERENCES\n")
    print(boot_summary,n=Inf,width=Inf)
  } else {
    warning("Do not rank kinetic orders from incomplete or ",
            "non-converged CV fits.")
  }
}

cat(
  "\nINTERPRETATION:\n",
  "- Two endmembers: river recharge and LSR, constrained by O18 and Cl.\n",
  "- River DO=7.5 mg/L fixed; rainfall DO initially saturated.\n",
  "- LSR loss is a DO DEFICIT before the water table, not a known rate.\n",
  "- Aquifer zero- and first-order fits have DIFFERENT rate units.\n",
  "- Mixing and reactions are conditional on assumed endmembers/MRT.\n",
  "- Sodium and excess air are excluded; DO is the only reactive fit.\n",
  "- Site 2801 DO is masked; 16758 excluded; original CSV untouched.\n",
  "- Well-average losses are not catchment-scale mass balances.\n",
  "- Convergence is necessary, not proof of chemical identifiability.\n",
  "- Use grouped DO CV and boundary/MRT sensitivity before publication.\n",
  sep=""
)
