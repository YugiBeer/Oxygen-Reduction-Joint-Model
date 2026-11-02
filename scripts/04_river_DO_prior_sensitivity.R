# 04_river_DO_prior_sensitivity.R
#
# Controlled comparison of a FIXED river aquifer-entry DO boundary against
# informative uncertain boundaries. Uses the two-endmember O18 + chloride
# mixing model and competing zero-/first-order aquifer DO kinetics.
#
# DEFAULT: two scenarios * two kinetic orders = four independent Stan fits.
#   A fixed:      river DO = 7.5 mg/L — preferred regional primary model
#   C prior_sd1p5:river DO ~ Normal(7.5,1.5^2) truncated [0,14.65]
# Optional B prior_sd1 can be restored in cfg04$scenarios if needed.
#
# Site 2801: measured groundwater DO=10.2 mg/L, MRT assumed zero.
# PRIMARY excludes its DO from fitting, but keeps O18 and chloride.
# Site 2130: shallow well next to water race; PRIMARY excludes its DO
# from fitting, but retains O18 and chloride in the mixing model.
# Neither site is removed from the groundwater input CSV, and no values
# are replaced in 00_setup.R. The observation mask alone changes.
#
# Results from this 39-DO-site regional fit go into a distinct output
# folder ("2801_masked_no2130"). Earlier 40- and 41-well results stay
# untouched, even when the experiment is rerun.
#
# Keeps primary 01 outputs, primary Stan file and 02 plots untouched.
# Data file remains local and unchanged. Samples saved before diagnostics.
# Run from repo root:
#   source("scripts/04_river_DO_prior_sensitivity.R")

if (!file.exists("scripts/00_setup.R"))
  stop("Set working directory to repository root.")
source("scripts/00_setup.R",local=TRUE)

cfg04 <- list(
  stan="stan/04_river_DO_prior_sensitivity.stan",
  include_proxy_DO=FALSE, # Site 2801: exclude DO; retain O18 and Cl
  exclude_DO_sites="2130", # Exclude DO; retain O18 and Cl
  scenarios=tibble::tibble(
    scenario=c("A_fixed_7p5","C_prior_sd1p5"),
    infer_river_DO=c(0L,1L),
    river_DO_prior_sd=c(0,1.5)
  ),
  chains=cfg$chains,
  parallel_chains=cfg$parallel_chains,
  iter_warmup=cfg$warmup,
  iter_sampling=cfg$sampling,
  adapt_delta=cfg$adapt_delta,
  max_treedepth=cfg$max_treedepth,
  seed=cfg$seed+4000L,
  max_rhat=1.05,
  min_bulk_ESS=400L
)

if (!file.exists(cfg04$stan))
  stop("Missing Stan model: ",cfg04$stan)
if (any(cfg04$scenarios$infer_river_DO==1L &
        cfg04$scenarios$river_DO_prior_sd<=0))
  stop("Inferred river DO requires a positive prior SD.")
if (any(cfg04$scenarios$infer_river_DO==0L &
        cfg04$scenarios$river_DO_prior_sd!=0))
  stop("Fixed scenario requires prior SD = 0.")
if (cfg$river_DO_boundary>=14.65)
  stop("River DO prior centre exceeds upper Stan constraint.")
if (!identical(d$DO[proxy_index],d$DO_original[proxy_index]))
  stop("The original DO for proxy 2801 was overwritten in setup. ",
       "Update and rerun scripts/00_setup.R.")
if (!is.finite(d$DO_original[proxy_index]))
  stop("2801 measured DO is missing; cannot run matched inclusion test.")
if (anyDuplicated(d$Site_ID[d$Site_ID==cfg$river_proxy_site]))
  stop("Duplicated 2801 proxy.")
# Explicit DO-only exclusion. The tracer records are not dropped.
exclude_DO_sites <- unique(as.character(cfg04$exclude_DO_sites))
if (anyNA(exclude_DO_sites) ||
    any(!grepl("^[A-Za-z0-9_]+$",exclude_DO_sites)))
  stop("exclude_DO_sites must contain valid Site_ID strings.")
if (any(!exclude_DO_sites %in% d$Site_ID))
  stop("Unknown excluded DO site(s): ",
       paste(setdiff(exclude_DO_sites,d$Site_ID),collapse=", "))
if (cfg$river_proxy_site %in% exclude_DO_sites)
  stop("2801 already has a separate DO inclusion control.")
excluded_ix <- which(d$Site_ID %in% exclude_DO_sites)
if (!("2130" %in% d$Site_ID))
  stop("Site 2130 missing from current prepared data.")

suite <- if (cfg04$include_proxy_DO) "2801_included" else "2801_masked"
if (length(exclude_DO_sites))
  suite <- paste0(suite,"_no",paste(sort(exclude_DO_sites),collapse="_"))
out04 <- file.path("outputs","04_river_DO_prior",suite)
dir.create(out04,recursive=TRUE,showWarnings=FALSE)
file04 <- function(name) file.path(out04,name)

# Start with the exact same chemistry, ages, and priors as the 01
# reference model. All comparisons use the SAME fixed DO mask.
xbase <- stan_data
xbase$has_DO[proxy_index] <- as.integer(cfg04$include_proxy_DO)
xbase$DO_obs <- as.numeric(d$DO_original)
# Make sure excluded DO cannot enter the Stan likelihood even indirectly.
xbase$has_DO[excluded_ix] <- 0L
xbase$DO_obs[excluded_ix] <- 0
if (cfg04$include_proxy_DO &&
    (xbase$has_DO[proxy_index]!=1L ||
     xbase$DO_obs[proxy_index]!=d$DO_original[proxy_index]))
  stop("Failed to include original 2801 DO.")
if (!cfg04$include_proxy_DO && xbase$has_DO[proxy_index]!=0L)
  stop("Failed to mask 2801 DO.")
stopifnot(
  xbase$N==nrow(d),
  sum(xbase$has_DO)==nrow(d)-as.integer(!cfg04$include_proxy_DO)-
    length(excluded_ix),
  all(xbase$has_DO[excluded_ix]==0L),
  all(xbase$DO_obs[excluded_ix]==0),
  identical(xbase$has_O18,stan_data$has_O18),
  identical(xbase$has_Cl,stan_data$has_Cl)
)

audit <- cfg04$scenarios %>%
  mutate(
    observed_2801_DO=d$DO_original[proxy_index],
    observed_2130_DO=d$DO_original[match("2130",d$Site_ID)],
    river_entry_prior_mean=cfg$river_DO_boundary,
    proxy_DO_used=cfg04$include_proxy_DO,
    excluded_DO_sites=paste(exclude_DO_sites,collapse=";"),
    excluded_DO_n=length(excluded_ix),
    n_independent_DO=sum(xbase$has_DO),
    rainfall_surface_DO=Csat,
    log_mrt_SD=cfg$log_mrt_SD,
    endmember_prior_scale=cfg$endmember_prior_scale
  )
readr::write_csv(audit,file04("scenario_design.csv"))
message("River DO prior sensitivity: ",suite,
        " | 2801 actual DO=",d$DO_original[proxy_index],
        " mg/L | DO likelihood wells=",sum(xbase$has_DO))
if (length(exclude_DO_sites)) {
  message("Excluded from DO likelihood: ",
          paste(exclude_DO_sites,collapse=", "),
          " (O18 and Cl preserved).")
}
message("Fit ",nrow(cfg04$scenarios),
        " scenarios * 2 aquifer orders, identical observation masks.")

# Check all active population parameters and per-well latents.
active <- c(
  "O18_R","O18_LSR","Cl_R","Cl_LSR",
  "f_LSR","z_mrt","sigma_O18","sigma_Cl","tau_LSR_O18",
  "survival_LSR","k_ref","sigma_DO"
)
quantiles95 <- function(fit,vars) {
  sm <- fit$summary(variables=vars)
  mat <- posterior::as_draws_matrix(fit$draws(variables=vars))
  sm %>%
    dplyr::select(variable,mean,median,sd,rhat,ess_bulk,ess_tail) %>%
    dplyr::mutate(
      lower95=vapply(variable,function(v)
        unname(stats::quantile(as.numeric(mat[,v]),.025)),numeric(1)),
      upper95=vapply(variable,function(v)
        unname(stats::quantile(as.numeric(mat[,v]),.975)),numeric(1))
    )
}
sampler_check <- function(fit,sc,order,infer_river) {
  parnames <- c(active,if (infer_river==1L) "river_DO_free")
  sm <- fit$summary(variables=parnames)
  hm <- fit$diagnostic_summary()
  n_div <- sum(hm$num_divergent,na.rm=TRUE)
  n_tree <- sum(hm$num_max_treedepth,na.rm=TRUE)
  min_ebfmi <- min(hm$ebfmi,na.rm=TRUE)
  max_rhat <- max(sm$rhat,na.rm=TRUE)
  min_ess <- min(sm$ess_bulk,na.rm=TRUE)
  bad_numeric <- any(!is.finite(sm$rhat)) ||
    any(!is.finite(sm$ess_bulk)) ||
    any(!is.finite(hm$ebfmi))
  passed <- !bad_numeric && n_div==0 && n_tree==0 &&
    is.finite(min_ebfmi) && min_ebfmi>=.3 &&
    max_rhat<=cfg04$max_rhat && min_ess>=cfg04$min_bulk_ESS
  list(
    status=tibble::tibble(
      scenario=sc,aquifer_order=order,
      max_rhat=max_rhat,min_bulk_ess=min_ess,
      divergences=n_div,max_treedepth_hits=n_tree,
      min_ebfmi=min_ebfmi,all_diagnostics_pass=passed
    ),
    chains=tibble::tibble(
      scenario=sc,aquifer_order=order,
      chain=seq_along(hm$num_divergent),
      num_divergent=as.integer(hm$num_divergent),
      num_max_treedepth=as.integer(hm$num_max_treedepth),
      ebfmi=as.numeric(hm$ebfmi)
    )
  )
}

parameters_of_interest <- c(
  "river_DO","LSR_DO_at_water_table","LSR_DO_deficit",
  "aquifer_rate_actual","k_ref","sigma_DO",
  "survival_LSR","mean_f_LSR",
  "mean_vadose_LSR_DO_loss","mean_aquifer_DO_loss",
  "O18_R","O18_LSR","Cl_R","Cl_LSR"
)
post <- list();status <- list();chain_rows <- list()
wells <- list();correlations <- list();predictions <- list()
mod04 <- cmdstanr::cmdstan_model(cfg04$stan,quiet=TRUE)
for (sc_i in seq_len(nrow(cfg04$scenarios))) {
  sc <- cfg04$scenarios[sc_i,]
  for (order in 0:1) {
    x <- xbase
    x$aq_order <- as.integer(order)
    x$infer_river_DO <- as.integer(sc$infer_river_DO)
    x$river_DO_prior_sd <- as.numeric(sc$river_DO_prior_sd)
    excluded_label <- if (length(exclude_DO_sites)) {
      paste(exclude_DO_sites,collapse=", ")
    } else {
      "none"
    }
    message("\nFIT ",sc$scenario," | aquifer order=",order,
            " | 2801 DO ",if (cfg04$include_proxy_DO) "included" else "masked",
            " | DO excluded: ",excluded_label)
    fit <- mod04$sample(
      data=x,
      seed=as.integer(cfg04$seed+201L*sc_i+10009L*order),
      chains=cfg04$chains,
      parallel_chains=cfg04$parallel_chains,
      iter_warmup=cfg04$iter_warmup,
      iter_sampling=cfg04$iter_sampling,
      adapt_delta=cfg04$adapt_delta,
      max_treedepth=cfg04$max_treedepth,
      refresh=250
    )
    fit$save_object(file=file04(paste0(
      "fit_",sc$scenario,"_order",order,".rds"
    )))
    chk <- sampler_check(fit,sc$scenario,order,sc$infer_river_DO)
    status[[length(status)+1L]] <- chk$status
    chain_rows[[length(chain_rows)+1L]] <- chk$chains
    readr::write_csv(dplyr::bind_rows(status),file04("fit_status.csv"))
    readr::write_csv(dplyr::bind_rows(chain_rows),
                     file04("sampler_chains.csv"))
    if (!chk$status$all_diagnostics_pass)
      warning("Sampler diagnostics FAIL for ",sc$scenario,
              " order ",order,". Do NOT interpret those posteriors.")

    post[[length(post)+1L]] <- quantiles95(
      fit,parameters_of_interest
    ) %>%
      dplyr::mutate(
        scenario=sc$scenario,
        aquifer_order=order,
        river_DO_prior_sd=sc$river_DO_prior_sd,
        rate_units=ifelse(
          variable=="aquifer_rate_actual",
          if (order==0) "mg/L/year" else "year^-1",NA_character_
        ),
        .before=1
      )
    readr::write_csv(dplyr::bind_rows(post),
                     file04("global_posteriors.csv"))

    # The held-out target sites are predicted for diagnostic purposes.
    # Neither site's DO enters the primary analysis likelihood.
    for (siteid in c(cfg$river_proxy_site,"2130")) {
      i <- match(siteid,d$Site_ID)
      v <- paste0(
        c("mean_DO","replicated_DO","f_LSR","true_MRT"),
        "[",i,"]"
      )
      z <- posterior::as_draws_matrix(fit$draws(variables=v))
      vals <- lapply(v,function(n) as.numeric(z[,n]))
      names(vals) <- c("mu","rep","f","age")
      predictions[[length(predictions)+1L]] <- tibble::tibble(
        scenario=sc$scenario,aquifer_order=order,
        site_id=siteid,DO_was_fitted=as.logical(x$has_DO[i]),
        DO_observed=d$DO_original[i],
        reported_MRT=d$MRT[i],
        fitted_DO_median=median(vals$mu),
        fitted_DO_lower95=unname(quantile(vals$mu,.025)),
        fitted_DO_upper95=unname(quantile(vals$mu,.975)),
        replicated_DO_median=median(vals$rep),
        replicated_DO_lower95=unname(quantile(vals$rep,.025)),
        replicated_DO_upper95=unname(quantile(vals$rep,.975)),
        P_replicated_DO_at_least_observed=
          mean(vals$rep>=d$DO_original[i]),
        f_LSR_median=median(vals$f),
        true_MRT_median=median(vals$age)
      )
    }
    readr::write_csv(dplyr::bind_rows(predictions),
                     file04("high_DO_site_predictions.csv"))

    lat <- posterior::as_draws_matrix(fit$draws(variables=c(
      "river_DO","LSR_DO_deficit","aquifer_rate_actual",
      "mean_f_LSR"
    )))
    correlations[[length(correlations)+1L]] <- tibble::tibble(
      scenario=sc$scenario,aquifer_order=order,
      corr_riverDO_aquiferRate=if(sc$infer_river_DO==1L) {
        cor(as.numeric(lat[,"river_DO"]),
            as.numeric(lat[,"aquifer_rate_actual"]))
      } else {
        NA_real_
      },
      corr_riverDO_LSRdeficit=if(sc$infer_river_DO==1L) {
        cor(as.numeric(lat[,"river_DO"]),
            as.numeric(lat[,"LSR_DO_deficit"]))
      } else {
        NA_real_
      },
      corr_LSRdeficit_aquiferRate=cor(
        as.numeric(lat[,"LSR_DO_deficit"]),
        as.numeric(lat[,"aquifer_rate_actual"])
      ),
      corr_riverDO_LSRfraction=if(sc$infer_river_DO==1L) {
        cor(as.numeric(lat[,"river_DO"]),
            as.numeric(lat[,"mean_f_LSR"]))
      } else {
        NA_real_
      }
    )
    readr::write_csv(dplyr::bind_rows(correlations),
                     file04("posterior_correlations.csv"))
  }
}
full <- dplyr::bind_rows(status)
long <- dplyr::bind_rows(post)
summary <- long %>%
  dplyr::filter(variable %in% c(
    "river_DO","LSR_DO_deficit","LSR_DO_at_water_table",
    "aquifer_rate_actual","mean_f_LSR","sigma_DO"
  )) %>%
  dplyr::select(
    scenario,aquifer_order,variable,
    median,lower95,upper95,rate_units
  )
readr::write_csv(summary,file04("comparison_summary.csv"))
message("\nSAMPLER QUALITY (all fitted models):")
print(full,n=Inf,width=Inf)
message("\nRIVER DO PRIOR SENSITIVITY: MEDIANS AND 95% INTERVALS")
print(summary,n=Inf,width=Inf)
if (!all(full$all_diagnostics_pass))
  warning("At least one sampler failed diagnostics. Do not interpret ",
          "its posterior or make kinetic-model comparisons.")
message("\nSAVED OUTPUTS: ",normalizePath(out04))
message("These are IN-SAMPLE scenario comparisons. No held-out ",
        "validation has been performed.")
message("No modification to scripts/01, scripts/02 or the primary Stan model.")
