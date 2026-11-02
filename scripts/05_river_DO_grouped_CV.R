# 05_river_DO_grouped_CV.R — rigorous held-out DO comparison
#
# Compares the BROAD river-DO prior Normal(7.5, 1.5^2), truncated
# to [0,14.65] mg/L, against the original FIXED 7.5-mg/L boundary.
# Both models use identical source chemistry, ages, DO masks, folds and
# priors except the river DO prior, and are fitted separately under
# zero- and first-order aquifer kinetics.
#
# PRIMARY ANALYSIS: site 2801 DO *masked* from every fit AND never scored.
# Its O18 and Cl observations remain in ALL fits. This avoids turning a
# 10.2-mg/L river-proximal groundwater observation with assumed MRT=0
# into a strong implicit measurement of the global river boundary.
#
# OPTIONAL SENSITIVITY: set include_2801_DO=TRUE to include original
# measured 2801 DO in each TRAINING fit, but never in the held-out
# validation scores. All 40 other wells are held out in matched folds.
# The output directories are separate for the two experiments.
#
# Four model variants (A-fixed order 0/1, C-broad order 0/1) * five
# matched group folds = 20 cross-validation fits for each experiment.
# Saves fit objects immediately. Cache + design checksum enable a
# stopped run to resume without repeating completed fits.
#
# No refitting of full-data models 01 or 04. Zero DO observations from
# held-out wells enter their Stan likelihood; conservative tracers and
# reported ages DO remain available to predict each held-out well.
#
# Run from repository root:
#   source("scripts/05_river_DO_grouped_CV.R")
# After CV:
#   source("scripts/02_DO_plotting.R")  # Turbo-colour results plots

if (!file.exists("scripts/00_setup.R"))
  stop("Set working directory to repository root.")
source("scripts/00_setup.R",local=TRUE)

cfg05 <- list(
  stan="stan/04_river_DO_prior_sensitivity.stan",
  include_2801_DO=FALSE, # PRIMARY: keep original 2801 DO masked
  scenarios=tibble::tibble(
    scenario=c("A_fixed_7p5","C_prior_sd1p5"),
    infer_river_DO=c(0L,1L),
    river_DO_prior_sd=c(0,1.5)
  ),
  folds=5L,
  run_only_fold=NULL,    # e.g. 1L for a one-fold test; NULL for all 5
  resume=TRUE,           # safely reuse existing fit_*.rds files
  refit_failed=FALSE,    # TRUE: refit ONLY cached failed fits
  failed_sampling=3500L, # longer sampling for selective retry
  seed=cfg$seed+5000L,
  chains=cfg$chains,
  parallel_chains=cfg$parallel_chains,
  warmup=cfg$warmup,
  sampling=cfg$sampling,
  adapt_delta=cfg$adapt_delta,
  max_treedepth=cfg$max_treedepth,
  max_rhat=1.05,
  min_bulk_ess=400,
  min_tail_ess=200,
  bootstrap_draws=2000L
)

if (!file.exists(cfg05$stan))
  stop("Missing Stan model: ",cfg05$stan)
if (cfg05$folds<2L)
  stop("CV requires at least two folds.")
if (cfg05$refit_failed &&
    cfg05$failed_sampling <= cfg05$sampling)
  stop("failed_sampling must exceed the original sampling iterations.")
if (!is.null(cfg05$run_only_fold) &&
    (!length(cfg05$run_only_fold) ||
     anyNA(cfg05$run_only_fold) ||
     any(!cfg05$run_only_fold %in% seq_len(cfg05$folds))))
  stop("run_only_fold must be NULL or fold numbers from 1 to ",
       cfg05$folds)
if (!identical(d$DO,d$DO_original))
  stop("Site DO changed in setup; 2801 must retain its actual measurement.")
if (d$DO_original[proxy_index]!=10.2)
  warning("Original 2801 DO differs from expected 10.2 mg/L: ",
          d$DO_original[proxy_index],
          ". Confirm local private data before interpreting sensitivity.")
if (!("2130" %in% d$Site_ID))
  stop("Site 2130 missing. Check the private input.")
if (!isTRUE(all.equal(
  sort(cfg05$scenarios$scenario),
  c("A_fixed_7p5","C_prior_sd1p5")
)))
  stop("Unexpected CV scenarios: update the paired comparison logic.")

suite <- if (cfg05$include_2801_DO) "2801_in_training" else "2801_masked"
cvdir <- file.path("outputs","05_river_DO_grouped_CV",suite)
dir.create(cvdir,recursive=TRUE,showWarnings=FALSE)
out05 <- function(file) file.path(cvdir,file)
scenarios <- cfg05$scenarios

# Score all non-proxy wells with DO and MRT. Site 2801 is NEVER scored,
# including in the "2801_in_training" sensitivity suite.
eligible05 <- which(stan_data$has_DO==1L &
                    d$Site_ID!=cfg$river_proxy_site)
stopifnot(
  length(eligible05)==nrow(d)-1L,
  !any(eligible05==proxy_index),
  all(is.finite(d$DO_original[eligible05]))
)

# Group samples from the same named family, exactly as in 00_setup.R.
# The same fold assignments are used for ALL four model variants,
# including both aquifer kinetic orders.
set.seed(cfg05$seed)
groups <- split(eligible05,d$cv_group[eligible05])
gperm <- sample(seq_along(groups))
gorder <- gperm[order(-lengths(groups)[gperm])]
folds <- vector("list",cfg05$folds)
for (i in gorder) {
  j <- which.min(lengths(folds))
  folds[[j]] <- c(folds[[j]],groups[[i]])
}
if (any(lengths(folds)==0) ||
    !identical(sort(unlist(folds)),eligible05))
  stop("Invalid grouped CV fold construction.")
foldmap <- purrr::map_dfr(seq_along(folds),function(j) {
  tibble::tibble(
    fold=j,sample_index=folds[[j]],
    Site_ID=d$Site_ID[folds[[j]]],
    cv_group=d$cv_group[folds[[j]]]
  )
})
if (any(vapply(split(foldmap$fold,foldmap$cv_group),
               function(f) length(unique(f))>1,logical(1))))
  stop("A CV group has been split between folds.")
if (!all(foldmap$sample_index %in% eligible05))
  stop("2801 would leak into CV scoring.")
stopifnot(nrow(foldmap)==length(eligible05))

xbase <- stan_data
xbase$DO_obs <- as.numeric(d$DO_original)
xbase$has_DO[proxy_index] <- as.integer(cfg05$include_2801_DO)

# Reproducibility/resume guard. Editing priors, input, folds, sampling,
# source code or Stan code requires a new suite directory. Prevent
# accidental reuse of cached fits against different assumptions.
filehash <- as.character(tools::md5sum(c(
  cfg$input,"scripts/00_setup.R",cfg05$stan
)))
if (anyNA(filehash))
  stop("Cannot compute hashes for CV input, setup and Stan model.")
design <- list(
  protocol="05_river_DO_grouped_CV_v1",
  hashes=filehash,
  proxy_training=cfg05$include_2801_DO,
  proxy_id=cfg$river_proxy_site,
  proxy_DO=d$DO_original[proxy_index],
  scenarios=scenarios,
  fold_assignment=foldmap,
  sample_n=nrow(d),
  obs_O18=xbase$O18,
  obs_Cl=xbase$Cl,
  obs_MRT=xbase$mrt_reported,
  obs_DO=xbase$DO_obs,
  seed=cfg05$seed,
  chains=cfg05$chains,
  parallel_chains=cfg05$parallel_chains,
  warmup=cfg05$warmup,
  sampling=cfg05$sampling,
  adapt_delta=cfg05$adapt_delta,
  max_treedepth=cfg05$max_treedepth
)
designfile <- out05("cv_design.rds")
if (file.exists(designfile)) {
  if (!identical(readRDS(designfile),design))
    stop("Existing saved CV fits were generated using a different ",
         "input, setup, Stan model or configuration. ",
         "Do not mix them: move/rename the old outputs/05 suite ",
         "before starting a new experiment.")
} else {
  existing <- list.files(cvdir,pattern="^fit_.*\\.rds$")
  if (length(existing))
    stop("Existing CV fit files without cv_design.rds; cannot safely ",
         "reuse them. Move the old suite directory.")
  saveRDS(design,designfile)
}
readr::write_csv(foldmap,out05("cv_folds.csv"))
readr::write_csv(
  d %>% dplyr::select(
    sample_index,Site_ID,DO_original,MRT,O18,Cl,cv_group
  ) %>%
    dplyr::mutate(
      proxy_2801_DO_used_for_training=cfg05$include_2801_DO,
      selected_for_CV=sample_index %in% eligible05
    ),
  out05("study_samples.csv")
)

message("\n05 GROUPED CV: ",suite,
        "; 40 non-proxy wells in ",cfg05$folds,
        " grouped folds; 4 models; 20 planned fits.")
message("BROAD river prior = Normal(",cfg$river_DO_boundary,
        ", 1.5^2) truncated to [0,14.65] mg/L.")
message("Site 2801 actual DO=",d$DO_original[proxy_index],
        " mg/L; use in TRAINING: ",cfg05$include_2801_DO,
        "; scored: NEVER.")
message("Site 2130 will be held out as part of its assigned fold.")

# Diagnostics: check all sampled latent parameters, source chemistry,
# global DO parameters and the key aggregate LSR loss quantity. The
# river boundary is not a sampled parameter in the fixed scenario.
diagnostic_vars <- function(infer) c(
  "O18_R","O18_LSR","Cl_R","Cl_LSR",
  "f_LSR","z_mrt","sigma_O18","tau_LSR_O18",
  "sigma_Cl","survival_LSR","k_ref","sigma_DO",
  "mean_f_LSR","LSR_DO_deficit","mean_vadose_LSR_DO_loss",
  "aquifer_rate_actual",
  if (infer==1L) "river_DO_free"
)
check05 <- function(fit,scenario,order,fold,infer) {
  summary <- fit$summary(variables=diagnostic_vars(infer))
  hm <- fit$diagnostic_summary()
  vals_ok <- nrow(summary)>0L &&
    all(is.finite(summary$rhat)) &&
    all(is.finite(summary$ess_bulk)) &&
    all(is.finite(summary$ess_tail)) &&
    all(is.finite(hm$ebfmi))
  max_rhat <- max(summary$rhat,na.rm=TRUE)
  min_bulk_ESS <- min(summary$ess_bulk,na.rm=TRUE)
  min_tail_ESS <- min(summary$ess_tail,na.rm=TRUE)
  min_ebfmi <- min(hm$ebfmi,na.rm=TRUE)
  divergences <- sum(hm$num_divergent)
  treedepth_hits <- sum(hm$num_max_treedepth)
  passed <- vals_ok &&
    max_rhat<=cfg05$max_rhat &&
    min_bulk_ESS>=cfg05$min_bulk_ess &&
    min_tail_ESS>=cfg05$min_tail_ess &&
    divergences==0 && treedepth_hits==0 &&
    min_ebfmi>=0.3
  list(
    status=tibble::tibble(
      scenario=scenario,aquifer_order=order,fold=fold,
      max_rhat=max_rhat,min_bulk_ess=min_bulk_ESS,
      min_tail_ess=min_tail_ESS,
      divergences=divergences,
      max_treedepth_hits=treedepth_hits,
      min_ebfmi=min_ebfmi,
      all_diagnostics_pass=passed
    ),
    parameters=summary %>%
      dplyr::mutate(
        scenario=scenario,aquifer_order=order,fold=fold,
        .before=1
      )
  )
}

# Fit object checkpointing + per-fit prediction/diagnostic output mean
# even an interruption after a single fold is recoverable.
model <- cmdstanr::cmdstan_model(cfg05$stan,quiet=TRUE)
selected_folds <- if (is.null(cfg05$run_only_fold)) {
  seq_len(cfg05$folds)
} else sort(unique(as.integer(cfg05$run_only_fold)))

for (fold in selected_folds) {
  held <- folds[[fold]]
  for (sc_i in seq_len(nrow(scenarios))) {
    scenario <- scenarios[sc_i,]
    for (aq_order in 0:1) {
      tag <- paste0(
        "fold",fold,"_",scenario$scenario,"_order",aq_order
      )
      fitfile <- out05(paste0("fit_",tag,".rds"))
      predfile <- out05(paste0("pred_",tag,".csv"))
      statusfile <- out05(paste0("status_",tag,".csv"))
      diagfile <- out05(paste0("parameter_diagnostics_",tag,".csv"))

      x <- xbase
      x$has_DO[held] <- 0L
      # Wipe held-out DO from the data field itself, not just mask.
      x$DO_obs[held] <- 0
      x$aq_order <- as.integer(aq_order)
      x$infer_river_DO <- as.integer(scenario$infer_river_DO)
      x$river_DO_prior_sd <- as.numeric(scenario$river_DO_prior_sd)
      stopifnot(
        all(x$has_DO[held]==0L),
        all(x$DO_obs[held]==0),
        all(x$has_O18==stan_data$has_O18),
        all(x$has_Cl==stan_data$has_Cl),
        !any(held==proxy_index),
        x$has_DO[proxy_index]==as.integer(cfg05$include_2801_DO)
      )
      # Permit a targeted retry with longer chains without changing
      # the base design or discarding the other 19 completed fits.
      # The previous failing RDS remains available as .previous_fail.rds.
      old_status <- if (file.exists(statusfile))
        readr::read_csv(statusfile,show_col_types=FALSE)
      else NULL
      retry_this <- cfg05$refit_failed &&
        file.exists(fitfile) &&
        !is.null(old_status) &&
        nrow(old_status)==1L &&
        !isTRUE(old_status$all_diagnostics_pass[[1]])
      if (file.exists(fitfile) && !retry_this) {
        if (!cfg05$resume)
          stop("Existing fit object ",fitfile,
               " and resume=FALSE: refusing to overwrite.")
        message("REUSE CV FIT: ",tag)
        fit <- readRDS(fitfile)
        iter_used <- if (!is.null(old_status) &&
                         "iter_sampling_used" %in% names(old_status))
          as.integer(old_status$iter_sampling_used[[1]])
        else cfg05$sampling
      } else {
        if (retry_this && !cfg05$resume)
          stop("Selective retry requires resume=TRUE.")
        iter_used <- if (retry_this)
          cfg05$failed_sampling
        else cfg05$sampling
        message(
          if (retry_this) "RETRY FAILED CV FIT: " else "SAMPLE CV: ",
          tag," | held out ",length(held),
          " DO observations | sampling iterations=",iter_used
        )
        fit <- model$sample(
          data=x,
          seed=as.integer(
            cfg05$seed+10001L*fold+1000L*aq_order+37L*sc_i+
              if (retry_this) 1000003L else 0L
          ),
          chains=cfg05$chains,
          parallel_chains=cfg05$parallel_chains,
          iter_warmup=cfg05$warmup,
          iter_sampling=iter_used,
          adapt_delta=cfg05$adapt_delta,
          max_treedepth=cfg05$max_treedepth,
          refresh=250
        )
        # Save BEFORE downstream diagnostics/scoring.
        if (retry_this) {
          staged <- paste0(fitfile,".retry.rds")
          archived <- paste0(fitfile,".previous_fail.rds")
          if (file.exists(staged) || file.exists(archived))
            stop("Retry staging/archive already exists for ",fitfile,
                 "; inspect before overwriting anything.")
          fit$save_object(file=staged)
          if (!file.rename(fitfile,archived))
            stop("Could not archive previous failed fit ",fitfile)
          if (!file.rename(staged,fitfile)) {
            file.rename(archived,fitfile)
            stop("Could not promote retried fit; old fit restored.")
          }
        } else fit$save_object(file=fitfile)
      }

      chk <- check05(
        fit,scenario$scenario,aq_order,fold,
        as.integer(scenario$infer_river_DO)
      )
      chk$status$iter_sampling_used <- iter_used
      chk$status$was_selectively_refit <- retry_this
      readr::write_csv(chk$status,statusfile)
      readr::write_csv(chk$parameters,diagfile)
      if (!chk$status$all_diagnostics_pass)
        warning(
          "MCMC diagnostics require review: ",tag,
          " (results retained; do NOT rank these models yet)"
        )

      mean_vars <- paste0("mean_DO[",held,"]")
      rep_vars <- paste0("replicated_DO[",held,"]")
      frac_vars <- paste0("f_LSR[",held,"]")
      dr <- posterior::as_draws_matrix(
        fit$draws(variables=c(mean_vars,rep_vars,frac_vars))
      )
      col_q <- function(names,q) vapply(names,function(v)
        unname(stats::quantile(as.numeric(dr[,v]),q)),numeric(1))
      means <- colMeans(dr[,mean_vars,drop=FALSE])
      frac <- vapply(frac_vars,function(v)
        median(as.numeric(dr[,v])),numeric(1))
      upper_tail <- vapply(seq_along(held),function(j)
        mean(as.numeric(dr[,rep_vars[j]])>=d$DO_original[held[j]]),
        numeric(1)
      )
      predictions <- tibble::tibble(
        scenario=scenario$scenario,
        aquifer_order=aq_order,
        fold=fold,
        sample_index=held,
        Site_ID=d$Site_ID[held],
        cv_group=d$cv_group[held],
        reported_MRT=d$MRT[held],
        observed_DO=d$DO_original[held],
        posterior_median_LSR_fraction=frac,
        predicted_DO=as.numeric(means),
        fitted_mean_lower95=col_q(mean_vars,.025),
        fitted_mean_upper95=col_q(mean_vars,.975),
        predictive_lower95=col_q(rep_vars,.025),
        predictive_upper95=col_q(rep_vars,.975),
        P_replicated_ge_observed=upper_tail
      ) %>%
        dplyr::mutate(
          residual=predicted_DO-observed_DO,
          absolute_error=abs(residual),
          covered95=observed_DO>=predictive_lower95 &
                    observed_DO<=predictive_upper95
        )
      readr::write_csv(predictions,predfile)
    }
  }
}

# Recover EVERY completed fold from its per-fit checkpoint (including
# folds completed during earlier partial runs). Never silently treat
# incomplete 1-fold tests as full five-fold cross-validation.
all_tags <- tidyr::expand_grid(
  fold=seq_len(cfg05$folds),
  scenario=scenarios$scenario,
  aquifer_order=0:1
) %>%
  dplyr::mutate(
    tag=paste0(
      "fold",fold,"_",scenario,"_order",aquifer_order
    ),
    pred_path=out05(paste0("pred_",tag,".csv")),
    status_path=out05(paste0("status_",tag,".csv"))
  )
have <- file.exists(all_tags$pred_path) &
        file.exists(all_tags$status_path)
completed <- all_tags[have,,drop=FALSE]
if (!nrow(completed))
  stop("No CV predictions available; check the selected folds.")

read_all <- function(paths) dplyr::bind_rows(lapply(
  paths,function(p) readr::read_csv(p,show_col_types=FALSE)
))
preds <- read_all(completed$pred_path)
statuses <- read_all(completed$status_path)
readr::write_csv(preds,out05("cv_predictions.csv"))
readr::write_csv(statuses,out05("cv_fit_status.csv"))

if (anyDuplicated(preds %>%
                  dplyr::select(scenario,aquifer_order,sample_index)))
  stop("Duplicate held-out predictions; inspect cached fold files.")
if (any(preds$Site_ID==cfg$river_proxy_site))
  stop("2801 must never be CV-scored.")
if (anyNA(preds$observed_DO) ||
    anyNA(preds$predicted_DO))
  stop("Missing held-out predictions.")

scores <- preds %>%
  dplyr::group_by(scenario,aquifer_order) %>%
  dplyr::summarise(
    n=n(),folds_completed=dplyr::n_distinct(fold),
    RMSE=sqrt(mean(residual^2)),
    MAE=mean(abs(residual)),
    bias=mean(residual),
    coverage95=mean(covered95),
    .groups="drop"
  )
# Readr may load sample_index as double while eligible05 is integer.
# identical() compares STORAGE TYPES and was incorrectly marking all
# folds incomplete despite valid 40/40 unique held-out predictions.
# Use an integer-normalised set check with counts and no duplication.
coverage <- preds %>%
  dplyr::group_by(scenario,aquifer_order) %>%
  dplyr::summarise(
    all_folds_complete=
      (dplyr::n()==length(eligible05)) &&
      (dplyr::n_distinct(fold)==cfg05$folds) &&
      !anyDuplicated(as.integer(sample_index)) &&
      setequal(as.integer(sample_index),eligible05),
    .groups="drop"
  )
quality <- statuses %>%
  dplyr::group_by(scenario,aquifer_order) %>%
  dplyr::summarise(
    all_MCMC_pass=
      (dplyr::n()==cfg05$folds) &&
      (dplyr::n_distinct(fold)==cfg05$folds) &&
      all(all_diagnostics_pass %in% TRUE),
    .groups="drop"
  )
scores <- scores %>%
  dplyr::left_join(coverage,by=c("scenario","aquifer_order")) %>%
  dplyr::left_join(quality,by=c("scenario","aquifer_order")) %>%
  dplyr::mutate(
    all_folds_complete=dplyr::coalesce(all_folds_complete,FALSE),
    all_MCMC_pass=dplyr::coalesce(all_MCMC_pass,FALSE)
  )
complete <- nrow(completed)==cfg05$folds*2L*nrow(scenarios) &&
  all(scores$all_folds_complete)
all_pass <- nrow(statuses)==nrow(completed) &&
  all(scores$all_MCMC_pass)
readr::write_csv(scores,out05("cv_scores.csv"))
cat("\n05 GROUPED CV SCORES (read flags before comparison)\n")
print(scores,n=Inf,width=Inf)

# Highlight the 2130 diagnostic using HELD-OUT predictions only.
site2130 <- preds %>% dplyr::filter(Site_ID=="2130")
if (nrow(site2130)>0) {
  readr::write_csv(site2130,out05("cv_2130_predictions.csv"))
  cat("\nSITE 2130 — GENUINELY HELD-OUT DO (not in likelihood)\n")
  print(site2130 %>% dplyr::select(
    scenario,aquifer_order,fold,observed_DO,predicted_DO,
    predictive_lower95,predictive_upper95,
    P_replicated_ge_observed
  ),n=Inf,width=Inf)
}

if (complete && all_pass) {
  # Construct 40 fully paired held-out observations, without imputation.
  paired <- preds %>% dplyr::select(
    sample_index,Site_ID,cv_group,reported_MRT,observed_DO,
    scenario,aquifer_order,predicted_DO
  ) %>% dplyr::mutate(
    model=paste0(scenario,"_order",aquifer_order)
  ) %>% dplyr::select(-scenario,-aquifer_order) %>%
    tidyr::pivot_wider(names_from=model,values_from=predicted_DO)
  expected_cols <- c(
    "A_fixed_7p5_order0","A_fixed_7p5_order1",
    "C_prior_sd1p5_order0","C_prior_sd1p5_order1"
  )
  if (!all(expected_cols %in% names(paired)) ||
      nrow(paired)!=length(eligible05) ||
      anyNA(paired[,expected_cols]))
    stop("Paired model predictions are incomplete or misaligned.")
  readr::write_csv(paired,out05("cv_paired_predictions.csv"))

  comparisons <- tibble::tibble(
    comparison=c(
      "C_minus_A_zero_order","C_minus_A_first_order",
      "first_minus_zero_fixed","first_minus_zero_broad"
    ),
    candidate=c(
      "C_prior_sd1p5_order0","C_prior_sd1p5_order1",
      "A_fixed_7p5_order1","C_prior_sd1p5_order1"
    ),
    comparator=c(
      "A_fixed_7p5_order0","A_fixed_7p5_order1",
      "A_fixed_7p5_order0","C_prior_sd1p5_order0"
    )
  )
  set.seed(cfg05$seed+909L)
  grouped <- split(paired,paired$cv_group)
  boot <- purrr::map_dfr(seq_len(nrow(comparisons)),function(i) {
    cmp <- comparisons[i,]
    purrr::map_dfr(seq_len(cfg05$bootstrap_draws),function(j) {
      chosen <- sample(
        seq_along(grouped),length(grouped),replace=TRUE
      )
      x <- dplyr::bind_rows(grouped[chosen])
      errA <- x[[cmp$candidate]]-x$observed_DO
      errB <- x[[cmp$comparator]]-x$observed_DO
      tibble::tibble(
        comparison=cmp$comparison,draw=j,
        delta_RMSE=sqrt(mean(errA^2))-sqrt(mean(errB^2)),
        delta_MAE=mean(abs(errA))-mean(abs(errB))
      )
    })
  })
  readr::write_csv(boot,out05("cv_group_bootstrap.csv"))
  boot_summary <- boot %>%
    dplyr::group_by(comparison) %>%
    dplyr::summarise(
      median_delta_RMSE=median(delta_RMSE),
      RMSE_lower95=unname(quantile(delta_RMSE,.025)),
      RMSE_upper95=unname(quantile(delta_RMSE,.975)),
      median_delta_MAE=median(delta_MAE),
      MAE_lower95=unname(quantile(delta_MAE,.025)),
      MAE_upper95=unname(quantile(delta_MAE,.975)),
      .groups="drop"
    )
  readr::write_csv(
    boot_summary,out05("cv_group_bootstrap_summary.csv")
  )
  cat("\nGROUPED BOOTSTRAP PAIRED DIFFERENCES\n")
  cat("Negative first-minus-zero favours first-order;",
      " negative C-minus-A favours broad river DO.\n")
  print(boot_summary,n=Inf,width=Inf)
} else {
  warning(
    "CV INCOMPLETE OR AT LEAST ONE SAMPLER FAILED DIAGNOSTICS: ",
    "scores are preliminary. No paired bootstrap/model ranking ",
    "is produced. Do not compare models until all 20 fits pass."
  )
}
cat("\nCV output folder: ",normalizePath(cvdir),"\n",sep="")
cat("Plotting only: source(\"scripts/02_DO_plotting.R\")\n")
