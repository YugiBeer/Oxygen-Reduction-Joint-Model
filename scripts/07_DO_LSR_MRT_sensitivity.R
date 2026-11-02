# 07_DO_LSR_MRT_sensitivity.R
#
# Targeted ONE-FACTOR-AT-A-TIME sensitivity analysis of the FINAL regional
# Heretaunga DO model. No changes to Stan, 00_setup, 04, or 05: their saved
# fit objects and the CV design fingerprint remain valid.
#
# PRIMARY BASELINE: 39 DO observations; 2801 and 2130 DO masked;
# delta18O and Cl at both sites remain in mixing likelihood;
# first-order aquifer kinetics; river-entry DO has a truncated
# Normal(mean=7.5, SD=1.5) mg/L prior.
#
# One sensitivity factor at a time:
#   * LSR surface DO saturation assumed at T = 10 or 20 C (baseline 15 C)
#   * LSR fraction surviving to WT ~ Beta(2,2), Beta(2,5), Beta(5,2)
#     (baseline Beta(1,1); these are deliberately contrasting
#      prior-sensitivity probes, NOT independently measured recharge priors)
#   * uncertainty in log groundwater MRT: SD=0.10 or 0.50
#     (baseline 0.25; each source still shares the same effective aquifer MRT)
#
# Baseline is reused from 04; SEVEN new MCMC fits, four chains each.
# Run in Positron at repository root:
#   source("scripts/07_DO_LSR_MRT_sensitivity.R")
#
# Outputs: outputs/07_DO_LSR_MRT_sensitivity/2801_masked_no2130/
#          plots/07_DO_LSR_MRT_sensitivity/2801_masked_no2130/
#
# Inference caveats:
# - Parameter sensitivity is not a new cross-validation exercise.
# - Posterior intervals are not confidence intervals for differences
#   between models fitted to the same data.
# - This analysis cannot resolve source-specific residence times or
#   progressive addition of LSR along groundwater flowpaths.

if (!file.exists("scripts/00_setup.R"))
  stop("Open the repository root as the R working directory.")
source("scripts/00_setup.R", local=TRUE)

required07 <- c("cmdstanr","posterior","ggplot2","dplyr","tidyr",
                "tibble","readr","viridisLite")
missing07 <- required07[!vapply(required07,requireNamespace,
                               logical(1),quietly=TRUE)]
if (length(missing07))
  stop("Install required R packages: ",paste(missing07,collapse=", "))

# ------------------- CONTROLS ------------------------------------------
cfg07 <- list(
  stan="stan/04_river_DO_prior_sensitivity.stan",
  baseline_dir="outputs/04_river_DO_prior/2801_masked_no2130",
  output_dir="outputs/07_DO_LSR_MRT_sensitivity/2801_masked_no2130",
  plot_dir="plots/07_DO_LSR_MRT_sensitivity/2801_masked_no2130",
  resume=TRUE,             # Reuse completed and matching fits
  refit_failed=FALSE,      # TRUE reruns cached fits that fail diagnostics
  run_only=NULL,           # NULL = all 7; or c("temp_10C") for a test
  warmup=cfg$warmup,
  sampling=cfg$sampling,
  chains=cfg$chains,
  parallel_chains=cfg$parallel_chains,
  adapt_delta=cfg$adapt_delta,
  max_treedepth=cfg$max_treedepth,
  seed=cfg$seed+7000L,
  max_rhat=1.05,
  min_bulk_ess=400,
  min_tail_ess=200,
  min_key_tail_ess=400
)
# -----------------------------------------------------------------------

cases07 <- tibble::tribble(
  ~scenario,            ~factor,          ~value_label,
  ~temp_C, ~beta_a, ~beta_b, ~mrt_logSD,
  "baseline",            "Baseline",        "15 C; Beta(1,1); MRT logSD 0.25",
  15,      1,       1,       0.25,
  "temp_10C",            "Surface DO",      "Saturation at 10 C",
  10,      1,       1,       0.25,
  "temp_20C",            "Surface DO",      "Saturation at 20 C",
  20,      1,       1,       0.25,
  "survival_beta2_2",    "LSR prior",        "Survival Beta(2,2)",
  15,      2,       2,       0.25,
  "survival_beta2_5",    "LSR prior",        "Survival Beta(2,5)",
  15,      2,       5,       0.25,
  "survival_beta5_2",    "LSR prior",        "Survival Beta(5,2)",
  15,      5,       2,       0.25,
  "mrt_logSD_0p10",      "Age uncertainty", "MRT logSD 0.10",
  15,      1,       1,       0.10,
  "mrt_logSD_0p50",      "Age uncertainty", "MRT logSD 0.50",
  15,      1,       1,       0.50
) %>% dplyr::mutate(
  Csat_rain=vapply(temp_C,do_sat_mgL,numeric(1),
                   pressure_ratio=cfg$atmospheric_pressure_ratio),
  scenario_order=dplyr::row_number()
)

if (anyDuplicated(cases07$scenario))
  stop("Scenario labels must be unique.")
if (!is.null(cfg07$run_only) &&
    !all(cfg07$run_only %in% cases07$scenario[-1]))
  stop("Unknown scenario in cfg07$run_only.")
if (!identical(cfg$river_DO_boundary,7.5) ||
    !isTRUE(all.equal(cfg$log_mrt_SD,0.25)) ||
    !isTRUE(all.equal(cfg$recharge_temperature_C,15)) ||
    !identical(cfg$survival_prior_a,1) ||
    !identical(cfg$survival_prior_b,1))
  stop("00_setup baseline settings changed. Reconcile with 04 baseline.")
if (!file.exists(cfg07$stan))
  stop("Stan file missing: ",cfg07$stan)

dir.create(cfg07$output_dir,recursive=TRUE,showWarnings=FALSE)
dir.create(cfg07$plot_dir,recursive=TRUE,showWarnings=FALSE)
of07 <- function(name) file.path(cfg07$output_dir,name)
pf07 <- function(name) file.path(cfg07$plot_dir,name)

# Build exactly the same 39-well DO mask used in script 04.
excluded07 <- which(d$Site_ID %in% c("2801","2130"))
if (length(excluded07)!=2L ||
    !setequal(as.character(d$Site_ID[excluded07]),c("2801","2130")))
  stop("Expected one 2801 and one 2130 site.")
xbase07 <- stan_data
xbase07$DO_obs <- as.numeric(d$DO_original)
xbase07$has_DO[excluded07] <- 0L
xbase07$DO_obs[excluded07] <- 0
xbase07$aq_order <- 1L
xbase07$infer_river_DO <- 1L
xbase07$river_DO_prior_sd <- 1.5
if (sum(xbase07$has_DO)!=39L ||
    !identical(xbase07$has_O18,stan_data$has_O18) ||
    !identical(xbase07$has_Cl,stan_data$has_Cl))
  stop("Unexpected observation mask: expected 39 DO plus retained tracers.")

# Verify baseline provenance rather than silently comparing different models.
baseline_fit07 <- file.path(cfg07$baseline_dir,
                            "fit_C_prior_sd1p5_order1.rds")
baseline_status07 <- file.path(cfg07$baseline_dir,"fit_status.csv")
baseline_design07 <- file.path(cfg07$baseline_dir,"scenario_design.csv")
need_base07 <- c(baseline_fit07,baseline_status07,baseline_design07)
if (any(!file.exists(need_base07)))
  stop("Missing original 39-well baseline from script 04:\n",
       paste(need_base07[!file.exists(need_base07)],collapse="\n"),
       "\nDo not rerun 04 unless those fitted outputs are actually missing.")

bd07 <- readr::read_csv(baseline_design07,show_col_types=FALSE)
br07 <- bd07 %>% dplyr::filter(
  .data$scenario=="C_prior_sd1p5",
  .data$infer_river_DO==1L,
  .data$river_DO_prior_sd==1.5
)
if (nrow(br07)!=1L ||
    !identical(as.integer(br07$n_independent_DO),39L) ||
    !identical(as.character(br07$excluded_DO_sites),"2130") ||
    !identical(as.logical(br07$proxy_DO_used),FALSE) ||
    !isTRUE(all.equal(as.numeric(br07$rainfall_surface_DO),
                      as.numeric(cases07$Csat_rain[1]),tolerance=1e-8)) ||
    !isTRUE(all.equal(as.numeric(br07$log_mrt_SD),0.25)))
  stop("Script-04 baseline design mismatch; aborting comparison.")

bs07 <- readr::read_csv(baseline_status07,show_col_types=FALSE)
bsrow07 <- bs07 %>% dplyr::filter(
  .data$scenario=="C_prior_sd1p5", .data$aquifer_order==1L)
if (nrow(bsrow07)!=1L || !isTRUE(bsrow07$all_diagnostics_pass))
  stop("Script-04 baseline failed diagnostics or has no pass record.")

# Fingerprint the actual model, setup and input data. The sensitivity
# design fingerprint is separate from script 05's existing cache!
hashes07 <- tools::md5sum(c(cfg$input,"scripts/00_setup.R",cfg07$stan))
if (anyNA(hashes07)) stop("Cannot fingerprint input/setup/Stan.")
base_hash07 <- unname(tools::md5sum(baseline_fit07))
if (is.na(base_hash07)) stop("Cannot fingerprint baseline fit.")
design07 <- list(
  protocol="07_DO_LSR_MRT_sensitivity_v1",
  hashes=as.character(hashes07),
  baseline_fit_hash=base_hash07,
  cases=cases07,
  DO_mask=xbase07$has_DO,
  tracer_O18=xbase07$O18,
  tracer_Cl=xbase07$Cl,
  MRT=xbase07$mrt_reported,
  river_DO_prior_mean=cfg$river_DO_boundary,
  river_DO_prior_sd=1.5,
  aq_order=1L,
  warmup=cfg07$warmup,sampling=cfg07$sampling,chains=cfg07$chains,
  adapt_delta=cfg07$adapt_delta,max_treedepth=cfg07$max_treedepth,
  seed=cfg07$seed
)
design_path07 <- of07("design.rds")
if (file.exists(design_path07)) {
  if (!identical(readRDS(design_path07),design07))
    stop("Existing sensitivity fits were generated from another design.",
         " Move the old output folder before changing assumptions.")
} else {
  cached07 <- list.files(cfg07$output_dir,pattern="^fit_.*[.]rds$")
  if (length(cached07))
    stop("Cached fits found without design.rds: cannot safely reuse.")
  saveRDS(design07,design_path07)
}
readr::write_csv(cases07,of07("scenario_design.csv"))

monitor07 <- c(
  "O18_R","O18_LSR","Cl_R","Cl_LSR","f_LSR","z_mrt",
  "sigma_O18","tau_LSR_O18","sigma_Cl","survival_LSR",
  "k_ref","sigma_DO","river_DO_free"
)
key07 <- c("aquifer_rate_actual","LSR_DO_at_water_table",
           "LSR_DO_deficit","river_DO","sigma_DO","mean_f_LSR")
param07 <- c(key07,"survival_LSR","mean_vadose_LSR_DO_loss",
             "mean_aquifer_DO_loss")

check07 <- function(fit,sc) {
  sm <- fit$summary(variables=monitor07)
  keysm <- fit$summary(variables=key07)
  hm <- fit$diagnostic_summary()
  div <- sum(hm$num_divergent)
  tree <- sum(hm$num_max_treedepth)
  min_ebfmi <- min(hm$ebfmi)
  max_rhat <- max(sm$rhat)
  min_bulk <- min(sm$ess_bulk)
  min_tail <- min(sm$ess_tail)
  min_key_tail <- min(keysm$ess_tail)
  valid <- all(is.finite(c(
    sm$rhat,sm$ess_bulk,sm$ess_tail,hm$ebfmi,
    keysm$rhat,keysm$ess_bulk,keysm$ess_tail
  )))
  passed <- valid &&
    div==0 && tree==0 &&
    min_ebfmi>=0.3 &&
    max_rhat<=cfg07$max_rhat &&
    min_bulk>=cfg07$min_bulk_ess &&
    min_tail>=cfg07$min_tail_ess &&
    min_key_tail>=cfg07$min_key_tail_ess
  tibble::tibble(
    scenario=sc,all_diagnostics_pass=passed,
    max_rhat=max_rhat,min_bulk_ess=min_bulk,
    min_tail_ess=min_tail,min_key_tail_ess=min_key_tail,
    divergences=div,max_treedepth_hits=tree,min_ebfmi=min_ebfmi
  )
}

summarise07 <- function(fit,case_row) {
  z <- posterior::as_draws_matrix(fit$draws(variables=param07))
  tibble::tibble(
    scenario=case_row$scenario,
    factor=case_row$factor,
    value_label=case_row$value_label,
    variable=param07,
    median=vapply(param07,function(v) median(z[,v]),numeric(1)),
    mean=vapply(param07,function(v) mean(z[,v]),numeric(1)),
    sd=vapply(param07,function(v) stats::sd(z[,v]),numeric(1)),
    lower95=vapply(param07,function(v)
      unname(stats::quantile(z[,v],0.025)),numeric(1)),
    upper95=vapply(param07,function(v)
      unname(stats::quantile(z[,v],0.975)),numeric(1))
  )
}

# Use saved baseline, without re-running script 04.
baseline07 <- readRDS(baseline_fit07)
status_rows07 <- list(check07(baseline07,"baseline"))
summary_rows07 <- list(summarise07(baseline07,cases07[1,]))
if (!status_rows07[[1]]$all_diagnostics_pass)
  warning("Baseline fails the stricter script-07 tail-ESS checks.",
          " Review before claiming robust inference.")

todo07 <- cases07 %>% dplyr::filter(.data$scenario!="baseline")
if (!is.null(cfg07$run_only))
  todo07 <- todo07 %>%
    dplyr::filter(.data$scenario %in% cfg07$run_only)

# Compile only when a new/refitted scenario is needed.
mod07 <- NULL
for (j in seq_len(nrow(todo07))) {
  sc <- todo07[j,]
  fit_path <- of07(paste0("fit_",sc$scenario,".rds"))
  rerun <- !file.exists(fit_path) || !cfg07$resume

  if (!rerun && cfg07$refit_failed) {
    cached_fit <- readRDS(fit_path)
    rerun <- !isTRUE(check07(cached_fit,sc$scenario)$all_diagnostics_pass)
    rm(cached_fit)
  }

  if (rerun) {
    if (is.null(mod07)) {
      message("Compiling model for first-order sensitivity fits...")
      mod07 <- cmdstanr::cmdstan_model(cfg07$stan,quiet=TRUE)
    }
    x <- xbase07
    x$Csat_rain <- as.numeric(sc$Csat_rain)
    x$beta_survival_a <- as.numeric(sc$beta_a)
    x$beta_survival_b <- as.numeric(sc$beta_b)
    x$log_mrt_SD <- as.numeric(sc$mrt_logSD)

    message("\nFIT ",sc$scenario,": ",sc$value_label,
            " | river prior SD=1.5 | first-order kinetics")
    fit <- mod07$sample(
      data=x,
      chains=cfg07$chains,
      parallel_chains=cfg07$parallel_chains,
      iter_warmup=cfg07$warmup,
      iter_sampling=cfg07$sampling,
      adapt_delta=cfg07$adapt_delta,
      max_treedepth=cfg07$max_treedepth,
      seed=as.integer(cfg07$seed+1009L*sc$scenario_order),
      refresh=250
    )
    # Write only AFTER sampling completes; old results remain safe
    # if an interrupted run fails before $save_object().
    if (file.exists(fit_path)) {
      backup <- sub("[.]rds$",".previous.rds",fit_path)
      file.copy(fit_path,backup,overwrite=TRUE)
    }
    fit$save_object(file=fit_path)
  } else {
    message("REUSE ",sc$scenario,": ",fit_path)
    fit <- readRDS(fit_path)
  }

  status_rows07[[length(status_rows07)+1L]] <-
    check07(fit,sc$scenario)
  summary_rows07[[length(summary_rows07)+1L]] <-
    summarise07(fit,sc)
  readr::write_csv(dplyr::bind_rows(status_rows07),
                   of07("fit_status.csv"))
  readr::write_csv(dplyr::bind_rows(summary_rows07),
                   of07("posterior_summary.csv"))
  if (!isTRUE(status_rows07[[length(status_rows07)]]$all_diagnostics_pass))
    warning("MCMC diagnostics need review for ",sc$scenario,
            "; conclusions should be provisional.")
  rm(fit)
}

status07 <- dplyr::bind_rows(status_rows07)
post07 <- dplyr::bind_rows(summary_rows07)
readr::write_csv(status07,of07("fit_status.csv"))
readr::write_csv(post07,of07("posterior_summary.csv"))

# Only make comparisons for diagnostic-approved fits. A failed run is
# retained for troubleshooting, not ranked against the primary model.
passed07 <- status07 %>%
  dplyr::filter(.data$all_diagnostics_pass) %>%
  dplyr::pull(.data$scenario)

rate07 <- post07 %>%
  dplyr::filter(.data$variable=="aquifer_rate_actual") %>%
  dplyr::left_join(
    status07 %>% dplyr::select(.data$scenario,.data$all_diagnostics_pass),
    by="scenario"
  ) %>%
  dplyr::mutate(
    delta_rate_pct=100*(.data$median/.data$median[
      .data$scenario=="baseline"]-1),
    baseline_median=.data$median[.data$scenario=="baseline"]
  )
readr::write_csv(rate07,of07("rate_sensitivity.csv"))

message("\nAQUIFER FIRST-ORDER RATE SENSITIVITY:")
print(rate07 %>% dplyr::select(.data$scenario,.data$median,
                              .data$lower95,.data$upper95,
                              .data$delta_rate_pct,
                              .data$all_diagnostics_pass),
      n=Inf,width=Inf)

if (!all(status07$all_diagnostics_pass)) {
  warning("Some MCMC checks failed; figures exclude flagged fits.")
}
good07 <- post07 %>% dplyr::filter(.data$scenario %in% passed07)
good_rate07 <- rate07 %>%
  dplyr::filter(.data$all_diagnostics_pass)
if (!nrow(good_rate07))
  stop("No diagnostic-approved models available to plot.")

# Turbo-compatible palette; grouped figures use scientific categories.
turbo07 <- viridisLite::turbo(256)
colors07 <- c(
  "Baseline"=turbo07[128],
  "Surface DO"=turbo07[230],
  "LSR prior"=turbo07[175],
  "Age uncertainty"=turbo07[32]
)
ordered07 <- cases07$scenario[cases07$scenario %in% good_rate07$scenario]
label_lookup07 <- setNames(cases07$value_label,cases07$scenario)
factor_lookup07 <- setNames(cases07$factor,cases07$scenario)
good07 <- good07 %>%
  dplyr::mutate(
    scenario=factor(.data$scenario,levels=rev(ordered07)),
    label=factor(.data$value_label,
                 levels=rev(cases07$value_label[cases07$scenario %in%
                                                 ordered07]))
  )
good_rate07 <- good_rate07 %>%
  dplyr::mutate(
    scenario=factor(.data$scenario,levels=rev(ordered07)),
    label=factor(.data$value_label,
                 levels=rev(cases07$value_label[cases07$scenario %in%
                                                 ordered07]))
  )

theme07 <- function() {
  ggplot2::theme_bw(base_size=12)+
    ggplot2::theme(
      panel.grid.minor=ggplot2::element_blank(),
      plot.title=ggplot2::element_text(face="bold"),
      legend.position="bottom"
    )
}
save07 <- function(p,name,width=10,height=6) {
  print(p)
  ggplot2::ggsave(pf07(name),plot=p,width=width,height=height,
                  dpi=300,bg="white",units="in")
  message("PLOT: ",normalizePath(pf07(name)))
}

base_k07 <- rate07$median[rate07$scenario=="baseline"]
p_rate07 <- ggplot2::ggplot(
  good_rate07,ggplot2::aes(x=.data$median,y=.data$label,
                          colour=.data$factor)
) +
  ggplot2::geom_vline(xintercept=base_k07,
                      linetype="dashed",colour="grey50") +
  ggplot2::geom_segment(
    ggplot2::aes(x=.data$lower95,xend=.data$upper95,
                 yend=.data$label),linewidth=1.15
  ) +
  ggplot2::geom_point(size=2.7) +
  ggplot2::scale_colour_manual(values=colors07,name="Assumption") +
  ggplot2::labs(
    title="First-order aquifer DO reaction-rate sensitivity",
    subtitle="Four-chain posterior medians and 95% credible intervals; dashed = baseline median",
    x="Aquifer first-order reaction rate (year^-1)",y=NULL,
    caption="One assumption varied at a time; river DO retains its Normal(7.5,1.5^2) prior."
  ) + theme07()
save07(p_rate07,"01_aquifer_rate_95CrI.png",width=10,height=6.5)

p_delta07 <- ggplot2::ggplot(
  good_rate07,ggplot2::aes(x=.data$delta_rate_pct,y=.data$label,
                          colour=.data$factor)
) +
  ggplot2::geom_vline(xintercept=0,colour="grey50",
                      linetype="dashed") +
  ggplot2::geom_point(size=3) +
  ggplot2::scale_colour_manual(values=colors07,name="Assumption") +
  ggplot2::labs(
    title="Shift in the aquifer-rate posterior median",
    subtitle="Descriptive percentage difference from baseline; not a credible interval for model differences",
    x="Change in median first-order rate (%)",y=NULL
  ) + theme07()
save07(p_delta07,"02_change_in_median_aquifer_rate.png",
       width=10,height=6.5)

p_wt07 <- ggplot2::ggplot(
  good07 %>% dplyr::filter(.data$variable=="LSR_DO_at_water_table"),
  ggplot2::aes(x=.data$median,y=.data$label,
               colour=.data$factor)
) +
  ggplot2::geom_segment(
    ggplot2::aes(x=.data$lower95,xend=.data$upper95,yend=.data$label),
    linewidth=1.1
  ) +
  ggplot2::geom_point(size=2.7) +
  ggplot2::scale_colour_manual(values=colors07,name="Assumption") +
  ggplot2::labs(
    title="Sensitivity of LSR DO at the water table",
    subtitle="Posterior medians and 95% credible intervals",
    x="LSR water-table DO (mg/L)",y=NULL
  ) + theme07()
save07(p_wt07,"03_LSR_water_table_DO_95CrI.png",width=10,height=6.5)

p_def07 <- ggplot2::ggplot(
  good07 %>% dplyr::filter(.data$variable=="LSR_DO_deficit"),
  ggplot2::aes(x=.data$median,y=.data$label,colour=.data$factor)
) +
  ggplot2::geom_segment(
    ggplot2::aes(x=.data$lower95,xend=.data$upper95,yend=.data$label),
    linewidth=1.1
  ) +
  ggplot2::geom_point(size=2.7) +
  ggplot2::scale_colour_manual(values=colors07,name="Assumption") +
  ggplot2::labs(
    title="Sensitivity of the LSR oxygen deficit",
    subtitle="Deficit is relative to each scenario's assumed land-surface saturation DO",
    x="Inferred LSR deficit (mg/L)",y=NULL
  ) + theme07()
save07(p_def07,"04_LSR_deficit_95CrI.png",width=10,height=6.5)

# The full grid is meaningful only when every scenario has finished.
expected07 <- setdiff(cases07$scenario,"baseline")
complete07 <- all(expected07 %in% status07$scenario)
all_pass07 <- complete07 && all(status07$all_diagnostics_pass)
if (all_pass07) {
  delta07 <- good_rate07 %>%
    dplyr::filter(as.character(.data$scenario)!="baseline") %>%
    dplyr::mutate(abs_change_pct=abs(.data$delta_rate_pct))
  max_shift07 <- max(delta07$abs_change_pct)
  max_case07 <- delta07 %>%
    dplyr::slice_max(.data$abs_change_pct,n=1,with_ties=FALSE)
  message("\nALL SEVEN SENSITIVITY FITS PASSED.")
  message("Maximum absolute change in first-order rate median: ",
          signif(max_shift07,3),"% (",max_case07$scenario,").")
  message("This is conditional sensitivity, not proof of structural identifiability.")
} else {
  warning("Study is incomplete and/or diagnostics failed.",
          " Do not report a final robustness conclusion yet.")
}
message("\nSAVED: ",normalizePath(cfg07$output_dir))
message("FIGURES: ",normalizePath(cfg07$plot_dir))
message("Source-specific recharge ages remain an UNTESTED structural limitation.")
