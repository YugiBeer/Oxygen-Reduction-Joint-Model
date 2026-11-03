# 07_export_tracer_only_fractions.R
#
# Estimate groundwater WATER fractions using O18 and Cl alone.
# Neither DO nor MRT nor excess N2 enters the likelihood.
# Appends posterior medians and 95% credible intervals to a COPY
# of the original input CSV; leaves the original file unchanged.
#
# From repository root:
#   source("scripts/07_export_tracer_only_fractions.R")

# ---- 1. Input and output paths ----------------------------------------
input_file <- "data/Uwe_HP_data_prepared.csv"
stan_file <- "stan/05_tracer_only_mixing.stan"
output_file <- "data/Uwe_HP_data_prepared_with_tracer_fractions.csv"
fit_dir <- "outputs/08_tracer_only_mixing"
fit_file <- file.path(fit_dir,"tracer_only_fit.rds")
design_file <- file.path(fit_dir,"fit_design.rds")
diagnostic_file <- file.path(fit_dir,"tracer_only_diagnostics.csv")
method_file <- file.path(fit_dir,"tracer_only_method.txt")

if (!file.exists(input_file) || !file.exists(stan_file))
  stop("Run from the repository root; input CSV or Stan file missing.")
for (pkg in c("readr","cmdstanr","posterior"))
  if (!requireNamespace(pkg,quietly=TRUE))
    stop("Please install R package: ",pkg)

# ---- 2. Read original data without changing its values ---------------
# Retain original data as character columns, including site identifiers,
# extra metadata and the literal NA values in the original CSV.
dat <- readr::read_csv(
  input_file,na=character(),
  col_types=readr::cols(.default=readr::col_character()),
  show_col_types=FALSE,progress=FALSE
)
if (!all(c("Site_ID","O18","Cl") %in% names(dat)))
  stop("Source CSV must include Site_ID, O18 and Cl.")
if (anyNA(dat$Site_ID) || any(!nzchar(dat$Site_ID)))
  stop("Blank Site_ID found in original input.")

new_columns <- c(
  "LSR_fraction_median","LSR_fraction_lower95","LSR_fraction_upper95",
  "river_fraction_median","river_fraction_lower95","river_fraction_upper95",
  "mixing_tracers_used","mixing_status"
)
if (any(new_columns %in% names(dat)))
  stop("Input has fraction columns already. Use the unmodified source CSV.")

O18 <- suppressWarnings(as.numeric(dat$O18))
Cl <- suppressWarnings(as.numeric(dat$Cl))
# Missing values may be written as blanks or literal NA.
missing_O18 <- trimws(dat$O18) %in% c("","NA")
missing_Cl <- trimws(dat$Cl) %in% c("","NA")
if (any(!missing_O18 & !is.finite(O18)) ||
    any(!missing_Cl & !is.finite(Cl)))
  stop("Unrecognised nonnumeric O18 or Cl value. Check the input CSV.")
if (any(Cl[is.finite(Cl)] < 0))
  stop("Negative chloride value in source CSV.")

has_O18 <- is.finite(O18)
has_Cl <- is.finite(Cl)
excluded <- dat$Site_ID=="16758"
selected <- which(!excluded & (has_O18 | has_Cl))
if (length(selected)<5)
  stop("Too few sites with an available conservative tracer.")

# Fit by ROW INDEX, not Site_ID: some source IDs occur more than once.
# Include 2801 and 2130 (their O18 and Cl are useful), and include
# samples lacking DO or MRT (neither quantity is needed for mixing).
stan_data <- list(
  N=length(selected),
  has_O18=as.integer(has_O18[selected]),
  has_Cl=as.integer(has_Cl[selected]),
  O18=ifelse(has_O18[selected],O18[selected],0),
  Cl=ifelse(has_Cl[selected],Cl[selected],0),
  Cl_sd=pmax(1,0.03*ifelse(has_Cl[selected],Cl[selected],0))
)

# ---- 3. Fit standalone tracer-only model; reuse exact matching fit ----
dir.create(fit_dir,recursive=TRUE,showWarnings=FALSE)
design <- list(
  input_md5=unname(tools::md5sum(input_file)),
  stan_md5=unname(tools::md5sum(stan_file)),
  selected_rows=selected,
  chains=4L,warmup=1400L,sampling=1400L,
  adapt_delta=0.98,max_treedepth=13L,seed=20261104L
)
if (xor(file.exists(fit_file),file.exists(design_file)))
  stop("Incomplete cache: inspect fit/design files in ",fit_dir)

if (file.exists(fit_file)) {
  if (!identical(readRDS(design_file),design))
    stop("Data, Stan model or configuration changed. Move the old ",
         "outputs/08_tracer_only_mixing folder before refitting.")
  message("Reusing saved tracer-only fit.")
  fit <- readRDS(fit_file)
} else {
  message("Fitting independent tracer-only O18 + Cl mixing model.")
  cmdstanr::cmdstan_path()
  model <- cmdstanr::cmdstan_model(stan_file,quiet=TRUE)
  fit <- model$sample(
    data=stan_data,chains=design$chains,
    parallel_chains=design$chains,seed=design$seed,
    iter_warmup=design$warmup,iter_sampling=design$sampling,
    adapt_delta=design$adapt_delta,
    max_treedepth=design$max_treedepth,refresh=250
  )
  fit$save_object(file=fit_file)
  saveRDS(design,design_file)
}

# ---- 4. Validate MCMC convergence before writing the shared CSV ------
parameters <- c("O18_R","O18_LSR","Cl_R","Cl_LSR",
                "f_LSR","sigma_O18","tau_LSR_O18","sigma_Cl")
diagnostics <- fit$summary(variables=parameters)
chains <- fit$diagnostic_summary()
if (!nrow(diagnostics) ||
    any(!is.finite(diagnostics$rhat)) ||
    any(!is.finite(diagnostics$ess_bulk)) ||
    any(!is.finite(diagnostics$ess_tail)) ||
    any(!is.finite(chains$ebfmi)))
  stop("Missing/nonfinite MCMC diagnostics; cannot export fractions.")

rhat_max <- max(diagnostics$rhat)
ess_bulk_min <- min(diagnostics$ess_bulk)
ess_tail_min <- min(diagnostics$ess_tail)
ebfmi_min <- min(chains$ebfmi)
divergences <- sum(chains$num_divergent)
tree_hits <- sum(chains$num_max_treedepth)
passed <- rhat_max<=1.05 && ess_bulk_min>=400 &&
  ess_tail_min>=400 && ebfmi_min>=0.3 &&
  divergences==0 && tree_hits==0

qc <- data.frame(
  metric=c("max_rhat","min_bulk_ESS","min_tail_ESS","min_EBFMI",
           "divergences","treedepth_hits","passed","input_rows",
           "fitted_rows","excluded_rows","one_tracer_rows"),
  value=as.character(c(rhat_max,ess_bulk_min,ess_tail_min,
                       ebfmi_min,divergences,tree_hits,passed,
                       nrow(dat),length(selected),sum(excluded),
                       sum(xor(has_O18,has_Cl) & !excluded)))
)
readr::write_csv(qc,diagnostic_file)
if (!passed)
  stop("Tracer-only fit failed MCMC diagnostic thresholds. See: ",
       diagnostic_file,". No shared CSV was written.")

# ---- 5. Summarise fractions; preserve all original input columns -----
variables <- paste0("f_LSR[",seq_along(selected),"]")
draws <- posterior::as_draws_matrix(fit$draws(variables=variables))
if (!all(variables %in% colnames(draws)))
  stop("Missing per-sample f_LSR draws in the saved tracer-only fit.")
draws <- draws[,variables,drop=FALSE]

lsr_med <- apply(draws,2,stats::median)
lsr_low <- apply(draws,2,stats::quantile,probs=0.025)
lsr_high <- apply(draws,2,stats::quantile,probs=0.975)

# Complement river and LSR FRACTIONS within each posterior draw.
# The 95% river limits are the inverse of the 95% LSR limits.
result <- dat
for (nm in new_columns[1:6])
  result[[nm]] <- rep(NA_real_,nrow(dat))
result$LSR_fraction_median[selected] <- lsr_med
result$LSR_fraction_lower95[selected] <- lsr_low
result$LSR_fraction_upper95[selected] <- lsr_high
result$river_fraction_median[selected] <- 1-lsr_med
result$river_fraction_lower95[selected] <- 1-lsr_high
result$river_fraction_upper95[selected] <- 1-lsr_low

result$mixing_tracers_used <- ifelse(
  has_O18 & has_Cl,"O18 + Cl",
  ifelse(has_O18,"O18 only",ifelse(has_Cl,"Cl only","none"))
)
result$mixing_status <- ifelse(
  excluded,"excluded (16758)",
  ifelse(has_O18 & has_Cl,"two tracers",
         ifelse(has_O18 | has_Cl,"one tracer - caution",
                "no available tracer"))
)

# Preserve original rows and input text values; append only new columns.
stopifnot(nrow(result)==nrow(dat),
          identical(result$Site_ID,dat$Site_ID),
          all(vapply(names(dat),function(nm)
            identical(result[[nm]],dat[[nm]]),logical(1))),
          all(is.na(result$LSR_fraction_median[excluded])))
# Never overwrite the SOURCE file. Only the derived COPY is replaced.
readr::write_csv(result,output_file,na="NA")

writeLines(c(
  "Tracer-only Bayesian two-endmember WATER mixing fractions.",
  paste("Input:",input_file," MD5:",design$input_md5),
  paste("Stan model:",stan_file," MD5:",design$stan_md5),
  paste("Fitted rows:",length(selected),"of",nrow(dat)),
  "Likelihood uses O18 and Cl only. No DO, MRT, N2 or source labels.",
  "Endmember priors, bounds and tracer errors match the primary",
  "DO model's tracer component (endmember_prior_scale=1).",
  "LSR and river fractions are proportions of WATER, not of DO.",
  "Intervals are marginal equal-tailed 95% Bayesian credible intervals.",
  "Site 16758 excluded; duplicate Site_ID rows treated separately.",
  "Samples with only one tracer have less constrained mixing fractions.",
  "The original input CSV remains unchanged, and all original rows are kept.",
  "For downstream uncertainty propagation use draws in tracer_only_fit.rds.",
  paste("Convergence: max Rhat",signif(rhat_max,4),
        "min bulk ESS",signif(ess_bulk_min,4),
        "min tail ESS",signif(ess_tail_min,4))
),method_file)

message("CREATED: ",normalizePath(output_file))
message("Fractions appended for ",length(selected)," of ",nrow(dat)," rows.")
message("Unmodified input: ",input_file)
message("Sampler diagnostics: ",diagnostic_file)
