# 03_high_DO_audit.R — investigate high dissolved oxygen wells
#
# Standalone audit: NO MCMC sampling, NO model compilation, no edits to data.
# Reuses saved 01 posterior draws and the exact 00 setup / boundaries.
# Run from repository root:
#   source("scripts/03_high_DO_audit.R")
#
# Tests:
# (1) Which wells exceed fixed river-entry DO of 7.5 mg/L?
# (2) Is their measured DO plausible relative to saturation at sample T?
# (3) Does site 2130 exceed a conservative physical DO upper bound at
#     its reported MRT and the fitted rate? This bound assumes 100% of
#     whichever entry source is most oxygenated, with ZERO LSR DO loss.
# (4) Does the FULL-DATA posterior reproduce 2130's observed DO?
#     This is an IN-SAMPLE diagnostic, NOT held-out validation.
#
# The proxy site 2801 is an observational reference ONLY. Its DO is not
# used in 01's DO likelihood, and MUST NOT be scored as an ordinary well.

if (!file.exists("scripts/00_setup.R"))
  stop("Set working directory to Oxygen-Reduction-Joint-Model root.")
source("scripts/00_setup.R",local=TRUE)

if (!"Temp" %in% names(d))
  stop("Private input is missing the measured groundwater Temp column.")
if (!requireNamespace("posterior",quietly=TRUE) ||
    !requireNamespace("viridisLite",quietly=TRUE))
  stop("Install the posterior and viridisLite packages to run this audit.")

dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
dir.create(plotdir,recursive=TRUE,showWarnings=FALSE)

high <- d %>%
  mutate(
    DO_measured=DO_original,
    DO_saturation_at_sample_T=if_else(
      is.finite(Temp)&Temp>=0&Temp<=30,
      do_sat_mgL(Temp,cfg$atmospheric_pressure_ratio),
      NA_real_
    ),
    fraction_saturation=DO_measured/DO_saturation_at_sample_T,
    percent_saturation=100*fraction_saturation,
    status=if_else(
      Site_ID==cfg$river_proxy_site,
      "river proxy: DO masked in Stan",
      "DO fitted in Stan"
    )
  ) %>%
  filter(DO_measured>cfg$river_DO_boundary) %>%
  arrange(desc(DO_measured)) %>%
  select(
    Site_ID,site,DO_measured,Temp,DO_saturation_at_sample_T,
    percent_saturation,MRT,O18,Cl,Na,source,status,
    sample_index
  )
readr::write_csv(high,output("high_DO_measurement_audit.csv"))
message("\nORIGINAL HIGH DO WELLS (threshold=",
        cfg$river_DO_boundary," mg/L):")
print(high,n=Inf,width=Inf)
if (!all(c("2130",cfg$river_proxy_site) %in% high$Site_ID))
  warning("Expected high DO records 2130 and 2801 not both found. ",
          "Inspect local private input version.")

# Summarise the physical limit for 2130. A common aquifer reaction rate
# and positive source losses can only LOWER DO below this maximum.
site_id <- "2130"
ix <- which(d$Site_ID==site_id)
if (length(ix)!=1L) stop("Expected exactly one site 2130 record.")
site <- d[ix,,drop=FALSE]
t_reported <- as.numeric(site$MRT[1])
obs <- as.numeric(site$DO_original[1])
Cmax_at_entry <- max(Csat,cfg$river_DO_boundary)

global_file <- output("global_posteriors.csv")
if (!file.exists(global_file))
  stop("Missing ",global_file,
       ". Run 01 first. No MCMC is started by this audit.")

g <- readr::read_csv(global_file,show_col_types=FALSE)
need_vars <- c("aquifer_rate_actual",
               "LSR_DO_at_water_table","mean_f_LSR")
if (!all(need_vars %in% g$variable))
  stop("Missing expected parameters in ",global_file)
parameters <- g %>%
  filter(variable %in% need_vars) %>%
  select(aquifer_order,variable,median) %>%
  tidyr::pivot_wider(names_from=variable,values_from=median)
if (!setequal(parameters$aquifer_order,0:1))
  stop("Need complete zero- and first-order posterior summaries.")

limit <- parameters %>%
  mutate(
    Site_ID=site_id,
    observed_DO=obs,
    reported_MRT=t_reported,
    Csat_assumed_recharge=Csat,
    river_DO_entry=cfg$river_DO_boundary,
    Cmax_entry_no_LSR_loss=Cmax_at_entry,
    MAX_DO_at_reported_MRT=if_else(
      aquifer_order==0L,
      pmax(0,Cmax_at_entry-aquifer_rate_actual*t_reported),
      Cmax_at_entry*exp(-aquifer_rate_actual*t_reported)
    ),
    DO_above_physical_bound=obs-MAX_DO_at_reported_MRT,
    highest_possible_MRT_for_observed_DO=if_else(
      aquifer_order==0L,
      pmax(0,(Cmax_at_entry-obs)/aquifer_rate_actual),
      pmax(0,log(Cmax_at_entry/obs)/aquifer_rate_actual)
    ),
    interpretation=paste(
      "Optimistic bound: 100% most oxygen-rich recharge,",
      "NO pre-water-table loss, at reported MRT; uses median aquifer rate"
    )
  )
readr::write_csv(limit,output("high_DO_physical_limits.csv"))
message("\nSITE 2130 — OPTIMISTIC PHYSICAL DO MAXIMA (mg/L):")
print(
  limit %>% select(
    aquifer_order,observed_DO,reported_MRT,
    aquifer_rate_actual,MAX_DO_at_reported_MRT,
    DO_above_physical_bound,
    highest_possible_MRT_for_observed_DO
  ),n=Inf,width=Inf
)

# Indicative TWO-TRACER mixing discrepancy at site 2130.
# These are point calculations from endmember PRIOR MEANS, NOT
# posterior mixing fractions or a third-source inference.
fmix <- tibble::tibble(
  Site_ID=site_id,
  method=c("delta18O prior mean","Cl prior mean"),
  approx_f_LSR=c(
    (site$O18[1]-(-7.79))/(-6.58-(-7.79)),
    (site$Cl[1]-4.9)/(60-4.9)
  ),
  note="Not a posterior estimate; assumes original tracer endmember prior means"
)
readr::write_csv(fmix,output("high_DO_tracer_fraction_check.csv"))
message("\nSITE 2130 — INDICATIVE SOURCE FRACTION BY EACH TRACER:")
print(fmix,n=Inf,width=Inf)

# Saved CmdStanR fit objects permit in-sample posterior checks WITHOUT
# another MCMC run. The 2130 observation was used to estimate these
# posteriors, so posterior predictive tail probabilities here are
# descriptive, not independent validation. True held-out testing requires
# refitting with has_DO[site2130] = 0.
pred <- list()
for (order in 0:1) {
  file <- output(paste0("fit_order",order,".rds"))
  if (!file.exists(file)) {
    warning("Saved fit missing: ",file,
            ". Skipping posterior check for order ",order)
    next
  }
  fit <- readRDS(file)
  vars <- paste0(
    c("mean_DO","replicated_DO","true_MRT","f_LSR"),"[",ix,"]"
  )
  draws <- posterior::as_draws_matrix(fit$draws(
    variables=c(vars,"aquifer_rate_actual")
  ))
  to_vec <- function(v) as.numeric(draws[,v])
  mu <- to_vec(vars[1])
  repdo <- to_vec(vars[2])
  age <- to_vec(vars[3])
  f <- to_vec(vars[4])
  rate <- to_vec("aquifer_rate_actual")
  # Draw-specific upper bound: use actual latent age and rate, with
  # 100% of the most oxygen-rich source and no LSR loss.
  upper <- if (order==0)
    pmax(0,Cmax_at_entry-rate*age)
  else Cmax_at_entry*exp(-rate*age)

  pred[[length(pred)+1L]] <- tibble::tibble(
    Site_ID=site_id,
    aquifer_order=order,
    observed_DO=obs,
    n_draws=length(mu),
    median_f_LSR=median(f),
    median_true_MRT=median(age),
    lower95_true_MRT=unname(quantile(age,.025)),
    upper95_true_MRT=unname(quantile(age,.975)),
    fitted_DO_median=median(mu),
    fitted_DO_lower95=unname(quantile(mu,.025)),
    fitted_DO_upper95=unname(quantile(mu,.975)),
    replicated_DO_median=median(repdo),
    replicated_DO_lower95=unname(quantile(repdo,.025)),
    replicated_DO_upper95=unname(quantile(repdo,.975)),
    P_fitted_mean_ge_measured=mean(mu>=obs),
    P_replicated_DO_ge_measured=mean(repdo>=obs),
    P_optimistic_bound_ge_measured=mean(upper>=obs)
  )
}
if (length(pred)) {
  chk <- dplyr::bind_rows(pred)
  readr::write_csv(chk,output("high_DO_2130_posterior_check.csv"))
  message("\nSITE 2130 — IN-SAMPLE POSTERIOR CHECK:")
  print(chk,n=Inf,width=Inf)

  # Turbo reference colours: blue=river, red=rainfall/LSR;
  # green for fitted/predictive MODEL statistics, not for a source.
  turb <- viridisLite::turbo(256)
  chk_plot <- chk %>%
    tidyr::pivot_longer(
      cols=c(
        fitted_DO_median,replicated_DO_median
      ),
      names_to="quantity",values_to="DO_median"
    ) %>%
    mutate(
      DO_lower95=dplyr::if_else(
        quantity=="fitted_DO_median",
        fitted_DO_lower95,replicated_DO_lower95
      ),
      DO_upper95=dplyr::if_else(
        quantity=="fitted_DO_median",
        fitted_DO_upper95,replicated_DO_upper95
      ),
      quantity=dplyr::recode(
        quantity,
        fitted_DO_median="Fitted mean (95% credible)",
        replicated_DO_median="Replicated DO (95% predictive)"
      ),
      aquifer_order=factor(
        aquifer_order,levels=0:1,
        labels=c("Zero order","First order")
      )
    )
  p <- ggplot2::ggplot(
    chk_plot,ggplot2::aes(x=aquifer_order,y=DO_median,
                         colour=quantity)
  )+
    ggplot2::geom_hline(
      yintercept=obs,linetype="dashed",
      colour="grey40",linewidth=.7
    )+
    ggplot2::geom_linerange(
      ggplot2::aes(ymin=DO_lower95,ymax=DO_upper95),
      position=ggplot2::position_dodge(width=.45),
      linewidth=.9
    )+
    ggplot2::geom_point(
      position=ggplot2::position_dodge(width=.45),size=2.9
    )+
    ggplot2::scale_colour_manual(
      values=c(
        "Fitted mean (95% credible)"=turb[100],
        "Replicated DO (95% predictive)"=turb[180]
      ),
      name=NULL
    )+
    ggplot2::labs(
      title=paste0("Site ",site_id,": high-DO model check"),
      subtitle=paste("Observed DO =",obs,
                     "mg/L (dashed); full-data, IN-SAMPLE posteriors"),
      x="Aquifer kinetic model",y="Dissolved oxygen (mg/L)",
      caption=paste(
        "Proxy 2801 DO was masked. Site 2130 DO was fitted,",
        "so this is not a held-out validation."
      )
    )+
    ggplot2::theme_bw(base_size=12)+
    ggplot2::theme(legend.position="bottom")
  print(p)
  path <- file.path(plotdir,"09_high_DO_site2130_audit.png")
  ggplot2::ggsave(
    path,p,width=9,height=5.5,dpi=300,bg="white"
  )
  message("Figure: ",normalizePath(path))
}
message("\nNo sampling or refitting performed by script 03.")
