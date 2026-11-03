# 05_DO_kinetic_trajectories and summary.R
#
# FOUR complementary DO vs MRT plots for the preferred model:
#  1. Observed mixed-groundwater DO against its sample-specific fitted DO
#     (mean_DO), with 95% posterior predictive intervals (replicated_DO).
#  2. RAW measured mixed-groundwater DO against the posterior river-only
#     100%-river reference (no sample-level fitted DO in this panel).
#  3. Sample-specific DO attributable to river recharge,
#     (1 - f_LSR) * river_DO_after_aquifer, alongside the same
#     100%-river reference trajectory.
#  4. STACKED fitted river and LSR DO contributions, overlaid with
#     measured mixed DO and a credible interval for the fitted total.
#
# Plot 4 uses POSTERIOR MEANS for its stacked source contributions and
# fitted total, so the components sum exactly; the first plot continues
# to use posterior medians for the fitted point values.
#
# The river-only trajectory is NOT the fitted mean of the mixed samples.
# No Stan sampling or modification of fitted objects occurs.
#
# Preferred model: 39 fitted DO observations, river DO
#   ~ Normal(7.5, 1.5^2) truncated [0, 14.65] mg/L;
# first-order aquifer DO kinetics; 2801 and 2130 DO excluded from
# the likelihood but their delta18O and chloride retained.
#
# Run from the repository root:
#   source("scripts/05_DO_kinetic_trajectories and summary.R")

if (!file.exists("scripts/00_setup.R"))
  stop("Set working directory to the repository root.")
source("scripts/00_setup.R", local=TRUE)

packages05 <- c("ggplot2","posterior","viridisLite")
not_found05 <- packages05[
  !vapply(packages05,requireNamespace,logical(1),quietly=TRUE)
]
if (length(not_found05))
  stop("Install R package(s): ",paste(not_found05,collapse=", "))

# ---- Configuration -----------------------------------------------------
suite05 <- "2801_masked_no2130"
scenario05 <- "C_prior_sd1p5"
order05 <- 1L
fit_dir05 <- file.path("outputs","04_river_DO_prior",suite05)
fit_path05 <- file.path(fit_dir05,
  paste0("fit_",scenario05,"_order",order05,".rds"))
status_path05 <- file.path(fit_dir05,"fit_status.csv")
plot_dir05 <- file.path("plots","final_39well")
dir.create(plot_dir05,recursive=TRUE,showWarnings=FALSE)
turbo05 <- viridisLite::turbo(256)
river_blue05 <- turbo05[32L]

if (!all(file.exists(c(fit_path05,status_path05))))
  stop("Preferred saved fit or fit_status.csv missing in ",fit_dir05,
       ". Run the full-data fitting script only if outputs are absent.")

fit_status05 <- utils::read.csv(status_path05,stringsAsFactors=FALSE)
needed05 <- c("scenario","aquifer_order","all_diagnostics_pass")
if (!all(needed05 %in% names(fit_status05)))
  stop("fit_status.csv missing required fields.")
passed05 <- fit_status05$all_diagnostics_pass[
  fit_status05$scenario==scenario05 &
    fit_status05$aquifer_order==order05
]
if (length(passed05)!=1L || !isTRUE(passed05))
  stop("Preferred model diagnostics are absent or failed; do not plot it.")

# ---- Obtain well-level fitted quantities from the same posterior ------
stopifnot(nrow(d)==41L)
N05 <- nrow(d)
idx05 <- seq_len(N05)
key05 <- c("river_DO","aquifer_rate_actual",
           "LSR_DO_at_water_table","LSR_DO_deficit")
v05 <- function(name) paste0(name,"[",idx05,"]")
vars05 <- c(key05, "k_ref",
            v05("f_LSR"),v05("mean_DO"),
            v05("replicated_DO"),v05("river_DO_after_aquifer"))

fit05 <- readRDS(fit_path05)
z05 <- posterior::as_draws_matrix(fit05$draws(variables=vars05))
if (!all(vars05 %in% colnames(z05)))
  stop("Saved model does not contain the expected well-level draws.")
if (any(!is.finite(z05[,vars05,drop=FALSE])))
  stop("Saved posterior contains non-finite draws.")

mat05 <- function(name) z05[,v05(name),drop=FALSE]
f05 <- mat05("f_LSR")
mu05 <- mat05("mean_DO")
rep05 <- mat05("replicated_DO")
river05 <- mat05("river_DO_after_aquifer")
river_contrib05 <- (1-f05)*river05   # mg/L contributed to each MIXED sample

qcol05 <- function(mat,p) apply(mat,2,stats::quantile,probs=p)

wells05 <- data.frame(
  Site_ID=d$Site_ID,
  MRT=d$MRT,
  DO_observed=d$DO_original,
  excluded_DO=d$Site_ID %in% c(cfg$river_proxy_site,"2130"),
  f_LSR=apply(f05,2,stats::median),
  DO_fitted=apply(mu05,2,stats::median),
  DO_pred_lower95=qcol05(rep05,0.025),
  DO_pred_upper95=qcol05(rep05,0.975),
  river_contrib=apply(river_contrib05,2,stats::median),
  river_contrib_lower95=qcol05(river_contrib05,0.025),
  river_contrib_upper95=qcol05(river_contrib05,0.975),
  stringsAsFactors=FALSE
)
if (sum(!wells05$excluded_DO)!=39L ||
    anyNA(wells05$DO_fitted))
  stop("Unexpected sample mask or missing predictions.")

# ---- Posterior river-only reference curve -----------------------------
# The curve uses the plotted (reported) MRT as a fixed horizontal value.
# Sample-specific contributions above use the model's latent true_MRT.
mrt_grid05 <- seq(0,max(d$MRT,na.rm=TRUE)*1.05,length.out=251L)
start05 <- z05[,"river_DO"]
rate05 <- z05[,"k_ref"]
river_curve05 <- vapply(
  mrt_grid05,
  function(t) start05*exp(-rate05*t),
  numeric(nrow(z05))
)
curve05 <- data.frame(
  MRT=mrt_grid05,
  median=apply(river_curve05,2,stats::median),
  lower95=qcol05(river_curve05,0.025),
  upper95=qcol05(river_curve05,0.975)
)

# ---- Shared appearance -------------------------------------------------
gradient05 <- ggplot2::scale_colour_gradientn(
  colours=turbo05[32:230],limits=c(0,1),
  name="LSR fraction",
  breaks=c(0,.25,.5,.75,1),
  labels=c("0%","25%","50%","75%","100%")
)
theme05 <- ggplot2::theme_bw(base_size=13) +
  ggplot2::theme(
    panel.grid.minor=ggplot2::element_blank(),
    plot.title=ggplot2::element_text(face="bold"),
    legend.position="right",
    legend.box="vertical"
  )
save05 <- function(plot,file,width=10,height=7) {
  path <- file.path(plot_dir05,file)
  print(plot)
  ggplot2::ggsave(path,plot=plot,width=width,height=height,
                  dpi=300,units="in",bg="white")
  message("PLOT: ",normalizePath(path))
}

excluded05 <- wells05[wells05$excluded_DO,,drop=FALSE]

# 1. Goodness-of-fit: compare OBSERVED mixed DO with fitted MIXED DO.
# 95% ranges are POSTERIOR PREDICTIVE for a sample observation,
# not just uncertainty on the conditional mean or a river reference.
p_fit05 <- ggplot2::ggplot(wells05,ggplot2::aes(x=MRT)) +
  ggplot2::geom_segment(
    ggplot2::aes(xend=MRT,y=DO_observed,yend=DO_fitted),
    colour="grey70",linewidth=.55
  ) +
  ggplot2::geom_linerange(
    ggplot2::aes(ymin=DO_pred_lower95,
                 ymax=DO_pred_upper95,colour=f_LSR),
    linewidth=.8,alpha=.60
  ) +
  ggplot2::geom_point(
    ggplot2::aes(y=DO_observed,shape=excluded_DO),
    colour="grey15",size=3,stroke=.85
  ) +
  ggplot2::geom_point(
    ggplot2::aes(y=DO_fitted,colour=f_LSR),
    size=2.5,alpha=.95
  ) +
  ggplot2::geom_text(
    data=excluded05,
    ggplot2::aes(y=DO_observed,label=Site_ID),
    colour="grey20",vjust=-.9,size=3.4,show.legend=FALSE
  ) +
  gradient05 +
  ggplot2::scale_shape_manual(
    values=c("FALSE"=1,"TRUE"=4),
    labels=c("FALSE"="DO used in fit","TRUE"="DO excluded from fit"),
    name="Observed DO"
  ) +
  ggplot2::scale_y_continuous(
    expand=ggplot2::expansion(mult=c(.03,.06))
  ) +
  ggplot2::coord_cartesian(ylim=c(0,11)) +
  ggplot2::labs(
    title="Observed and model-fitted groundwater DO",
    subtitle="Open circles / crosses: measured DO; coloured points: fitted mixed-groundwater DO",
    x="Reported groundwater MRT (years)",
    y="Dissolved oxygen (mg/L)",
    caption="Vertical ranges: 95% in-sample posterior predictive intervals; 2801 and 2130 DO excluded from fitting."
  ) + theme05
save05(p_fit05,
       paste0("13_DO_observed_vs_fitted_mixed_",suite05,".png"),
       width=11,height=7.5)

# 2. Observed mixed samples against the 100%-river-only reference.
# This is intentionally a PROCESS REFERENCE, not a model-observation
# goodness-of-fit plot. All observed points include their LSR fraction.
p_raw05 <- ggplot2::ggplot() +
  ggplot2::geom_ribbon(
    data=curve05,
    ggplot2::aes(x=MRT,ymin=lower95,ymax=upper95),
    fill=river_blue05,alpha=.16
  ) +
  ggplot2::geom_line(
    data=curve05,
    ggplot2::aes(x=MRT,y=median),
    colour="grey20",linewidth=.9,linetype="longdash"
  ) +
  ggplot2::geom_point(
    data=wells05,
    ggplot2::aes(x=MRT,y=DO_observed,
                 colour=f_LSR,shape=excluded_DO),
    size=2.8,alpha=.95
  ) +
  ggplot2::geom_text(
    data=excluded05,
    ggplot2::aes(x=MRT,y=DO_observed,label=Site_ID),
    colour="grey20",vjust=-.9,size=3.4,show.legend=FALSE
  ) +
  gradient05 +
  ggplot2::scale_shape_manual(
    values=c("FALSE"=16,"TRUE"=1),
    labels=c("FALSE"="DO used in fit","TRUE"="DO excluded from fit"),
    name="Observed DO"
  ) +
  ggplot2::scale_y_continuous(
    limits=c(0,NA),expand=ggplot2::expansion(mult=c(0,.05))
  ) +
  ggplot2::labs(
    title="First-order river reference with uncertain river DO",
    subtitle="Raw groundwater DO (coloured by LSR fraction) versus 100% river-only first-order reference",
    x="Reported groundwater MRT (years)",
    y="Dissolved oxygen (mg/L)",
    caption="The 95% blue band describes river-only model uncertainty, NOT a predictive interval for mixed groundwater samples."
  ) + theme05
save05(p_raw05,
       paste0("14_DO_river_reference_raw_samples_",suite05,".png"),
       width=11,height=7.5)

# 3. River-derived contribution to each mixed groundwater sample.
# Posterior contributions use JOINT draws of source fraction and river
# concentration AFTER aquifer reaction, not products of separate medians.
# The dashed line is the 100%-river REFERENCE, never the expected mixed DO.
p_river05 <- ggplot2::ggplot() +
  ggplot2::geom_ribbon(
    data=curve05,
    ggplot2::aes(x=MRT,ymin=lower95,ymax=upper95),
    fill=river_blue05,alpha=.16
  ) +
  ggplot2::geom_line(
    data=curve05,
    ggplot2::aes(x=MRT,y=median),
    colour="grey20",linewidth=.9,linetype="longdash"
  ) +
  ggplot2::geom_linerange(
    data=wells05,
    ggplot2::aes(x=MRT,ymin=river_contrib_lower95,
                 ymax=river_contrib_upper95,colour=f_LSR),
    alpha=.70,linewidth=.8
  ) +
  ggplot2::geom_point(
    data=wells05,
    ggplot2::aes(x=MRT,y=river_contrib,colour=f_LSR),
    size=2.8,alpha=.95
  ) +
  ggplot2::geom_text(
    data=excluded05,
    ggplot2::aes(x=MRT,y=river_contrib,label=Site_ID),
    colour="grey20",vjust=-.9,size=3.4,show.legend=FALSE
  ) +
  gradient05 +
  ggplot2::scale_y_continuous(
    limits=c(0,NA),
    expand=ggplot2::expansion(mult=c(0,.05))
  ) +
  ggplot2::labs(
    title="River-derived DO contribution in each groundwater sample",
    subtitle="Dashed curve: 100% river-recharge reference; points: modelled river contribution after mixing",
    x="Reported groundwater MRT (years)",
    y="River-attributable DO in the mixed sample (mg/L)",
    caption="Sample points = (1 - f_LSR) x river DO after aquifer reaction; vertical ranges and reference band are 95% posterior credible intervals."
  ) + theme05
save05(p_river05,
       paste0("15_DO_river_contribution_",suite05,".png"),
       width=11,height=7.5)

# 4. Fitted source decomposition: how river and LSR add to total DO.
# Plot the posterior MEAN of each contribution, because posterior
# means add exactly: E[river] + E[LSR] = E[total fitted DO].
# A sum of separate posterior medians would not necessarily equal
# the posterior median of the total, so do not stack those medians.
lsr_contrib05 <- mu05 - river_contrib05

wells05$river_component_mean <- colMeans(river_contrib05)
wells05$lsr_component_mean <- colMeans(lsr_contrib05)
wells05$fitted_total_mean <- colMeans(mu05)
wells05$fitted_lower95 <- qcol05(mu05,0.025)
wells05$fitted_upper95 <- qcol05(mu05,0.975)

stopifnot(
  max(abs(wells05$river_component_mean +
          wells05$lsr_component_mean -
          wells05$fitted_total_mean)) < 1e-8
)

# Stacked vertical segments: blue from zero to river contribution,
# red from the river contribution up to the full fitted concentration.
components05 <- rbind(
  data.frame(
    MRT=wells05$MRT,
    bottom=0,
    top=wells05$river_component_mean,
    source="River"
  ),
  data.frame(
    MRT=wells05$MRT,
    bottom=wells05$river_component_mean,
    top=wells05$fitted_total_mean,
    source="LSR"
  )
)

p_components05 <- ggplot2::ggplot() +
  ggplot2::geom_linerange(
    data=wells05,
    ggplot2::aes(x=MRT,ymin=fitted_lower95,ymax=fitted_upper95),
    colour="grey60",linewidth=.75,alpha=.80
  ) +
  ggplot2::geom_segment(
    data=components05,
    ggplot2::aes(x=MRT,xend=MRT,y=bottom,yend=top,colour=source),
    linewidth=2.3,alpha=.90
  ) +
  ggplot2::geom_point(
    data=wells05,
    ggplot2::aes(x=MRT,y=DO_observed,shape=excluded_DO),
    colour="grey15",size=3.2,stroke=.9
  ) +
  ggplot2::geom_point(
    data=wells05,
    ggplot2::aes(x=MRT,y=fitted_total_mean),
    colour="black",size=2.2
  ) +
  ggplot2::geom_text(
    data=excluded05,
    ggplot2::aes(x=MRT,y=DO_observed,label=Site_ID),
    colour="grey20",vjust=-.9,size=3.4,show.legend=FALSE
  ) +
  ggplot2::scale_colour_manual(
    values=c("River"=river_blue05,"LSR"=turbo05[230L]),
    breaks=c("River","LSR"),
    name="Fitted DO contribution"
  ) +
  ggplot2::scale_shape_manual(
    values=c("FALSE"=1,"TRUE"=4),
    labels=c("FALSE"="DO used in fit","TRUE"="DO excluded from fit"),
    name="Observed DO"
  ) +
  ggplot2::coord_cartesian(ylim=c(0,11)) +
  ggplot2::labs(
    title="River and LSR contributions to fitted groundwater DO",
    subtitle="Blue + red = fitted total (black dot); open circle / cross = observed DO",
    x="Reported groundwater MRT (years)",
    y="Dissolved oxygen (mg/L)",
    caption=paste(
      "Stacked contributions and black totals are posterior means;",
      "grey whiskers are 95% credible intervals for fitted mean DO.",
      "2801 and 2130 DO were excluded from calibration."
    )
  ) + theme05

save05(p_components05,
       paste0("16_DO_fitted_components_",suite05,".png"),
       width=11,height=7.5)

# ---- Four-parameter summary from the SAME preferred posterior --------
labels05 <- c("River-entry DO","First-order aquifer rate",
              "LSR DO at water table","LSR DO deficit")
units05 <- c("mg/L","yr\u207b\u00b9","mg/L","mg/L")
digits05 <- c(2L,4L,2L,2L)
summary05 <- data.frame(
  parameter=labels05,variable=key05,
  median=vapply(key05,function(v) stats::median(z05[,v]),numeric(1)),
  lower95=vapply(key05,function(v)
    unname(stats::quantile(z05[,v],.025)),numeric(1)),
  upper95=vapply(key05,function(v)
    unname(stats::quantile(z05[,v],.975)),numeric(1))
)
fmt05 <- function(x) paste0(
  "| ",paste(x,collapse=" | ")," |"
)
summary_table05 <- c(
  "| Parameter | Posterior median | 95% credible interval |",
  "|---|---:|---:|",
  vapply(seq_along(key05),function(i) fmt05(c(
    labels05[i],
    paste(sprintf("%.*f",digits05[i],summary05$median[i]),units05[i]),
    paste0(
      sprintf("%.*f",digits05[i],summary05$lower95[i]),"\u2013",
      sprintf("%.*f",digits05[i],summary05$upper95[i])
    )
  )),character(1))
)
summary_path05 <- file.path(
  plot_dir05,paste0("17_DO_primary_posterior_summary_",suite05,".md")
)
writeLines(enc2utf8(c(
  "# Preferred first-order DO model: posterior summary",
  "",
  "River DO ~ Normal(7.5, 1.5\u00b2) mg/L, truncated to [0, 14.65].",
  paste0("Saved fit: \u0060",fit_path05,"\u0060"),
  "",
  summary_table05,
  "",
  "The LSR DO deficit is a concentration difference, not a rate.",
  "The DO vs MRT figures are in-sample posterior diagnostics, not held-out CV.",
  "Panel 2 plots raw sample DO against the river-only reference; its band is not a mixed-sample predictive interval.",
  "Panel 4 uses additive POSTERIOR MEAN source contributions; its whiskers describe uncertainty in fitted mean DO, not replicated observations."
)),summary_path05,useBytes=TRUE)
print(summary05,row.names=FALSE,digits=5)
message("SUMMARY: ",normalizePath(summary_path05))
