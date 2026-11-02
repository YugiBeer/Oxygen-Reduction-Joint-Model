# 02_DO_plotting.R — FINAL 39-well plots only
#
# NO CmdStan compilation, saved fit loading, or MCMC.
# Reads final 04 full-data summaries and 05 held-out CV predictions.
# Original script-01 diagnostic plots preserved in:
#   archive/scripts/02_DO_plotting_legacy.R
#
# Run from repo root:
#   source("scripts/02_DO_plotting.R")
#
# The separate script 06 plots actual posterior samples and correlations.

if (!file.exists("scripts/00_setup.R"))
  stop("Set the working directory to repository root.")
source("scripts/00_setup.R", local=TRUE)

if (!requireNamespace("viridisLite", quietly=TRUE))
  stop("Install the viridisLite package to use Turbo colours.")
if (!requireNamespace("ggplot2", quietly=TRUE))
  stop("Install ggplot2 to make plots.")

plotdir <- file.path("plots","final_39well")
dir.create(plotdir,recursive=TRUE,showWarnings=FALSE)
manifest <- character()

turbo <- viridisLite::turbo(256L)
river_blue <- turbo[32L]
rain_red <- turbo[230L]
aquifer_green <- turbo[128L]
turbo_gradient <- turbo[32:230]
ref_colour <- "grey55"

fraction_colour <- function() {
  ggplot2::scale_colour_gradientn(
    colours=turbo_gradient,limits=c(0,1),
    name="LSR fraction",breaks=c(0,.25,.5,.75,1)
  )
}
order_label <- function(x) {
  factor(as.integer(x),levels=0:1,
         labels=c("Zero-order aquifer","First-order aquifer"))
}
base_theme <- function() {
  ggplot2::theme_bw(base_size=12)+
    ggplot2::theme(
      panel.grid.minor=ggplot2::element_blank(),
      legend.position="bottom",
      strip.text=ggplot2::element_text(face="bold"),
      plot.title=ggplot2::element_text(face="bold"),
      plot.subtitle=ggplot2::element_text(size=10),
      plot.caption=ggplot2::element_text(hjust=0)
    )
}
require_cols <- function(df, columns, filename) {
  missing <- setdiff(columns,names(df))
  if (length(missing))
    stop("Unexpected format in ",filename,": missing ",
         paste(missing,collapse=", "))
}
save_plot <- function(p,filename,width=10,height=6) {
  print(p)
  destination <- file.path(plotdir,filename)
  ggplot2::ggsave(
    filename=destination,plot=p,width=width,height=height,
    units="in",dpi=300,bg="white"
  )
  manifest <<- c(manifest,destination)
  message("PLOT: ",normalizePath(destination))
  invisible(destination)
}

# Both final analyses must be present; no legacy baseline output required.
needed02 <- c(
  file.path("outputs","04_river_DO_prior",
            "2801_masked_no2130","comparison_summary.csv"),
  file.path("outputs","05_river_DO_grouped_CV",
            "2801_masked_no2130","cv_predictions.csv"),
  file.path("outputs","05_river_DO_grouped_CV",
            "2801_masked_no2130","cv_scores.csv")
)
if (any(!file.exists(needed02))) {
  stop("Missing final 39-well output(s):\n",
       paste(needed02[!file.exists(needed02)],collapse="\n"),
       "\nRun scripts 04/05 (or copy their saved results) first.")
}

# 01–02: Final 39-well river DO prior sensitivity (posterior summaries).
# River and rainfall source colours use the same Turbo palette as script 06.
for (proxy_suite in "2801_masked_no2130") {
  sensitivity_path <- file.path(
    "outputs","04_river_DO_prior",proxy_suite,"comparison_summary.csv"
  )
  if (!file.exists(sensitivity_path)) next
  sens <- readr::read_csv(
    sensitivity_path,show_col_types=FALSE
  ) %>%
    mutate(
      order_label=order_label(aquifer_order),
      scenario=factor(
        scenario,levels=c(
          "A_fixed_7p5","B_prior_sd1","C_prior_sd1p5"
        )
      )
    )
  if (anyNA(sens$scenario)) {
    warning("Unexpected scenario names in ",sensitivity_path,
            "; skipping sensitivity plots.")
    next
  }
  selected <- sens %>% filter(variable %in% c(
    "river_DO","LSR_DO_at_water_table"
  )) %>% mutate(
    source=recode(
      variable,
      river_DO="River entry",
      LSR_DO_at_water_table="Rainfall / LSR at WT"
    )
  )
  if (nrow(selected)) {
    ps <- ggplot(
      selected,aes(
        x=scenario,y=median,colour=source,
        ymin=lower95,ymax=upper95
      )
    )+
      geom_pointrange(
        position=position_dodge(width=.4),
        linewidth=.65
      )+
      facet_wrap(~order_label)+
      scale_colour_manual(
        values=c(
          "River entry"=river_blue,
          "Rainfall / LSR at WT"=rain_red
        ),name="Recharge DO"
      )+
      labs(
        title=paste("River DO boundary sensitivity —",proxy_suite),
        subtitle="Posterior medians and central 95% credible intervals",
        x="River boundary prior specification",
        y="Oxygen concentration (mg/L)",
        caption=paste(
          "River point is a fitted GLOBAL boundary;",
          "site 2801 is a groundwater observation,",
          "not an independent river-water measurement."
        )
      )+base_theme()+
      theme(axis.text.x=element_text(
        angle=20,hjust=1
      ))
    save_plot(
      ps,paste0("09_boundary_prior_",proxy_suite,".png"),
      width=11,height=6
    )
  }

  rate <- sens %>% filter(
    variable %in% c("aquifer_rate_actual","LSR_DO_deficit")
  ) %>% mutate(
    quantity=recode(
      variable,
      aquifer_rate_actual="Aquifer rate (unit depends on order)",
      LSR_DO_deficit="LSR DO deficit (mg/L)"
    )
  )
  if (nrow(rate)) {
    pr <- ggplot(
      rate,aes(
        x=scenario,y=median,
        ymin=lower95,ymax=upper95,colour=quantity
      )
    )+
      geom_pointrange(linewidth=.65)+
      facet_grid(
        quantity~order_label,
        scales="free_y"
      )+
      scale_colour_manual(
        values=c(
          "Aquifer rate (unit depends on order)"=aquifer_green,
          "LSR DO deficit (mg/L)"=rain_red
        ),name=NULL
      )+
      labs(
        title=paste("Aquifer rates and LSR deficit —",proxy_suite),
        subtitle="Check parameter stability and identifiability as river DO prior widens",
        x="River boundary prior specification",
        y="Posterior median and central 95% interval"
      )+base_theme()+
      theme(axis.text.x=element_text(angle=20,hjust=1))
    save_plot(
      pr,paste0("10_kinetics_prior_",proxy_suite,".png"),
      width=11,height=8
    )
  }
}

# 03–04: Final 39-well grouped 5-fold CV of fixed river DO versus broad
# river DO Normal(7.5,1.5^2) prior. The study samples are paired across
# both kinetic orders; site 2801 is NEVER scored. Plot-only.
for (proxy_suite in "2801_masked_no2130") {
  cvdir05 <- file.path(
    "outputs","05_river_DO_grouped_CV",proxy_suite
  )
  pred_path05 <- file.path(cvdir05,"cv_predictions.csv")
  score_path05 <- file.path(cvdir05,"cv_scores.csv")
  if (!file.exists(pred_path05) ||
      !file.exists(score_path05)) next
  cv05 <- readr::read_csv(
    pred_path05,show_col_types=FALSE
  )
  sc05 <- readr::read_csv(
    score_path05,show_col_types=FALSE
  )
  required_cv05 <- c(
    "scenario","aquifer_order","Site_ID","sample_index",
    "observed_DO","predicted_DO",
    "predictive_lower95","predictive_upper95",
    "posterior_median_LSR_fraction","reported_MRT"
  )
  require_cols(cv05,required_cv05,pred_path05)
  if (!nrow(cv05) || any(cv05$Site_ID==cfg$river_proxy_site))
    stop("CV plotting requires non-empty held-out predictions ",
         "that exclude site 2801 from scoring.")
  is_complete05 <- all(sc05$all_folds_complete) &&
    all(sc05$all_MCMC_pass)
  if (!is_complete05)
    warning("CV results for ",proxy_suite,
            " incomplete or sampler diagnostics failed: ",
            "figures are provisional; do not rank models.")
  cv05_nfolds <- max(sc05$folds_completed,na.rm=TRUE)
  excluded_2130 <- !("2130" %in% cv05$Site_ID)
  cv05 <- cv05 %>% mutate(
    order_label=order_label(aquifer_order),
    scenario=factor(
      scenario,levels=c("A_fixed_7p5","C_prior_sd1p5"),
      labels=c("River fixed 7.5","River prior SD 1.5")
    )
  )
  p05 <- ggplot(
    cv05,aes(
      x=observed_DO,y=predicted_DO,
      colour=posterior_median_LSR_fraction
    )
  )+
    geom_abline(
      slope=1,intercept=0,
      colour=ref_colour,linetype="dashed"
    )+
    geom_linerange(
      aes(ymin=predictive_lower95,
          ymax=predictive_upper95),
      alpha=.45,linewidth=.65
    )+
    geom_point(size=2.5)+
    geom_text(
      data=cv05 %>% filter(Site_ID=="2130"),
      aes(label="2130"),
      colour="grey25",vjust=-.8,size=3,
      show.legend=FALSE
    )+
    facet_grid(order_label~scenario)+
    coord_equal()+fraction_colour()+
    labs(
      title=paste(
        "Grouped held-out DO cross-validation —",proxy_suite
      ),
      subtitle=paste(
        cv05_nfolds,"/5 matched folds: predicted mean and 95%",
        "predictive interval.",
         if (excluded_2130) "2130 excluded." else "2130 labelled.",
        if (is_complete05) "" else "PRELIMINARY."
      ),
      x="Observed held-out DO (mg/L)",
      y="Predicted held-out DO (mg/L)",
      caption=paste(
        "2801 DO",if(startsWith(proxy_suite,"2801_masked"))
        "masked in every training fit" else
        "included in training but never scored",
        "; O18 and Cl retained for all sites."
      )
    )+base_theme()
  save_plot(
    p05,paste0(
      "11_river_prior_CV_",proxy_suite,".png"
    ),width=11,height=10
  )

  # Per-well prediction improvement of broad prior versus fixed boundary.
  # Negative delta absolute error means the broad prior predicts better.
  aa <- cv05 %>% filter(
    scenario=="River fixed 7.5"
  ) %>% select(
    sample_index,Site_ID,order_label,
    reported_MRT,observed_DO,
    predicted_fixed=predicted_DO
  )
  bb <- cv05 %>% filter(
    scenario=="River prior SD 1.5"
  ) %>% select(
    sample_index,order_label,
    predicted_broad=predicted_DO,
    f_LSR=posterior_median_LSR_fraction
  )
  paired05 <- inner_join(
    aa,bb,by=c("sample_index","order_label")
  ) %>% mutate(
    delta_abs_error=
      abs(predicted_broad-observed_DO)-
      abs(predicted_fixed-observed_DO)
  )
  if (nrow(paired05)==2L*length(unique(cv05$sample_index))) {
    pdelta <- ggplot(
      paired05,aes(
        reported_MRT,delta_abs_error,
        colour=f_LSR
      )
    )+
      geom_hline(
        yintercept=0,colour=ref_colour,
        linetype="dashed"
      )+
      geom_point(size=2.8)+
      geom_text(
        data=paired05 %>% filter(Site_ID=="2130"),
        aes(label="2130"),
        colour="grey25",vjust=-.9,
        size=3,show.legend=FALSE
      )+
      facet_wrap(~order_label)+fraction_colour()+
      labs(
        title=paste(
          "Predictive effect of uncertain river DO —",
          proxy_suite
        ),
        subtitle=paste(
          "Change in held-out absolute DO error:",
          "broad prior minus fixed boundary"
        ),
        x="Reported MRT (years)",
        y="Change in absolute error (mg/L)",
        caption=paste(
          "Negative = lower error with prior.",
          "Every well scored once, independent of its own DO."
        )
      )+base_theme()
    save_plot(
      pdelta,
      paste0("12_river_prior_CV_delta_",proxy_suite,".png"),
      width=11,height=6
    )
  }
}

cat("\nPLOTS WRITTEN (no MCMC was run):\n")
cat(paste0("  ",manifest,collapse="\n"),"\n",sep="")
message("Plots should appear in Positron's Plots pane when sourced ",
        "from the R console. Saved PNGs are available regardless.")
