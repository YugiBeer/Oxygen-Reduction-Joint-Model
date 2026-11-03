# 05_DO_kinetic_trajectories and summary.R — observed DO and river reference trajectories
#
# Plot-only companion to scripts/01_DO_plotting.R. It loads the saved
# script-04 CmdStanR fits and shows observed groundwater DO against reported
# MRT, posterior river endmember reference trajectories, and posterior median
# LSR fraction as the point colour.
# The reference uses the jointly fitted shared aquifer rate and assumed
# effective aquifer age. An LSR trajectory is omitted because its travel time
# and independent loss rate are not determined by this model.
# Prints a four-parameter preferred-model summary at the end and saves one
# Markdown table alongside the figures, using the same saved posterior.
#
# Run from the repository root:
#   source("scripts/05_DO_kinetic_trajectories and summary.R")

# 1. Load inputs and calculate the plotting data and summary table.
# Plotting starts at section 2, after all calculations are complete.
if (!file.exists("scripts/00_setup.R"))
  stop("Set the working directory to the repository root.")
source("scripts/00_setup.R",local=TRUE)

required07 <- c("ggplot2","dplyr","tibble","posterior","viridisLite")
missing07 <- required07[
  !vapply(required07,requireNamespace,logical(1),quietly=TRUE)
]
if (length(missing07))
  stop("Install R package(s): ",paste(missing07,collapse=", "))

suite07 <- "2801_masked_no2130"
scenario07 <- c("A_fixed_7p5","C_prior_sd1p5")
input07 <- file.path("outputs","04_river_DO_prior",suite07)
plotdir07 <- file.path("plots","final_39well")
dir.create(plotdir07,recursive=TRUE,showWarnings=FALSE)

fit_files07 <- expand.grid(
  scenario=scenario07,aquifer_order=0:1,
  stringsAsFactors=FALSE
)
fit_files07$filename <- file.path(
  input07,
  paste0("fit_",fit_files07$scenario,
         "_order",fit_files07$aquifer_order,".rds")
)
if (any(!file.exists(fit_files07$filename)))
  stop("Missing saved fit object(s):\n",
       paste(fit_files07$filename[!file.exists(fit_files07$filename)],
             collapse="\n"),
       "\nRun scripts/02_river_DO_prior_sensitivity.R first.")

status07_path <- file.path(input07,"fit_status.csv")
if (!file.exists(status07_path))
  stop("Missing ",status07_path,"; cannot check sampler diagnostics.")
status07 <- readr::read_csv(status07_path,show_col_types=FALSE)
if (!all(c("scenario","aquifer_order","all_diagnostics_pass") %in%
         names(status07)))
  stop("Unexpected format in ",status07_path)
if (any(!status07$all_diagnostics_pass))
  warning("At least one saved fit failed sampler diagnostics; the plotted ",
          "posterior bands are provisional.")

turbo07 <- viridisLite::turbo(256L)
order_labels07 <- c(`0`="Zero-order aquifer",`1`="First-order aquifer")
scenario_labels07 <- c(
  A_fixed_7p5="River DO fixed at 7.5 mg/L",
  C_prior_sd1p5="River DO prior SD = 1.5 mg/L"
)

dat07 <- dplyr::filter(d,is.finite(DO_original),is.finite(MRT)) %>%
  dplyr::mutate(
    excluded_DO=Site_ID %in% c(cfg$river_proxy_site,"2130"),
    label_hjust=ifelse(MRT==0,-.1,.5)
  )

draw_matrix07 <- function(fit,variables) {
  posterior::as_draws_matrix(fit$draws(variables=variables))
}
quantile_curve07 <- function(x,river) {
  tibble::tibble(
    MRT=x,
    median=stats::median(river),
    lower95=stats::quantile(river,.025),
    upper95=stats::quantile(river,.975)
  )
}

all_points07 <- list()
all_curves07 <- list()
for (i in seq_len(nrow(fit_files07))) {
  fit <- readRDS(fit_files07$filename[i])
  order07 <- fit_files07$aquifer_order[i]
  scenario07_i <- fit_files07$scenario[i]
  n07 <- nrow(dat07)
  fvars07 <- paste0("f_LSR[",seq_len(n07),"]")
  vars07 <- c("river_DO","k_ref",fvars07)
  z07 <- draw_matrix07(fit,vars07)

  fmat07 <- z07[,fvars07,drop=FALSE]
  points07 <- dat07 %>%
    dplyr::mutate(
      scenario=scenario07_i,
      aquifer_order=order07,
      lsr_fraction=apply(fmat07,2,stats::median)
    )
  all_points07[[length(all_points07)+1L]] <- points07

  grid07 <- seq(0,max(dat07$MRT,na.rm=TRUE)*1.05,length.out=250L)
  river07 <- vector("list",length(grid07))
  for (j in seq_along(grid07)) {
    t07 <- grid07[j]
    if (order07==0L) {
      river07[[j]] <- pmax(0,z07[,"river_DO"]-
        z07[,"k_ref"]*cfg$C_ref*t07)
    } else {
      river07[[j]] <- z07[,"river_DO"]*exp(-z07[,"k_ref"]*t07)
    }
  }
  curves07 <- dplyr::bind_rows(lapply(seq_along(grid07),function(j)
    quantile_curve07(grid07[j],river07[[j]])
  )) %>%
    dplyr::mutate(scenario=scenario07_i,aquifer_order=order07)
  all_curves07[[length(all_curves07)+1L]] <- curves07
}

points07 <- dplyr::bind_rows(all_points07) %>%
  dplyr::mutate(
    order_label=factor(aquifer_order,levels=0:1,
                       labels=unname(order_labels07)),
    scenario_label=factor(scenario,levels=scenario07,
                          labels=unname(scenario_labels07))
  )
curves07 <- dplyr::bind_rows(all_curves07) %>%
  dplyr::mutate(
    order_label=factor(aquifer_order,levels=0:1,
                       labels=unname(order_labels07)),
    scenario_label=factor(scenario,levels=scenario07,
                          labels=unname(scenario_labels07))
  )

# Prepare the preferred-model data and site labels before plotting.
curves07_variable <- dplyr::filter(
  curves07,aquifer_order==1,scenario=="C_prior_sd1p5"
)
points07_variable <- dplyr::filter(
  points07,aquifer_order==1,scenario=="C_prior_sd1p5"
)
labels07 <- dplyr::filter(points07,excluded_DO)
labels07_variable <- dplyr::filter(points07_variable,excluded_DO)

# Calculate the reporting table from the same preferred-model posterior.
summary_fit_file07 <- fit_files07$filename[
  fit_files07$scenario=="C_prior_sd1p5" &
    fit_files07$aquifer_order==1L
]
summary_variables07 <- c(
  "river_DO","aquifer_rate_actual",
  "LSR_DO_at_water_table","LSR_DO_deficit"
)
summary_draws07 <- draw_matrix07(
  readRDS(summary_fit_file07),summary_variables07
)
if (!all(summary_variables07 %in% colnames(summary_draws07)) ||
    any(!is.finite(summary_draws07[,summary_variables07,drop=FALSE])))
  stop("Preferred-model summary requires finite draws for all four parameters.")

summary07 <- data.frame(
  parameter=c("River-entry DO","First-order aquifer rate",
              "LSR DO at water table","LSR DO deficit"),
  variable=summary_variables07,
  unit=c("mg/L","yr\u207b\u00b9","mg/L","mg/L"),
  digits=c(2L,4L,2L,2L),
  median=vapply(summary_variables07,function(v)
    stats::median(summary_draws07[,v]),numeric(1)),
  lower95=vapply(summary_variables07,function(v)
    unname(stats::quantile(summary_draws07[,v],.025)),numeric(1)),
  upper95=vapply(summary_variables07,function(v)
    unname(stats::quantile(summary_draws07[,v],.975)),numeric(1)),
  row.names=NULL
)
summary_table07 <- data.frame(
  Parameter=summary07$parameter,
  `Posterior median`=paste(
    sprintf("%.*f",summary07$digits,summary07$median),summary07$unit
  ),
  `95% credible interval`=paste0(
    sprintf("%.*f",summary07$digits,summary07$lower95),"\u2013",
    sprintf("%.*f",summary07$digits,summary07$upper95)
  ),
  check.names=FALSE
)

# 2. Plot the four-model comparison.
# All appearance settings are directly in this ggplot block.
p07 <- ggplot2::ggplot() +
  ggplot2::geom_ribbon(
    data=curves07,
    ggplot2::aes(x=MRT,ymin=lower95,ymax=upper95),
    fill=turbo07[32L],alpha=.16
  ) +
  ggplot2::geom_line(
    data=curves07,
    ggplot2::aes(x=MRT,y=median),
    colour="grey15",linewidth=.8,linetype="longdash"
  ) +
  ggplot2::geom_point(
    data=points07,
    ggplot2::aes(x=MRT,y=DO_original,colour=lsr_fraction,shape=excluded_DO),
    size=2.8,alpha=.85
  ) +
  ggplot2::geom_text(
    data=labels07,
    ggplot2::aes(x=MRT,y=DO_original,label=Site_ID,hjust=label_hjust),
    colour="grey20",vjust=-.9,size=3.2,show.legend=FALSE
  ) +
  ggplot2::facet_grid(order_label~scenario_label) +
  ggplot2::scale_colour_gradientn(
    colours=turbo07[32:230],limits=c(0,1),name="LSR fraction",
    breaks=c(0,.25,.5,.75,1),labels=c("0%","25%","50%","75%","100%")
  ) +
  ggplot2::scale_shape_manual(
    values=c(`FALSE`=16,`TRUE`=1),name=NULL,
    labels=c(`FALSE`="DO used in fit",`TRUE`="DO excluded from fit")
  ) +
  ggplot2::scale_y_continuous(
    limits=c(0,NA),expand=ggplot2::expansion(mult=c(0,.04))
  ) +
  ggplot2::labs(
    title="Observed groundwater DO and river endmember reference",
    subtitle="River reference: posterior median and central 95% credible band",
    x="Mean residence time (years)",y="Dissolved oxygen (mg/L)",
    caption=NULL
  ) +
  ggplot2::theme_bw(base_size=12) +
  ggplot2::theme(
    panel.grid.minor=ggplot2::element_blank(),
    plot.title=ggplot2::element_text(face="bold"),
    plot.subtitle=ggplot2::element_text(size=10),
    strip.text=ggplot2::element_text(face="bold"),
    legend.position="inside",
    legend.position.inside=c(.98,.98),
    legend.justification=c(1,1),
    legend.direction="vertical",
    legend.box="vertical",
    legend.box.just="left",
    legend.box.background=ggplot2::element_rect(fill="white",colour="grey80"),
    legend.background=ggplot2::element_blank(),
    legend.box.margin=ggplot2::margin(6,6,6,6)
  ) +
  ggplot2::guides(
    colour=ggplot2::guide_colourbar(
      order=1,title.position="top",direction="vertical",
      barwidth=grid::unit(.35,"cm"),barheight=grid::unit(3,"cm")
    ),
    shape=ggplot2::guide_legend(order=2,ncol=1)
  )

outfile07 <- file.path(
  plotdir07,paste0("13_DO_kinetic_trajectories_",suite07,".png")
)
ggplot2::ggsave(
  outfile07,p07,width=12,height=10,units="in",dpi=300,bg="white"
)
print(p07)
message("PLOT: ",normalizePath(outfile07))

# 3. Plot the preferred first-order model.
# This is a separate, self-contained plot: edit its appearance here.
p07_variable <- ggplot2::ggplot() +
  ggplot2::geom_ribbon(
    data=curves07_variable,
    ggplot2::aes(x=MRT,ymin=lower95,ymax=upper95),
    fill=turbo07[32L],alpha=.16
  ) +
  ggplot2::geom_line(
    data=curves07_variable,
    ggplot2::aes(x=MRT,y=median),
    colour="grey15",linewidth=.8,linetype="longdash"
  ) +
  ggplot2::geom_point(
    data=points07_variable,
    ggplot2::aes(x=MRT,y=DO_original,colour=lsr_fraction,shape=excluded_DO),
    size=3.2,alpha=.88
  ) +
  ggplot2::geom_text(
    data=labels07_variable,
    ggplot2::aes(x=MRT,y=DO_original,label=Site_ID,hjust=label_hjust),
    colour="grey20",vjust=-.9,size=3.5,show.legend=FALSE
  ) +
  ggplot2::scale_colour_gradientn(
    colours=turbo07[32:230],limits=c(0,1),name="LSR fraction",
    breaks=c(0,.25,.5,.75,1),labels=c("0%","25%","50%","75%","100%")
  ) +
  ggplot2::scale_shape_manual(
    values=c(`FALSE`=16,`TRUE`=1),name=NULL,
    labels=c(`FALSE`="DO used in fit",`TRUE`="DO excluded from fit")
  ) +
  ggplot2::scale_x_continuous(
    limits=c(0,55),expand=ggplot2::expansion(c(0,0))
  ) +
  ggplot2::scale_y_continuous(
    limits=c(0,NA),expand=ggplot2::expansion(mult=c(0,.04))
  ) +
  ggplot2::labs(
    title="First-order river reference with uncertain river DO",
    x="Mean residence time (years)",y="Dissolved oxygen (mg/L)",
    caption=NULL
  ) +
  ggplot2::theme_bw(base_size=13) +
  ggplot2::theme(
    panel.grid.minor=ggplot2::element_blank(),
    plot.title=ggplot2::element_text(face="bold"),
    legend.position="inside",
    legend.position.inside=c(.98,.98),
    legend.justification=c(1,1),
    legend.direction="vertical",
    legend.box="vertical",
    legend.box.just="left",
    legend.box.background=ggplot2::element_rect(fill="white",colour="grey80"),
    legend.background=ggplot2::element_blank(),
    legend.box.margin=ggplot2::margin(6,6,6,6)
  ) +
  ggplot2::guides(
    colour=ggplot2::guide_colourbar(
      order=1,title.position="top",direction="vertical",
      barwidth=grid::unit(.35,"cm"),barheight=grid::unit(3,"cm")
    ),
    shape=ggplot2::guide_legend(order=2,ncol=1)
  )

outfile07_variable <- file.path(
  plotdir07,
  paste0("14_DO_first_order_variable_river_prior_",suite07,".png")
)
ggplot2::ggsave(
  outfile07_variable,p07_variable,
  width=10,height=8,units="in",dpi=300,bg="white"
)
print(p07_variable)
message("PLOT: ",normalizePath(outfile07_variable))

# 4. Print and save the summary table calculated above.
summary_markdown07 <- c(
  "# Preferred first-order model: posterior summary",
  "",
  "Uncertain river-entry DO: Normal(7.5, 1.5\u00b2) mg/L prior, truncated to [0, 14.65].",
  paste0("Saved fit: `",summary_fit_file07,"`."),
  "",
  paste0("| ",paste(names(summary_table07),collapse=" | ")," |"),
  "| --- | ---: | ---: |",
  apply(summary_table07,1,function(row)
    paste0("| ",paste(row,collapse=" | ")," |")),
  "",
  "Central 95% credible intervals are the 2.5th and 97.5th posterior percentiles; units match the median column.",
  "The aquifer rate is shared by both pathways. The LSR deficit is relative to assumed surface saturation; it is not a loss rate."
)
summary_path07 <- file.path(
  plotdir07,paste0("15_DO_primary_posterior_summary_",suite07,".md")
)
writeLines(enc2utf8(summary_markdown07),summary_path07,useBytes=TRUE)
cat("\nPreferred first-order model with uncertain river-entry DO\n\n")
print(summary_table07,row.names=FALSE,right=FALSE)
message("SUMMARY: ",normalizePath(summary_path07))
