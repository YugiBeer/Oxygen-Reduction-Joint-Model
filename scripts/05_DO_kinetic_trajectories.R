# 07_DO_kinetic_trajectories.R — observed DO and posterior kinetic trajectories
#
# Plot-only companion to scripts/02_DO_plotting.R. It loads the saved
# script-04 CmdStanR fits and shows observed groundwater DO against reported
# MRT, posterior river and LSR endmember trajectories, and posterior median
# LSR fraction as the point colour.
#
# Run from the repository root:
#   source("scripts/07_DO_kinetic_trajectories.R")

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
       "\nRun scripts/04_river_DO_prior_sensitivity.R first.")

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
fraction_colour07 <- function() {
  ggplot2::scale_colour_gradientn(
    colours=turbo07[32:230],limits=c(0,1),
    name="LSR fraction",breaks=c(0,.25,.5,.75,1),
    labels=function(x) paste0(round(100*x),"%")
  )
}
source_colours07 <- c(
  "River endmember"=turbo07[32L],
  "LSR endmember"=turbo07[230L]
)
order_labels07 <- c(`0`="Zero-order aquifer",`1`="First-order aquifer")
scenario_labels07 <- c(
  A_fixed_7p5="River DO fixed at 7.5 mg/L",
  C_prior_sd1p5="River DO prior SD = 1.5 mg/L"
)

dat07 <- dplyr::filter(d,is.finite(DO_original),is.finite(MRT)) %>%
  dplyr::mutate(
    excluded_DO=Site_ID %in% c(cfg$river_proxy_site,"2130")
  )

draw_matrix07 <- function(fit,variables) {
  posterior::as_draws_matrix(fit$draws(variables=variables))
}
quantile_curve07 <- function(x,river,lsr) {
  tibble::tibble(
    MRT=x,
    source=c("River endmember","LSR endmember"),
    median=c(stats::median(river),stats::median(lsr)),
    lower95=c(stats::quantile(river,.025),stats::quantile(lsr,.025)),
    upper95=c(stats::quantile(river,.975),stats::quantile(lsr,.975))
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
  tvars07 <- paste0("true_MRT[",seq_len(n07),"]")
  vars07 <- c("river_DO","LSR_DO_at_water_table","k_ref",
              fvars07,tvars07)
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
  lsr07 <- vector("list",length(grid07))
  for (j in seq_along(grid07)) {
    t07 <- grid07[j]
    if (order07==0L) {
      river07[[j]] <- pmax(0,z07[,"river_DO"]-
        z07[,"k_ref"]*cfg$C_ref*t07)
      lsr07[[j]] <- pmax(0,z07[,"LSR_DO_at_water_table"]-
        z07[,"k_ref"]*cfg$C_ref*t07)
    } else {
      river07[[j]] <- z07[,"river_DO"]*exp(-z07[,"k_ref"]*t07)
      lsr07[[j]] <- z07[,"LSR_DO_at_water_table"]*
        exp(-z07[,"k_ref"]*t07)
    }
  }
  curves07 <- dplyr::bind_rows(lapply(seq_along(grid07),function(j)
    quantile_curve07(grid07[j],river07[[j]],lsr07[[j]])
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
    source=factor(source,levels=names(source_colours07)),
    order_label=factor(aquifer_order,levels=0:1,
                       labels=unname(order_labels07)),
    scenario_label=factor(scenario,levels=scenario07,
                          labels=unname(scenario_labels07))
  )

p07 <- ggplot2::ggplot() +
  ggplot2::geom_ribbon(
    data=curves07,
    ggplot2::aes(x=MRT,ymin=lower95,ymax=upper95,fill=source),
    alpha=.16,colour=NA
  )+
  ggplot2::geom_line(
    data=curves07,
    ggplot2::aes(x=MRT,y=median,linetype=source),
    colour="grey15",linewidth=.9
  )+
  ggplot2::geom_point(
    data=points07,
    ggplot2::aes(x=MRT,y=DO_original,colour=lsr_fraction,
                 shape=excluded_DO),
    size=2.8,alpha=.85
  )+
  ggplot2::geom_text(
    data=points07 %>% dplyr::filter(
      Site_ID %in% c(cfg$river_proxy_site,"2130")
    ),
    ggplot2::aes(x=MRT,y=DO_original,label=Site_ID),
    colour="grey20",vjust=-.9,size=3.2,show.legend=FALSE
  )+
  ggplot2::facet_grid(order_label~scenario_label)+
  ggplot2::scale_fill_manual(values=source_colours07,name="Kinetic source")+
  ggplot2::scale_linetype_manual(
    values=c("River endmember"="solid","LSR endmember"="dashed"),
    name="Kinetic source"
  )+
  ggplot2::scale_shape_manual(
    values=c(`FALSE`=16,`TRUE`=1),
    labels=c(`FALSE`="DO used in fit",`TRUE`="DO excluded from fit"),
    name=NULL
  )+
  fraction_colour07()+
  ggplot2::scale_y_continuous(
    limits=c(0,NA),expand=ggplot2::expansion(mult=c(0,.04))
  )+
  ggplot2::labs(
    title="Observed groundwater DO and fitted two-endmember kinetics",
    subtitle="Curves show posterior median with central 95% credible bands; points use reported MRT",
    x="Mean residence time (years)",
    y="Dissolved oxygen (mg/L)",
    caption="River and LSR curves are endmember DO after aquifer reaction. Point colour is posterior median LSR fraction; open points were excluded from the DO likelihood."
  )+
  ggplot2::theme_bw(base_size=12)+
  ggplot2::theme(
    panel.grid.minor=ggplot2::element_blank(),
    legend.position="bottom",
    strip.text=ggplot2::element_text(face="bold"),
    plot.title=ggplot2::element_text(face="bold"),
    plot.subtitle=ggplot2::element_text(size=10),
    plot.caption=ggplot2::element_text(hjust=0)
  )+
  ggplot2::guides(
    fill=ggplot2::guide_legend(order=1),
    linetype=ggplot2::guide_legend(order=1),
    colour=ggplot2::guide_colourbar(order=2),
    shape=ggplot2::guide_legend(order=3)
  )

outfile07 <- file.path(
  plotdir07,paste0("13_DO_kinetic_trajectories_",suite07,".png")
)
ggplot2::ggsave(
  outfile07,p07,width=12,height=10,units="in",dpi=300,bg="white"
)
print(p07)
message("PLOT: ",normalizePath(outfile07))

# A focused version for the preferred first-order model with the variable
# river-DO prior. This is intentionally a separate file so it can be used
# directly in figures without the other three model panels.
curves07_variable <- curves07 %>%
  dplyr::filter(
    aquifer_order==1,
    scenario=="C_prior_sd1p5"
  )
points07_variable <- points07 %>%
  dplyr::filter(
    aquifer_order==1,
    scenario=="C_prior_sd1p5"
  )

p07_variable <- ggplot2::ggplot() +
  ggplot2::geom_ribbon(
    data=curves07_variable,
    ggplot2::aes(x=MRT,ymin=lower95,ymax=upper95,fill=source),
    alpha=.16,colour=NA
  )+
  ggplot2::geom_line(
    data=curves07_variable,
    ggplot2::aes(x=MRT,y=median,linetype=source),
    colour="grey15",linewidth=1
  )+
  ggplot2::geom_point(
    data=points07_variable,
    ggplot2::aes(x=MRT,y=DO_original,colour=lsr_fraction,
                 shape=excluded_DO),
    size=3.2,alpha=.88
  )+
  ggplot2::geom_text(
    data=points07_variable %>% dplyr::filter(
      Site_ID %in% c(cfg$river_proxy_site,"2130")
    ),
    ggplot2::aes(x=MRT,y=DO_original,label=Site_ID),
    colour="grey20",vjust=-.9,size=3.5,show.legend=FALSE
  )+
  ggplot2::scale_fill_manual(values=source_colours07,name="Kinetic source")+
  ggplot2::scale_linetype_manual(
    values=c("River endmember"="solid","LSR endmember"="dashed"),
    name="Kinetic source"
  )+
  ggplot2::scale_shape_manual(
    values=c(`FALSE`=16,`TRUE`=1),
    labels=c(`FALSE`="DO used in fit",`TRUE`="DO excluded from fit"),
    name=NULL
  )+
  fraction_colour07()+
  ggplot2::scale_x_continuous(
    limits=c(0,55),expand=ggplot2::expansion(c(0,0))
  )+
  ggplot2::scale_y_continuous(
    limits=c(0,NA),expand=ggplot2::expansion(mult=c(0,.04))
  )+
  ggplot2::labs(
    title="First-order aquifer kinetics with variable river DO",
    #subtitle="River DO prior: Normal(7.5, 1.5²)",
    x="Mean residence time (years)",
    y="Dissolved oxygen (mg/L)",
    caption="Curves show posterior median with central 95% credible bands after aquifer reaction. Point colour is posterior median LSR fraction."
  )+
  ggplot2::theme_bw(base_size=13)+
  ggplot2::theme(
    panel.grid.minor=ggplot2::element_blank(),
    legend.position="bottom",
    plot.title=ggplot2::element_text(face="bold"),
    plot.subtitle=ggplot2::element_text(size=10),
    plot.caption=ggplot2::element_text(hjust=0)
  )+
  ggplot2::guides(
    fill=ggplot2::guide_legend(order=1),
    linetype=ggplot2::guide_legend(order=1),
    colour=ggplot2::guide_colourbar(order=2),
    shape=ggplot2::guide_legend(order=3)
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

