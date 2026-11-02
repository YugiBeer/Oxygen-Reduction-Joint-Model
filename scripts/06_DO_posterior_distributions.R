# 06_DO_posterior_distributions.R
#
# Plot posterior distributions from script 04's SAVED CmdStanMCMC fits.
# NO CmdStan compilation, fitting, sampling, or modification of model outputs.
#
# Primary analysis: 39 DO observations (2801 and 2130 DO excluded).
# Four models: fixed/broad river-DO boundary x zero-/first-order aquifer.
#
# Run from repository root:
#   source("scripts/06_DO_posterior_distributions.R")
#
# Produces 300-dpi PNGs plus posterior summary/correlation CSVs in
# plots/04_river_DO_prior/2801_masked_no2130/posterior_distributions/

# --- USER SETTINGS -----------------------------------------------------
suite06 <- "2801_masked_no2130"
scenario06 <- c("A_fixed_7p5", "C_prior_sd1p5")
orders06 <- 0:1
base_size06 <- 12
dpi06 <- 300
# Thin only the joint scatterplot to keep its PNG legible; all available
# posterior draws are retained for density plots, intervals and correlations.
joint_scatter_every06 <- 3L
# ----------------------------------------------------------------------

if (!file.exists("scripts/04_river_DO_prior_sensitivity.R"))
  stop("Set the working directory to the repository root.")

required06 <- c("ggplot2", "dplyr", "tibble", "viridisLite")
missing06 <- required06[
  !vapply(required06, requireNamespace, logical(1), quietly=TRUE)
]
if (length(missing06))
  stop("Please install R packages: ", paste(missing06, collapse=", "))

library(ggplot2)
library(dplyr)

input_dir06 <- file.path("outputs", "04_river_DO_prior", suite06)
plot_dir06 <- file.path("plots", "04_river_DO_prior",
                        suite06, "posterior_distributions")
dir.create(plot_dir06, recursive=TRUE, showWarnings=FALSE)

design06 <- expand.grid(
  scenario=scenario06, aquifer_order=orders06,
  stringsAsFactors=FALSE
)
design06$filename <- file.path(
  input_dir06,
  paste0("fit_",design06$scenario,
         "_order",design06$aquifer_order,".rds")
)
if (any(!file.exists(design06$filename))) {
  stop("Missing saved fit object(s):\n",
       paste(design06$filename[!file.exists(design06$filename)],
             collapse="\n"),
       "\nRun script 04 once with the matching scenario/DO-exclusion settings.")
}

# Refuse to depict unvalidated posteriors without an explicit warning.
diagnostics06 <- file.path(input_dir06, "fit_status.csv")
if (!file.exists(diagnostics06))
  stop("Missing fit_status.csv. Cannot confirm sampler diagnostics.")
status06 <- utils::read.csv(diagnostics06, stringsAsFactors=FALSE)
needed_status06 <- c("scenario", "aquifer_order", "all_diagnostics_pass")
if (!all(needed_status06 %in% names(status06)))
  stop("fit_status.csv does not contain required diagnostic columns.")
for (i in seq_len(nrow(design06))) {
  flag <- status06$all_diagnostics_pass[
    status06$scenario==design06$scenario[i] &
      status06$aquifer_order==design06$aquifer_order[i]
  ]
  if (length(flag)!=1L || is.na(flag) || !flag)
    stop("Sampler diagnostics are incomplete/failed for ",
         design06$scenario[i], " order ",design06$aquifer_order[i],
         ". Review fit_status.csv.")
}

# Stan parameter names; all are scalar and available in the script 04 fits.
vars06 <- c(
  "aquifer_rate_actual",
  "LSR_DO_at_water_table",
  "LSR_DO_deficit",
  "river_DO",
  "sigma_DO",
  "mean_f_LSR"
)

# Important: fit$save_object() writes the CmdStanR object with draws.
# readRDS() + $draws() uses these SAVED draws and does not run Stan.
draws06 <- lapply(seq_len(nrow(design06)), function(i) {
  f <- readRDS(design06$filename[i])
  z <- tryCatch(
    as.data.frame(f$draws(variables=vars06, format="draws_df")),
    error=function(e) stop(
      "Could not load posterior draws from ",
      design06$filename[i], ": ", conditionMessage(e))
  )
  if (!all(vars06 %in% names(z)))
    stop("Missing posterior variable(s) in ",design06$filename[i],
         ": ",paste(setdiff(vars06,names(z)),collapse=", "))
  z <- z[,vars06,drop=FALSE]
  z$scenario <- design06$scenario[i]
  z$aquifer_order <- design06$aquifer_order[i]
  z
})
post06 <- dplyr::bind_rows(draws06)

if (nrow(post06)==0L ||
    any(!is.finite(as.matrix(post06[,vars06,drop=FALSE]))))
  stop("Posterior draws are missing or non-finite.")

# Colours consistent with 02_DO_plotting.R: Turbo blue for fixed river
# boundary, Turbo red for uncertain river boundary.
turbo06 <- viridisLite::turbo(256L)
col06 <- c(A_fixed_7p5=turbo06[32L],
           C_prior_sd1p5=turbo06[230L])
if (!all(scenario06 %in% names(col06)))
  stop("Define a colour for every configured scenario in col06.")

post06 <- post06 %>%
  mutate(
    scenario_label=factor(
      scenario, levels=scenario06,
      labels=c("River DO fixed: 7.5 mg/L",
               "River DO uncertain: SD 1.5 mg/L")
    ),
    order_label=factor(
      aquifer_order,levels=0:1,
      labels=c("Zero-order aquifer",
               "First-order aquifer")
    )
  )

theme06 <- function() {
  theme_bw(base_size=base_size06) +
    theme(
      panel.grid.minor=element_blank(),
      legend.position="none",
      strip.text=element_text(face="bold"),
      plot.title=element_text(face="bold"),
      plot.subtitle=element_text(size=10),
      plot.caption=element_text(hjust=0)
    )
}

save06 <- function(p, filename, width=10.5, height=7.0) {
  print(p) # Display in Positron Plots pane when sourced from console.
  destination <- file.path(plot_dir06, filename)
  ggsave(
    filename=destination,plot=p,width=width,height=height,
    units="in",dpi=dpi06,bg="white"
  )
  message("PLOT: ", normalizePath(destination))
  invisible(destination)
}

# Marginal posterior median and equal-tail 95% CREDIBLE intervals.
summary06 <- post06 %>%
  tidyr::pivot_longer(
    cols=dplyr::all_of(vars06),
    names_to="variable", values_to="value"
  ) %>%
  group_by(scenario, aquifer_order, variable) %>%
  summarise(
    median=median(value),
    lower95=unname(quantile(value, 0.025)),
    upper95=unname(quantile(value, 0.975)),
    .groups="drop"
  ) %>%
  mutate(
    scenario_label=factor(
      scenario,levels=scenario06,
      labels=c("River DO fixed: 7.5 mg/L",
               "River DO uncertain: SD 1.5 mg/L")
    ),
    order_label=factor(
      aquifer_order,levels=0:1,
      labels=c("Zero-order aquifer","First-order aquifer")
    )
  )
utils::write.csv(
  summary06,
  file.path(plot_dir06,"posterior_summary_95CI.csv"),
  row.names=FALSE
)

density06 <- function(variable,title,xlab,filename) {
  dat <- post06
  # River DO is exactly constant (7.5) in scenario A and cannot be
  # represented with kernel density. Handle that case separately below.
  if (identical(variable,"river_DO")) {
    dat <- dat %>% filter(scenario=="C_prior_sd1p5")
    sm <- summary06 %>%
      filter(variable=="river_DO",scenario=="C_prior_sd1p5")
    p <- ggplot(dat,aes(x=river_DO)) +
      geom_density(
        fill=unname(col06["C_prior_sd1p5"]),
        colour=unname(col06["C_prior_sd1p5"]),
        alpha=0.18,linewidth=0.8
      ) +
      geom_vline(
        xintercept=7.5,colour=unname(col06["A_fixed_7p5"]),
        linetype="dashed",linewidth=0.8
      ) +
      geom_segment(
        data=sm,
        aes(x=lower95,xend=upper95,y=0,yend=0),
        inherit.aes=FALSE,linewidth=2.5,
        colour=unname(col06["C_prior_sd1p5"])
      ) +
      geom_point(
        data=sm,aes(x=median,y=0),
        inherit.aes=FALSE,shape=21,size=2.5,
        fill="white",colour=unname(col06["C_prior_sd1p5"])
      ) +
      facet_wrap(~order_label,ncol=2,scales="free_x") +
      labs(title=title,
           subtitle="Red: uncertain river DO posterior; blue dashed: fixed 7.5 mg/L boundary",
           x=xlab,y="Posterior density",
           caption="Baseline: equal-tail 95% credible interval and posterior median. No density is possible for a fixed parameter.")
  } else {
    sm <- summary06 %>% filter(.data$variable==variable)
    p <- ggplot(dat,aes(x=.data[[variable]],
                        colour=scenario,fill=scenario)) +
      geom_density(alpha=0.18,linewidth=0.8,adjust=1.0) +
      geom_segment(
        data=sm,
        aes(x=lower95,xend=upper95,y=0,yend=0,colour=scenario),
        inherit.aes=FALSE,linewidth=2.5
      ) +
      geom_point(
        data=sm,
        aes(x=median,y=0,colour=scenario),
        inherit.aes=FALSE,shape=21,fill="white",size=2.5
      ) +
      scale_colour_manual(values=col06,drop=FALSE) +
      scale_fill_manual(values=col06,drop=FALSE) +
      facet_grid(scenario_label~order_label,scales="free_x") +
      labs(title=title,
           subtitle="Fixed versus uncertain river-entry DO; zero- versus first-order kinetics",
           x=xlab,y="Posterior density",
           caption="Thick baseline: equal-tail 95% credible interval; white circle: posterior median.")
  }
  save06(p+theme06(),filename)
}

density06("aquifer_rate_actual",
          "Aquifer oxygen-consumption rate",
          "Zero-order: mg/L/year; first-order: year^-1",
          "01_aquifer_rate_posterior.png")
density06("LSR_DO_at_water_table",
          "Land-surface recharge DO at the water table",
          "DO (mg/L)",
          "02_LSR_DO_at_water_table_posterior.png")
density06("LSR_DO_deficit",
          "Oxygen deficit along land-surface recharge pathway",
          "DO deficit (mg/L)",
          "03_LSR_DO_deficit_posterior.png")
density06("river_DO",
          "River-entry DO posterior distribution",
          "DO (mg/L)",
          "04_river_DO_posterior.png")
density06("mean_f_LSR",
          "Mean land-surface recharge fraction",
          "Mean fraction of land-surface recharge",
          "05_mean_LSR_fraction_posterior.png")
density06("sigma_DO",
          "Residual groundwater DO dispersion",
          "Model residual scale (mg/L)",
          "06_sigma_DO_posterior.png")

# Joint posterior: examine whether greater LSR loss can be offset by
# a smaller aquifer reaction rate. Correlations use ALL stored draws.
cor06 <- post06 %>%
  group_by(scenario,aquifer_order,scenario_label,order_label) %>%
  summarise(
    r_LSRdeficit_aquiferRate=cor(
      LSR_DO_deficit,aquifer_rate_actual
    ),
    r_riverDO_aquiferRate=if (isTRUE(dplyr::first(scenario)=="A_fixed_7p5"))
      NA_real_
    else
      cor(river_DO,aquifer_rate_actual),
    .groups="drop"
  )
# Fixed river DO has zero variance; its correlation with rate is undefined.
utils::write.csv(
  cor06,
  file.path(plot_dir06,"posterior_correlations.csv"),
  row.names=FALSE
)

joint06 <- post06 %>%
  group_by(scenario,aquifer_order) %>%
  slice(seq(1L,dplyr::n(),by=joint_scatter_every06)) %>%
  ungroup()
joint_fig06 <- ggplot(
  joint06,
  aes(x=LSR_DO_deficit,y=aquifer_rate_actual,colour=scenario)
) +
  geom_point(alpha=0.095,size=0.6) +
  scale_colour_manual(values=col06) +
  facet_grid(order_label~scenario_label,scales="free_y") +
  geom_text(
    data=cor06,
    aes(x=-Inf,y=Inf,
        label=sprintf("r = %.2f",r_LSRdeficit_aquiferRate)),
    inherit.aes=FALSE,
    hjust=-0.12,vjust=1.3,size=4
  ) +
  labs(
    title="Joint posterior: LSR oxygen deficit vs aquifer reaction rate",
    subtitle="Posterior correlations shown in each panel; points are thinned for display only",
    x="LSR DO deficit (mg/L)",
    y="Aquifer rate (mg/L/year for zero-order; year^-1 for first-order)",
    caption="Posterior dependence is not causal evidence. Intervals and correlations use all saved draws."
  ) +
  theme06()
save06(joint_fig06,"07_joint_LSRdeficit_aquifer_rate.png",
       width=10.5,height=7.5)

message("\nPosterior summary (median and 95% CrI):")
print(summary06 %>%
        filter(variable %in% c(
          "aquifer_rate_actual","LSR_DO_deficit",
          "LSR_DO_at_water_table","river_DO"
        )) %>%
        select(scenario,aquifer_order,variable,
               median,lower95,upper95),
      n=Inf,width=Inf)
message("\nPosterior correlations:")
print(cor06 %>%
        select(scenario,aquifer_order,
               r_LSRdeficit_aquiferRate,
               r_riverDO_aquiferRate),
      n=Inf,width=Inf)
message("\nSeven PNGs and two CSVs saved in: ",
        normalizePath(plot_dir06))
message("Plotting only. No Stan sampling or model fitting performed.")
