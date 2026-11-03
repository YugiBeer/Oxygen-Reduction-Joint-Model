# 06_DO_posterior_distributions.R
#
# ONE presentation-ready figure: posterior of first-order aquifer DO rate
# for the final 39-DO-well model with uncertain river-entry DO.
# Reuses saved 04 draws; does NOT run CmdStan, refit, or touch cached CV.
#
# Also writes a CSV of posterior summaries directly from the saved draws,
# which is the authoritative check on potentially stale summary CSVs.
#
# Run at repository root:
#   source("scripts/06_DO_posterior_distributions.R")
#
# Legacy 7-figure comparison code is kept in
# archive/scripts/06_DO_posterior_distributions_legacy.R

if (!file.exists("scripts/04_river_DO_prior_sensitivity.R"))
  stop("Set the working directory to the repository root.")
if (!requireNamespace("ggplot2",quietly=TRUE) ||
    !requireNamespace("posterior",quietly=TRUE) ||
    !requireNamespace("viridisLite",quietly=TRUE))
  stop("Required packages: ggplot2, posterior, viridisLite.")

suite06 <- "2801_masked_no2130"
fit_dir06 <- file.path("outputs","04_river_DO_prior",suite06)
plot_dir06 <- file.path("plots","04_river_DO_prior",
                        suite06,"posterior_distributions")
fit_file06 <- file.path(fit_dir06,"fit_C_prior_sd1p5_order1.rds")
status_file06 <- file.path(fit_dir06,"fit_status.csv")
if (!all(file.exists(c(fit_file06,status_file06))))
  stop("Missing fitted 39-well preferred model or diagnostics. Run script 04 only if the fits are genuinely absent.")

status06 <- utils::read.csv(status_file06,stringsAsFactors=FALSE)
needed06 <- c("scenario","aquifer_order","all_diagnostics_pass")
if (!all(needed06 %in% names(status06)))
  stop("fit_status.csv lacks required diagnostic columns.")
ok06 <- status06[
  status06$scenario=="C_prior_sd1p5" &
    status06$aquifer_order==1L,
  "all_diagnostics_pass"
]
if (length(ok06)!=1L || !isTRUE(ok06))
  stop("Preferred first-order fit has missing/failed sampler diagnostics.")

# Use the actual saved draws, not a normal approximation to the posterior.
fit06 <- readRDS(fit_file06)
vars06 <- c("aquifer_rate_actual","river_DO",
            "LSR_DO_at_water_table","LSR_DO_deficit",
            "sigma_DO","mean_f_LSR")
z06 <- posterior::as_draws_matrix(fit06$draws(variables=vars06))
if (!all(vars06 %in% colnames(z06)))
  stop("Saved posterior lacks expected model quantities.")
if (any(!is.finite(z06[,vars06,drop=FALSE])))
  stop("Non-finite posterior draws detected.")
results06 <- data.frame(
  variable=vars06,
  median=vapply(vars06,function(v) stats::median(z06[,v]),numeric(1)),
  sd=vapply(vars06,function(v) stats::sd(z06[,v]),numeric(1)),
  lower95=vapply(vars06,function(v)
    unname(stats::quantile(z06[,v],0.025)),numeric(1)),
  upper95=vapply(vars06,function(v)
    unname(stats::quantile(z06[,v],0.975)),numeric(1))
)

dir.create(plot_dir06,recursive=TRUE,showWarnings=FALSE)
summary_path06 <- file.path(
  plot_dir06,"preferred_first_order_posterior_summary.csv"
)
utils::write.csv(results06,summary_path06,row.names=FALSE)
print(results06,row.names=FALSE,digits=5)

# Highlight a discrepancy rather than propagating stale CSV numbers.
prior_csv06 <- file.path(fit_dir06,"comparison_summary.csv")
if (file.exists(prior_csv06)) {
  old06 <- utils::read.csv(prior_csv06,stringsAsFactors=FALSE)
  if (all(c("scenario","aquifer_order","variable","median") %in% names(old06))) {
    old06 <- old06[
      old06$scenario=="C_prior_sd1p5" &
        old06$aquifer_order==1L &
        old06$variable=="aquifer_rate_actual",,
      drop=FALSE
    ]
    new_rate06 <- results06$median[
      results06$variable=="aquifer_rate_actual"
    ]
    if (nrow(old06)==1L &&
        abs(old06$median-new_rate06)>0.0002)
      warning(
        "Saved 04 comparison_summary.csv and the saved posterior RDS disagree ",
        "for the first-order rate (CSV=",signif(old06$median,4),
        ", RDS=",signif(new_rate06,4),"). ",
        "Investigate provenance before final publication."
      )
  }
}

k06 <- as.numeric(z06[,"aquifer_rate_actual"])
q06 <- as.numeric(stats::quantile(k06,c(.025,.5,.975)))
density06 <- stats::density(k06,n=2048)
dd06 <- data.frame(x=density06$x,y=density06$y)
ci06 <- dd06[dd06$x>=q06[1] & dd06$x<=q06[3],,drop=FALSE]
primary_color06 <- viridisLite::turbo(256)[230]
p06 <- ggplot2::ggplot(dd06,ggplot2::aes(x=x,y=y)) +
  ggplot2::geom_area(
    data=ci06,fill=primary_color06,alpha=0.28
  ) +
  ggplot2::geom_line(colour=primary_color06,linewidth=1.1) +
  ggplot2::geom_vline(xintercept=q06[2],linewidth=0.8) +
  ggplot2::geom_vline(xintercept=q06[c(1,3)],
                      linetype="dashed",colour="grey40") +
  ggplot2::labs(
    title="First-order aquifer DO reduction rate",
    subtitle=sprintf("Posterior median %.4f yr^-1; 95%% credible interval %.4f-%.4f",
                     q06[2],q06[1],q06[3]),
    x=expression("Aquifer first-order rate (yr"^-1*")"),
    y="Posterior density",
    caption="Uncertain river-entry DO; 39 fitted DO observations. Dashed lines: 95% credible limits."
  ) +
  ggplot2::theme_classic(base_size=13) +
  ggplot2::theme(
    plot.title=ggplot2::element_text(face="bold"),
    plot.caption=ggplot2::element_text(colour="grey35")
  )
print(p06)
plot_path06 <- file.path(plot_dir06,"01_aquifer_rate_posterior.png")
ggplot2::ggsave(
  plot_path06,plot=p06,width=9,height=5.5,
  dpi=300,units="in",bg="white"
)
message("Preferred first-order posterior figure: ",normalizePath(plot_path06))
message("Posterior summary from saved RDS: ",normalizePath(summary_path06))
message("Earlier figures in the same folder may remain from old runs; this script creates one figure only.")
