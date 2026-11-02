# 02_DO_plotting.R — plot-only results and model diagnostics
#
# Run independently AFTER script 01 has generated CSV summaries:
#   source("scripts/02_DO_plotting.R")
#
# NO CmdStan compilation, MCMC sampling, fitting or fit-object loading.
# Uses existing results and setup inputs. It is safe to rerun repeatedly.
# Figures PRINT in the Positron R Plots pane and SAVE as 300-dpi PNGs.
#
# Palette: Turbo throughout. River endmember = Turbo blue;
# rainfall/land-surface-recharge (LSR) endmember = Turbo red.
# For mixed samples, Turbo blue -> red represents LSR fraction (0 -> 1).
# Green is used for aquifer loss, which is NOT a river endmember.

if (!file.exists("scripts/00_setup.R"))
  stop("Run from the repository root (missing scripts/00_setup.R).")
# 00_setup only validates and loads data; it never runs MCMC.
source("scripts/00_setup.R", local=TRUE)

if (!requireNamespace("viridisLite",quietly=TRUE))
  stop("Install viridisLite for the Turbo palette.")

post_file <- output("global_posteriors.csv")
well_file <- output("well_posteriors.csv")
draw_file <- output("joint_draws.csv")
required_files <- c(post_file,well_file,draw_file)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files))
  stop(
    "Missing saved posterior CSV files:\n",
    paste(missing_files,collapse="\n"),
    "\nRun scripts/01_DO_dual_pathway_kinetics.R first, ",
    "then rerun 02. No MCMC is needed to regenerate existing plots."
  )

global <- readr::read_csv(post_file,show_col_types=FALSE)
well <- readr::read_csv(well_file,show_col_types=FALSE)
joint <- readr::read_csv(draw_file,show_col_types=FALSE)

require_cols <- function(df,columns,filename) {
  missed <- setdiff(columns,names(df))
  if (length(missed))
    stop("Unexpected format in ",filename,": missing ",
         paste(missed,collapse=", "))
}
require_cols(global,c("aquifer_order","variable","median",
                      "lower95","upper95"),post_file)
# fit$summary() in 01 writes q5/q95 for wells. These are 90% credible
# intervals, NOT 95%: do not mislabel them or treat them as prediction intervals.
require_cols(well,c("aquifer_order","sample_index","Site_ID",
                    "variable","median","q5","q95",
                    "DO_original","MRT"),well_file)
require_cols(joint,c("aquifer_order","LSR_DO_deficit_mgL",
                     "aquifer_rate","rate_units"),draw_file)

if (!setequal(unique(as.integer(global$aquifer_order)),0:1) ||
    !setequal(unique(as.integer(well$aquifer_order)),0:1) ||
    !setequal(unique(as.integer(joint$aquifer_order)),0:1))
  stop("Posterior summaries must contain BOTH aquifer kinetic orders (0 and 1).")

status_path <- output("fit_status.csv")
if (file.exists(status_path)) {
  fit_status <- readr::read_csv(status_path,show_col_types=FALSE)
  if ("all_diagnostics_pass" %in% names(fit_status) &&
      any(!fit_status$all_diagnostics_pass,na.rm=TRUE))
    warning("A full model failed MCMC checks. Review fit_status.csv ",
            "before interpreting any of these figures.")
}

# Turbo: blue at ~12% of the scale, green at 50%, red at ~90%.
turbo <- viridisLite::turbo(256L)
river_blue <- turbo[32L]
rain_red <- turbo[230L]
aquifer_green <- turbo[128L]
# Restrict the continuous gradient to actual blue -> red endmember colours,
# rather than the dark purple extrema of the full Turbo palette.
turbo_gradient <- turbo[32:230]
ref_colour <- "grey55"   # reference lines only, not scientific data

fraction_colour <- function() {
  scale_colour_gradientn(
    colours=turbo_gradient,limits=c(0,1),
    name="LSR fraction",breaks=c(0,.25,.5,.75,1)
  )
}
fraction_fill <- function() {
  scale_fill_gradientn(
    colours=turbo_gradient,limits=c(0,1),
    name="LSR fraction",breaks=c(0,.25,.5,.75,1)
  )
}
order_label <- function(x)
  factor(as.integer(x),levels=0:1,
         labels=c("Zero-order aquifer","First-order aquifer"))
base_theme <- function() {
  theme_bw(base_size=12)+
    theme(
      panel.grid.minor=element_blank(),
      legend.position="bottom",
      strip.text=element_text(face="bold"),
      plot.title=element_text(face="bold"),
      plot.subtitle=element_text(size=10),
      plot.caption=element_text(hjust=0)
    )
}

dir.create(plotdir,recursive=TRUE,showWarnings=FALSE)
manifest <- character()
save_plot <- function(p,filename,width=10,height=6) {
  # Explicit print is needed when source() is run non-interactively.
  print(p)
  path <- file.path(plotdir,filename)
  ggplot2::ggsave(
    filename=path,plot=p,width=width,height=height,
    units="in",dpi=300,bg="white"
  )
  manifest <<- c(manifest,path)
  message("PLOT: ",normalizePath(path))
  invisible(path)
}

# Pivot well summaries to one row per well per kinetic order.
# Values are posterior summaries of fitted MEAN DO, not replicated DO.
well_vars <- c(
  "f_LSR","mean_DO","true_MRT",
  "vadose_LSR_DO_loss","aquifer_DO_loss"
)
# Stan/CmdStanR reports array elements as "f_LSR[1]", "mean_DO[1]",
# etc. The *base* variable name is the key for pivot_wider(), while the
# bracket index identifies the well. Filtering on the raw indexed name
# silently returns zero rows, which was the source of the missing columns.
indexed_well <- well %>%
  mutate(
    parameter=sub("\\[[0-9]+\\]$", "", variable),
    index_from_variable=suppressWarnings(as.integer(
      stringr::str_match(variable, "\\[([0-9]+)\\]$")[,2]
    ))
  ) %>%
  filter(parameter %in% well_vars)

if (nrow(indexed_well)==0L)
  stop("No per-well posterior variables found. Expected indexed names ",
       "such as f_LSR[1] and mean_DO[1] in ",well_file)
if (anyNA(indexed_well$index_from_variable) ||
    any(indexed_well$index_from_variable!=indexed_well$sample_index))
  stop("Inconsistent well index: variable brackets and sample_index ",
       "disagree in ",well_file)
if (anyDuplicated(indexed_well %>%
                  select(aquifer_order,sample_index,parameter)))
  stop("Duplicate well-parameter records in ",well_file)

wide <- indexed_well %>%
  select(aquifer_order,sample_index,Site_ID,DO_original,
         MRT,parameter,median,q5,q95) %>%
  tidyr::pivot_wider(
    names_from=parameter,values_from=c(median,q5,q95),
    names_glue="{.value}_{parameter}"
  ) %>%
  mutate(order_label=order_label(aquifer_order))
expected <- c(
  "median_f_LSR","q5_f_LSR","q95_f_LSR",
  "median_mean_DO","q5_mean_DO","q95_mean_DO",
  "median_vadose_LSR_DO_loss","median_aquifer_DO_loss",
  "q5_vadose_LSR_DO_loss","q95_vadose_LSR_DO_loss",
  "q5_aquifer_DO_loss","q95_aquifer_DO_loss"
)
require_cols(wide,expected,"pivoted well posterior summary")
if (nrow(wide)!=2*nrow(d))
  stop("Well posterior table does not contain exactly one record ",
       "per site and kinetic order. Check prior script 01 output.")

modelled <- wide %>% filter(Site_ID!=cfg$river_proxy_site)
if (nrow(modelled)!=2*length(eligible))
  stop("Unexpected number of fitted wells in plotting data.")

# 01: observed versus fitted DO, including marginal 90% credible
# intervals for the fitted mean. Never treat the proxy 2801 as a target.
p_do <- ggplot(
  modelled,
  aes(x=DO_original,y=median_mean_DO,colour=median_f_LSR)
)+
  geom_abline(
    slope=1,intercept=0,linetype="dashed",
    colour=ref_colour,linewidth=.6
  )+
  geom_linerange(
    aes(ymin=q5_mean_DO,ymax=q95_mean_DO),
    alpha=.5,linewidth=.6
  )+
  geom_point(size=2.7,alpha=.9)+
  facet_wrap(~order_label)+fraction_colour()+
  coord_equal()+
  labs(
    title="Observed versus fitted groundwater DO",
    subtitle="Points = posterior mean-DO medians; bars = 90% credible intervals for fitted mean",
    x="Observed DO (mg/L)",y="Fitted mean DO (mg/L)",
    caption=paste(
      "Site",cfg$river_proxy_site,"DO excluded from fitting;",
      "in-sample fit, not held-out prediction."
    )
  )+base_theme()
save_plot(p_do,"01_observed_vs_fitted_DO.png",width=11,height=6)

# 02: the two-endmember mixing line in delta18O/chloride space.
# Physical endmembers are explicitly marked in Turbo blue and red.
chem <- modelled %>%
  left_join(d %>% select(sample_index,O18,Cl),
            by="sample_index") %>%
  filter(is.finite(O18),is.finite(Cl))
chem_vars <- c("O18_R","O18_LSR","Cl_R","Cl_LSR")
ends <- global %>%
  filter(variable %in% chem_vars) %>%
  select(aquifer_order,variable,median) %>%
  tidyr::pivot_wider(names_from=variable,values_from=median) %>%
  mutate(order_label=order_label(aquifer_order))
require_cols(ends,chem_vars,"global posterior endmember table")

p_mix <- ggplot()+
  geom_segment(
    data=ends,
    aes(x=Cl_R,y=O18_R,xend=Cl_LSR,yend=O18_LSR),
    colour=aquifer_green,linetype="dashed",linewidth=.9
  )+
  geom_point(
    data=chem,
    aes(x=Cl,y=O18,colour=median_f_LSR),
    size=2.7,alpha=.87
  )+
  geom_point(
    data=ends,
    aes(x=Cl_R,y=O18_R),inherit.aes=FALSE,
    colour=river_blue,size=4,shape=18
  )+
  geom_point(
    data=ends,
    aes(x=Cl_LSR,y=O18_LSR),inherit.aes=FALSE,
    colour=rain_red,size=4,shape=18
  )+
  geom_text(
    data=ends,
    aes(x=Cl_R,y=O18_R),inherit.aes=FALSE,
    label="River",hjust=0,vjust=1.6,
    colour=river_blue,fontface="bold",size=3.5
  )+
  geom_text(
    data=ends,
    aes(x=Cl_LSR,y=O18_LSR),inherit.aes=FALSE,
    label="Rainfall / LSR",hjust=1,vjust=-.8,
    colour=rain_red,fontface="bold",size=3.5
  )+
  facet_wrap(~order_label)+fraction_colour()+
  labs(
    title="Two-endmember conservative-tracer mixing",
    subtitle="River endmember = blue; rainfall/LSR = red; line connects posterior median endmembers",
    x="Chloride (mg/L)",y=expression(delta^18*O~("\u2030")),
    caption="Well colours indicate inferred LSR fraction; endmember positions are posterior medians."
  )+base_theme()
save_plot(p_mix,"02_O18_chloride_mixing.png",width=11,height=6)

# 03: source fractions for individual wells, sorted by posterior
# mean LSR fraction across both kinetic orders. q5/q95 = 90% interval.
rank <- modelled %>%
  group_by(sample_index,Site_ID) %>%
  summarise(f=mean(median_f_LSR),.groups="drop") %>%
  arrange(f,sample_index) %>%
  mutate(site_label=paste0(Site_ID," [",sample_index,"]"))
frac <- modelled %>%
  mutate(site_label=paste0(Site_ID," [",sample_index,"]")) %>%
  left_join(rank %>% select(sample_index,rank_fraction=f),
            by="sample_index")
frac$site_label <- factor(frac$site_label,levels=rank$site_label)
p_fraction <- ggplot(
  frac,aes(y=site_label,x=median_f_LSR,
           colour=median_f_LSR)
)+
  geom_vline(
    xintercept=.5,colour=ref_colour,linetype="dashed"
  )+
  geom_segment(
    aes(x=q5_f_LSR,xend=q95_f_LSR,
        y=site_label,yend=site_label),
    alpha=.65,linewidth=.7
  )+
  geom_point(size=2.4)+
  facet_wrap(~order_label)+fraction_colour()+
  scale_x_continuous(limits=c(0,1),
                     breaks=seq(0,1,.2))+
  labs(
    title="Well-specific land-surface recharge fractions",
    subtitle="Posterior median and 90% credible interval; sorted consistently across kinetic orders",
    x="Fraction of recharge from rainfall / LSR",
    y="Sampling site [sample index]"
  )+base_theme()+
  theme(axis.text.y=element_text(size=7.7))
save_plot(
  p_fraction,"03_LSR_fraction_by_well.png",
  width=12,height=max(9,min(17,nrow(rank)*.32))
)

# 04: joint posterior identifiability. Turbocolour used on the sampled
# posterior points, with a neutral conditional smoother.
id_draws <- joint %>%
  mutate(
    order_label=order_label(aquifer_order),
    LSR_DO_deficit_mgL=as.numeric(LSR_DO_deficit_mgL),
    aquifer_rate=as.numeric(aquifer_rate)
  )
p_ident <- ggplot(
  id_draws,
  aes(x=LSR_DO_deficit_mgL,y=aquifer_rate,
      colour=LSR_DO_deficit_mgL)
)+
  geom_point(alpha=.18,size=.8)+
  geom_smooth(
    aes(group=1),method="loess",se=FALSE,
    colour=aquifer_green,linewidth=1.1
  )+
  facet_wrap(~order_label,scales="free_y")+
  scale_colour_gradientn(
    colours=turbo_gradient,name="LSR DO deficit (mg/L)"
  )+
  labs(
    title="Joint posterior: LSR oxygen deficit versus aquifer rate",
    subtitle="Correlation and trade-off between inferred reaction components",
    x="LSR oxygen deficit to water table (mg/L)",
    y="Aquifer rate (units differ by kinetic order)"
  )+base_theme()
save_plot(p_ident,"04_LSR_aquifer_identifiability.png",width=11,height=6)

# 05: partition of per-well DO loss by stage. Aquifer loss combines
# BOTH recharge sources: use Turbo green, NOT river blue, to avoid
# wrongly implying all aquifer loss occurs in the river endmember.
losses <- bind_rows(
  modelled %>% transmute(
    aquifer_order,order_label,Site_ID,MRT,
    pathway="LSR: pre-water-table",
    loss=median_vadose_LSR_DO_loss,
    lower=q5_vadose_LSR_DO_loss,
    upper=q95_vadose_LSR_DO_loss
  ),
  modelled %>% transmute(
    aquifer_order,order_label,Site_ID,MRT,
    pathway="Aquifer: both sources",
    loss=median_aquifer_DO_loss,
    lower=q5_aquifer_DO_loss,
    upper=q95_aquifer_DO_loss
  )
)
p_losses <- ggplot(
  losses,aes(x=MRT,y=loss,colour=pathway)
)+
  geom_linerange(
    aes(ymin=lower,ymax=upper),
    alpha=.3,linewidth=.6
  )+
  geom_point(size=2.2,alpha=.85)+
  facet_wrap(~order_label)+
  scale_colour_manual(
    values=c(
      "LSR: pre-water-table"=rain_red,
      "Aquifer: both sources"=aquifer_green
    ),name=NULL
  )+
  labs(
    title="Modelled oxygen loss by pathway",
    subtitle="Well-wise posterior medians and 90% credible intervals; aquifer loss includes both sources",
    x="Reported groundwater MRT (years)",
    y="Nominal oxygen loss contribution (mg/L)",
    caption="Unweighted sampled-well contributions, not a catchment oxygen budget."
  )+base_theme()
save_plot(p_losses,"05_vadose_vs_aquifer_loss.png",width=11,height=6)

# 06: physically interpretable source-pathway curves, evaluated
# using POSTERIOR MEDIAN parameters (not predictive credible bands).
curve_par <- global %>%
  filter(variable %in% c(
    "aquifer_rate_actual","LSR_DO_at_water_table"
  )) %>%
  select(aquifer_order,variable,median) %>%
  tidyr::pivot_wider(names_from=variable,values_from=median)
require_cols(curve_par,c("aquifer_rate_actual",
                         "LSR_DO_at_water_table"),
             "curve parameter posterior table")
tmax <- max(d$MRT[is.finite(d$MRT)],na.rm=TRUE)
if (!is.finite(tmax) || tmax<=0) stop("No positive MRT for source curves.")
curves <- tidyr::crossing(
  aquifer_order=0:1,
  years=seq(0,tmax,length.out=350L)
) %>%
  left_join(curve_par,by="aquifer_order") %>%
  mutate(
    order_label=order_label(aquifer_order),
    river=if_else(
      aquifer_order==0L,
      pmax(0,cfg$river_DO_boundary-aquifer_rate_actual*years),
      cfg$river_DO_boundary*exp(-aquifer_rate_actual*years)
    ),
    lsr=if_else(
      aquifer_order==0L,
      pmax(0,LSR_DO_at_water_table-aquifer_rate_actual*years),
      LSR_DO_at_water_table*exp(-aquifer_rate_actual*years)
    )
  ) %>%
  select(order_label,years,river,lsr) %>%
  tidyr::pivot_longer(
    c(river,lsr),names_to="source",
    values_to="DO_mgL"
  ) %>%
  mutate(
    source=recode(
      source,river="River recharge",
      lsr="Rainfall / LSR"
    )
  )

p_decay <- ggplot()+
  geom_line(
    data=curves,
    aes(x=years,y=DO_mgL,colour=source),
    linewidth=1.3
  )+
  geom_point(
    data=modelled,
    aes(x=MRT,y=DO_original,fill=median_f_LSR),
    shape=21,colour=aquifer_green,
    alpha=.80,size=2.5,stroke=.35
  )+
  facet_wrap(~order_label)+
  scale_colour_manual(
    values=c("River recharge"=river_blue,
             "Rainfall / LSR"=rain_red),
    name="Recharge pathway"
  )+
  fraction_fill()+
  labs(
    title="DO depletion along the two recharge pathways",
    subtitle="Blue river and red rainfall/LSR curves use marginal posterior median parameters",
    x="Aquifer residence time (years)",
    y="Dissolved oxygen (mg/L)",
    caption=paste(
      "Points are observed mixed-groundwater DO, coloured by inferred",
      "LSR fraction. Source curves are illustrative, not prediction bands."
    )
  )+base_theme()
save_plot(p_decay,"06_source_pathway_DO_curves.png",width=11,height=6)

# 07: fitted residuals against groundwater MRT for systematic biases.
resid <- modelled %>%
  mutate(fitted_minus_observed=median_mean_DO-DO_original)
p_resid <- ggplot(
  resid,aes(x=MRT,y=fitted_minus_observed,
            colour=median_f_LSR)
)+
  geom_hline(yintercept=0,linetype="dashed",
             colour=ref_colour)+
  geom_point(size=2.6,alpha=.88)+
  facet_wrap(~order_label)+fraction_colour()+
  labs(
    title="DO fit residuals versus groundwater residence time",
    subtitle="In-sample fitted median minus observed DO; check for age-dependent bias",
    x="Reported groundwater MRT (years)",
    y="Fitted minus observed DO (mg/L)"
  )+base_theme()
save_plot(p_resid,"07_DO_residuals_vs_MRT.png",width=11,height=6)

# Optional 08: if grouped cross-validation was run, plot genuinely
# held-out DO predictions (95% PREDICTIVE interval from script 01).
cv_path <- output("cv_predictions.csv")
if (file.exists(cv_path)) {
  cv <- readr::read_csv(cv_path,show_col_types=FALSE)
  require_cols(
    cv,c("aquifer_order","observed_DO","predicted_DO",
         "predictive_lower95","predictive_upper95",
         "fold"),cv_path
  )
  if (nrow(cv)>0) {
    cv <- cv %>% mutate(order_label=order_label(aquifer_order))
    p_cv <- ggplot(cv,aes(x=observed_DO,y=predicted_DO))+
      geom_abline(slope=1,intercept=0,
                  linetype="dashed",colour=ref_colour)+
      geom_linerange(
        aes(ymin=predictive_lower95,
            ymax=predictive_upper95),
        colour=turbo[80L],alpha=.5
      )+
      geom_point(colour=aquifer_green,size=2.4,alpha=.9)+
      facet_wrap(~order_label)+coord_equal()+
      labs(
        title="Grouped cross-validated groundwater DO",
        subtitle="Genuinely held-out DO: mean prediction and 95% predictive interval",
        x="Observed held-out DO (mg/L)",
        y="Predicted held-out DO (mg/L)",
        caption=paste(
          "CV results only; incomplete folds or failed convergence",
          "must not be used to rank kinetic orders."
        )
      )+base_theme()
    save_plot(
      p_cv,"08_grouped_CV_DO_predictions.png",
      width=11,height=6
    )
  }
}

# 09–10: Optional river DO prior sensitivity. Only displayed if the
# experiment 04 has been run. These preserve the blue river / red rainfall
# endmember semantics of the main Turbo-colour plots. No refitting.
for (proxy_suite in c("2801_included","2801_masked")) {
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

# 11–12: Optional grouped 5-fold CV of fixed river DO versus broad
# river DO Normal(7.5,1.5^2) prior. The study samples are paired across
# both kinetic orders; site 2801 is NEVER scored. Plot-only.
for (proxy_suite in c("2801_masked","2801_in_training")) {
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
        "predictive interval. 2130 labelled.",
        if (is_complete05) "" else "PRELIMINARY."
      ),
      x="Observed held-out DO (mg/L)",
      y="Predicted held-out DO (mg/L)",
      caption=paste(
        "2801 DO",if(proxy_suite=="2801_masked")
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
  if (nrow(paired05)==2L*(nrow(d)-1L)) {
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
