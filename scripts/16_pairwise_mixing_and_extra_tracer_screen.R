# 16_pairwise_mixing_and_extra_tracer_screen.R
#
# Audit the original mixing model and identify other available tracers.
# Main pairwise comparisons:
#   Cl + excess air     -> predict held-out delta18O
#   delta18O + excess air -> predict held-out chloride
#   delta18O + Cl       -> predict held-out excess air
#
# The observed withheld tracer is NEVER used to estimate that site's
# two-tracer fraction. However, endmembers are empirical values derived
# from four study samples, so this is a DIAGNOSTIC, not genuinely external
# or independent cross-validation. For Bayesian held-out-tracer CV with
# refitted endmembers, run script 15 instead.
#
# H2 has negative values and is LIKELY delta2H but needs confirmation
# from the original analytical report. There is no river δ2H endpoint
# in this dataset, so it is NOT added blindly to the source-mixing fit.
# Sodium, magnesium, calcium, sulfate and EC may reflect water chemistry
# and reaction as well as mixing. Redox-reactive analytes are surveyed
# descriptively only; they are NOT source-mixing tracers here.
#
# DO is ONLY used in a secondary sensitivity fit of the apparent k,
# AFTER fixing mixing fractions from the tracer observations. It does
# not feed back into the fractions.
#
# Input: private local data/groundwater.csv. Outputs are gitignored.
library(tidyverse)
setwd(here::here("Oxygen-Reduction-Joint-Model"))

cfg <- list(
  input = "data/groundwater.csv",
  river_endmember_site = "2801",
  land_surface_sites = c("1032", "4294"),
  land_surface_prison_site = "5988 HB Prison",
  land_surface_prison_depth_m = 108,
  max_three_tracer_disagreement = 0.35,
  river_DO_mg_L = 7.5,
  max_LSR_entry_DO_mg_L = 9.5,
  max_k_per_year = 0.5
)

if (!file.exists(cfg$input))
  stop("Missing private data/groundwater.csv")
dir.create("outputs", showWarnings = FALSE, recursive = TRUE)
dir.create("plots", showWarnings = FALSE, recursive = TRUE)

dat <- read_csv(
  cfg$input, show_col_types = FALSE,
  col_types = cols(
    Site_ID = col_character(),
    site = col_character(),
    .default = col_guess()
  )
) %>%
  mutate(sample_index = row_number(), .before = 1)

tracers <- c("O18", "Cl", "xs_Air")
stopifnot(all(c("Site_ID", "Depth", tracers) %in% names(dat)))
river <- dat %>% filter(Site_ID == cfg$river_endmember_site)
lsr <- dat %>% filter(
  Site_ID %in% cfg$land_surface_sites |
  (Site_ID == cfg$land_surface_prison_site &
   Depth == cfg$land_surface_prison_depth_m)
)
stopifnot(
  nrow(river) == 1, nrow(lsr) == 3,
  all(is.finite(unlist(river[tracers]))),
  all(is.finite(as.matrix(lsr[tracers])))
)

R <- vapply(tracers, function(x) river[[x]][1], numeric(1))
L <- vapply(tracers, function(x) mean(lsr[[x]]), numeric(1))
stopifnot(all(abs(L - R) > 1e-9))
endmembers <- tibble(
  tracer = tracers,
  river_endmember = unname(R),
  LSR_endmember = unname(L),
  difference = unname(L - R)
)
write_csv(endmembers, "outputs/mixing_pairwise_endmembers.csv")

is_endpoint <- (
  dat$Site_ID == cfg$river_endmember_site |
  dat$sample_index %in% lsr$sample_index
)

# Fraction from EACH tracer, on the raw scale and after physical clipping.
raw_f <- sapply(tracers, function(t) {
  (dat[[t]] - R[[t]]) / (L[[t]] - R[[t]])
})
clipped_f <- pmax(0, pmin(1, raw_f))
colnames(raw_f) <- tracers
colnames(clipped_f) <- tracers
median_nonmissing <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) median(x) else NA_real_
}
f3 <- apply(clipped_f, 1, median_nonmissing)
num_available <- rowSums(is.finite(raw_f))
disagreement <- apply(raw_f, 1, function(x) {
  x <- x[is.finite(x)]
  if (length(x) >= 2) diff(range(x)) else NA_real_
})
quality <- num_available == 3L & is.finite(disagreement) &
  disagreement <= cfg$max_three_tracer_disagreement

fraction_details <- dat %>%
  select(sample_index, Site_ID, Depth, any_of("source")) %>%
  mutate(
    fLSR_O18_raw = raw_f[, "O18"],
    fLSR_Cl_raw = raw_f[, "Cl"],
    fLSR_air_raw = raw_f[, "xs_Air"],
    fLSR_O18_clipped = clipped_f[, "O18"],
    fLSR_Cl_clipped = clipped_f[, "Cl"],
    fLSR_air_clipped = clipped_f[, "xs_Air"],
    fLSR_three_tracer_median = f3,
    unclipped_tracer_fraction_range = disagreement,
    three_tracer_quality_pass = quality,
    reference_endmember_sample = is_endpoint
  )
write_csv(fraction_details, "outputs/mixing_pairwise_all_fractions.csv")

# Pairwise fractions are equal-weight means of the TWO separately
# clipped fractions: equivalent to their median when there are two.
pairs <- tribble(
  ~pair_name,   ~first, ~second,  ~withheld,
  "O18_Cl",      "O18",  "Cl",      "xs_Air",
  "O18_air",     "O18",  "xs_Air",  "Cl",
  "Cl_air",      "Cl",   "xs_Air",  "O18"
)
pairwise <- pmap_dfr(pairs, function(pair_name, first,
                                    second, withheld) {
  f <- (clipped_f[, first] + clipped_f[, second]) / 2
  predicted <- R[[withheld]] +
    f * (L[[withheld]] - R[[withheld]])
  valid <- is.finite(f) & is.finite(dat[[withheld]])
  tibble(
    sample_index = dat$sample_index[valid],
    Site_ID = dat$Site_ID[valid],
    Depth = dat$Depth[valid],
    pair = pair_name,
    predictor_1 = first,
    predictor_2 = second,
    withheld_tracer = withheld,
    pair_LSR_fraction = f[valid],
    three_tracer_LSR_fraction = f3[valid],
    reference_endmember_sample = is_endpoint[valid],
    observed = dat[[withheld]][valid],
    predicted = predicted[valid],
    prediction_minus_observation = predicted[valid] -
      dat[[withheld]][valid]
  )
})
write_csv(pairwise, "outputs/mixing_pairwise_predictions.csv")

pairwise_scores <- map_dfr(c(TRUE, FALSE), function(exclude_endpoint) {
  z <- pairwise
  if (exclude_endpoint)
    z <- z %>% filter(!reference_endmember_sample)
  z %>%
    group_by(pair, withheld_tracer) %>%
    summarise(
      n = n(),
      RMSE = sqrt(mean(prediction_minus_observation^2)),
      MAE = mean(abs(prediction_minus_observation)),
      bias = mean(prediction_minus_observation),
      observed_sd = sd(observed),
      normalised_RMSE = RMSE / observed_sd,
      median_abs_change_in_LSR_fraction =
        median(abs(pair_LSR_fraction -
                   three_tracer_LSR_fraction)),
      .groups = "drop"
    ) %>%
    mutate(
      exclude_reference_endmembers = exclude_endpoint,
      .before = 1
    )
})
write_csv(pairwise_scores, "outputs/mixing_pairwise_scores.csv")

# Test how pair choice affects the apparent age-related DO coefficient,
# with SAME 3-tracer-screened wells in every fit. LSR fractions are
# fixed; k and common effective LSR DO are fitted jointly to DO.
# This part is a sensitivity analysis, NOT an independent mixing test.
fit_DO_given_f <- function(DO, MRT, f, river_DO, Pmax, kmax) {
  profiled <- function(k) {
    a <- exp(-k * MRT)
    river_component <- (1-f) * river_DO * a
    w <- f * a
    P <- if (sum(w^2) > 0) {
      sum(w * (DO-river_component)) / sum(w^2)
    } else 0
    P <- max(0, min(Pmax, P))
    yhat <- river_component + w * P
    c(P = P, k = k, SSE = sum((DO-yhat)^2))
  }
  # Dense grid protects against profile local minima; refine nearby.
  grid <- seq(0, kmax, length.out = 401L)
  val <- lapply(grid, profiled)
  err <- vapply(val, function(x) x[["SSE"]], numeric(1))
  i <- which.min(err)
  low <- grid[max(1L, i-1L)]
  high <- grid[min(length(grid), i+1L)]
  optimum <- if (low < high) {
    optimize(function(x) profiled(x)[["SSE"]],
             interval = c(low, high), tol = 1e-9)$minimum
  } else grid[i]
  candidates <- lapply(
    unique(c(0, grid[i], optimum, kmax)), profiled
  )
  winner <- which.min(vapply(
    candidates, function(x) x[["SSE"]], numeric(1)
  ))
  candidates[[winner]]
}

k_base <- dat %>%
  filter(
    quality,
    is.finite(DO), DO > 0,
    is.finite(MRT), MRT >= 0
  ) %>%
  select(sample_index, Site_ID, DO, MRT) %>%
  mutate(fLSR = f3[match(sample_index, dat$sample_index)])

DO_sensitivity <- map_dfr(
  c("three_tracer_median", pairs$pair_name),
  function(method) {
    z <- k_base
    if (method != "three_tracer_median") {
      zz <- pairwise %>%
        filter(pair == method) %>%
        select(sample_index, pair_LSR_fraction)
      z <- z %>%
        left_join(zz, by = "sample_index") %>%
        mutate(fLSR = pair_LSR_fraction) %>%
        select(-pair_LSR_fraction)
    }
    stopifnot(all(is.finite(z$fLSR)))
    fit <- fit_DO_given_f(
      z$DO, z$MRT, z$fLSR,
      cfg$river_DO_mg_L, cfg$max_LSR_entry_DO_mg_L,
      cfg$max_k_per_year
    )
    predict <- function(held, model) {
      ((1-held$fLSR)*cfg$river_DO_mg_L +
         held$fLSR*model[["P"]]) *
        exp(-model[["k"]]*held$MRT)
    }
    LOOWO <- vapply(seq_len(nrow(z)), function(i) {
      train <- z[-i, ]
      m <- fit_DO_given_f(
        train$DO, train$MRT, train$fLSR,
        cfg$river_DO_mg_L, cfg$max_LSR_entry_DO_mg_L,
        cfg$max_k_per_year
      )
      predict(z[i, ], m)
    }, numeric(1))
    tibble(
      fraction_method = method,
      n_same_wells = nrow(z),
      apparent_k_per_year = unname(fit[["k"]]),
      fitted_effective_LSR_DO = unname(fit[["P"]]),
      fitted_LSR_DO_loss = cfg$max_LSR_entry_DO_mg_L -
        unname(fit[["P"]]),
      fitted_RMSE_mg_L = sqrt(fit[["SSE"]] / nrow(z)),
      leave_one_well_out_RMSE_mg_L =
        sqrt(mean((z$DO - LOOWO)^2))
    )
  }
)
write_csv(
  DO_sensitivity,
  "outputs/mixing_pairwise_effect_on_DO_model.csv"
)

# Other analytes: assess COVERAGE and empirical redundancy, without
# asserting they behave as conservative source-mixing tracers.
possible <- c(
  "H2", "Na", "Mg", "Ca", "K", "SiO2", "SO4",
  "EC", "Condy", "CaCO3", "Alk", "TDS",
  "B", "F", "NO3N", "Fe", "Mn", "xs_N2",
  "NH3", "CH4_\u00b5molL"
)
possible <- intersect(possible, names(dat))
safe_cor <- function(x, y, method = "pearson") {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 5L || sd(x[ok]) == 0 || sd(y[ok]) == 0) {
    return(NA_real_)
  }
  cor(x[ok], y[ok], method = method)
}
candidates <- map_dfr(possible, function(nm) {
  x <- dat[[nm]]
  tibble(
    candidate = nm,
    n_available = sum(is.finite(x)),
    n_paired_O18 = sum(is.finite(x) & is.finite(dat$O18)),
    n_paired_Cl = sum(is.finite(x) & is.finite(dat$Cl)),
    n_paired_air = sum(is.finite(x) & is.finite(dat$xs_Air)),
    correlation_O18 = safe_cor(x, dat$O18),
    correlation_Cl = safe_cor(x, dat$Cl),
    correlation_air = safe_cor(x, dat$xs_Air),
    Spearman_Cl = safe_cor(x, dat$Cl, "spearman")
  )
})
write_csv(candidates, "outputs/mixing_additional_tracer_screen.csv")

# Existing binary 'source' label is useful only if its provenance is
# independent of the current three tracers: DON'T use as a prior yet.
if ("source" %in% names(dat)) {
  label_summary <- tibble(
    source_label = dat$source,
    fLSR_3tracer = f3,
    quality = quality
  ) %>%
    filter(!is.na(source_label)) %>%
    group_by(source_label) %>%
    summarise(
      n = n(),
      median_fLSR = median(fLSR_3tracer, na.rm = TRUE),
      n_LSR_dominated = sum(fLSR_3tracer >= .5, na.rm = TRUE),
      n_three_tracer_QC_pass = sum(quality),
      .groups = "drop"
    )
  write_csv(
    label_summary, "outputs/mixing_source_label_comparison.csv"
  )
}

theme_set(theme_bw(base_size = 12))
p1 <- ggplot(
  pairwise,
  aes(observed, predicted)
) +
  geom_abline(intercept = 0, slope = 1,
              linetype = "dashed", colour = "grey50") +
  geom_point(
    aes(shape = reference_endmember_sample),
    alpha = .8, size = 2.2
  ) +
  facet_wrap(~ withheld_tracer, scales = "free", ncol = 3) +
  labs(
    title = "Two-tracer source fractions predict the third tracer",
    subtitle = "Fixed empirical endpoints; diagnostic, NOT external CV",
    x = "Observed tracer", y = "Predicted from other two tracers",
    shape = "Endpoint sample"
  )
print(p1)
ggsave(
  "plots/mixing_pairwise_predict_third_tracer.png",
  p1, width = 12, height = 5, dpi = 300
)

p2 <- fraction_details %>%
  select(
    Site_ID, fLSR_O18_clipped,
    fLSR_Cl_clipped, fLSR_air_clipped
  ) %>%
  pivot_longer(
    starts_with("fLSR_"),
    names_to = "tracer", values_to = "fraction"
  ) %>%
  ggplot(aes(tracer, fraction, group = Site_ID)) +
  geom_line(alpha = .2) +
  geom_point(alpha = .55) +
  coord_flip() +
  labs(
    title = "Inconsistency between separate tracer mixing fractions",
    x = "Tracer used alone", y = "LSR recharge fraction"
  )
print(p2)
ggsave(
  "plots/mixing_pairwise_source_fraction_disagreement.png",
  p2, width = 9, height = 6, dpi = 300
)

cat("\nSOURCE-ENDMEMBER COMPOSITIONS:\n")
print(endmembers)
cat("\nPAIRWISE HELD-OUT TRACER DIAGNOSTICS:\n")
print(pairwise_scores, n = Inf)
cat("\nEFFECT ON APPARENT DO AGE COEFFICIENT (SAME WELLS):\n")
print(DO_sensitivity, n = Inf)
cat("\nADDITIONAL ANALYTE COVERAGE / CORRELATION SCREEN:\n")
print(candidates, n = Inf)
cat("\nIMPORTANT INTERPRETIVE LIMITS:\n")
cat(
  "* Endmember samples are used to set the very mixing line being\n",
  "  tested. Excluding them from scores reduces circularity but does\n",
  "  not make this diagnostic independently validated.\n",
  "* Binary 'source' labels are NOT independent ground truth unless\n",
  "  their assignment method is confirmed.\n",
  "* H2 is probably delta2H, not dissolved H2; confirm field units.\n",
  "* No river H2 endpoint is available here.\n",
  "* Strong Na/Cl or EC/Cl correlations imply limited independent\n",
  "  mixing information and do NOT prove conservative behavior.\n",
  "* Excess N2, nitrate, Fe, Mn, methane and ammonium are reactive.\n",
  "* Three-tracer agreement and the fitted oxygen age rate depend\n",
  "  on assumed groundwater endmembers, not measured pure recharge.\n",
  "* For proper Bayesian matched held-out predictions, use script 15.\n",
  sep = ""
)
