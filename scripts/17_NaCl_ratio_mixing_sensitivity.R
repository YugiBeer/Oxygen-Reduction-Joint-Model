# 17_NaCl_ratio_mixing_sensitivity.R
#
# Test whether replacing chloride concentration with the sodium-to-
# chloride MASS RATIO improves conservative source mixing.
# Tracer endmembers from site 2801 (river) and mean of sites
# 1032,4294,5988 HB Prison depth 108m (LSR). The endmembers are study
# groundwater samples, not external pure-recharge samples.
#
# CRITICAL: Na/Cl is a RATIO OF MIXED CONCENTRATIONS. The endmember
# ratios cannot be linearly interpolated. Exact binary mass balance:
#
# q = Na_mix/Cl_mix
#   = [Na_R + f*(Na_L-Na_R)]/[Cl_R + f*(Cl_L-Cl_R)].
#
# f = [Na_R - q*Cl_R] / [q*(Cl_L-Cl_R)-(Na_L-Na_R)].
#
# Ratio-derived fractions can be unstable and can fall outside [0,1].
# Report RAW ratios/fractions and physically clipped versions separately.
# Ratios cannot identify exchange or chloride enrichment uniquely.
#
# Independent tracer-only pair comparisons:
#    Cl + excess air       -> predicts isotope
#    Na/Cl + excess air    -> predicts isotope
#    isotope + Cl          -> predicts excess air
#    isotope + Na/Cl       -> predicts excess air
#
# Values of a withheld tracer are excluded from THAT WELL'S mixing
# estimate, but the source endpoints come from the same dataset.
# Therefore these are empirical sensitivity tests, NOT external CV.
#
# Finally fit a first-order DO model under different FIXED source
# fractions, always using the SAME previously screened wells.
# DO is not used to infer source fractions.
#
# To assess uncertainty in endmembers, mixing and source DO jointly,
# a new Bayesian model should use Na and Cl as two observed ions
# with correlated process errors. Do not treat Na/Cl and Cl as two
# statistically independent tracers (shared denominator).
library(tidyverse)
setwd(here::here("Oxygen-Reduction-Joint-Model"))

cfg <- list(
  input = "data/groundwater.csv",
  river_site = "2801",
  LSR_sites = c("1032", "4294"),
  LSR_prison_site = "5988 HB Prison",
  LSR_prison_depth = 108,
  three_tracer_disagreement_threshold = 0.35,
  river_DO = 7.5,
  maximum_effective_LSR_DO = 9.5,
  maximum_k = .5
)
if (!file.exists(cfg$input)) stop("Missing private data/groundwater.csv")
dir.create("outputs", showWarnings = FALSE, recursive = TRUE)
dir.create("plots", showWarnings = FALSE, recursive = TRUE)

x <- read_csv(
  cfg$input, show_col_types = FALSE,
  col_types = cols(Site_ID = col_character(), .default = col_guess())
) %>%
  mutate(sample_index = row_number(), .before = 1)
analytes <- c("Na", "Cl", "O18", "xs_Air")
stopifnot(all(c(analytes, "Site_ID", "Depth", "DO", "MRT") %in% names(x)))

river <- x %>% filter(Site_ID == cfg$river_site)
lsr <- x %>% filter(
  Site_ID %in% cfg$LSR_sites |
  (Site_ID == cfg$LSR_prison_site &
     Depth == cfg$LSR_prison_depth)
)
stopifnot(nrow(river) == 1, nrow(lsr) == 3)
R <- vapply(analytes, function(t) river[[t]][1], numeric(1))
L <- vapply(analytes, function(t) mean(lsr[[t]]), numeric(1))
stopifnot(all(is.finite(R)), all(is.finite(L)),
          all(L != R), all(x$Cl[is.finite(x$Cl)] > 0))
endmembers <- tibble(
  analyte = analytes, river = unname(R), LSR = unname(L)
)
endmembers <- bind_rows(
  endmembers,
  tibble(
    analyte = "Na/Cl",
    river = R[["Na"]] / R[["Cl"]],
    LSR = L[["Na"]] / L[["Cl"]]
  )
)
write_csv(endmembers, "outputs/NaCl_mixing_endmembers.csv")

# Source fractions from original THREE individually measured tracers.
tracers <- c("O18", "Cl", "xs_Air")
raw_3 <- sapply(tracers, function(t) {
  (x[[t]] - R[[t]]) / (L[[t]] - R[[t]])
})
stopifnot(is.matrix(raw_3), nrow(raw_3) == nrow(x))
clip <- function(a) pmax(pmin(a, 1), 0) # preserve matrix dim
f_individual <- clip(raw_3)
n_available <- rowSums(is.finite(raw_3))
range_3 <- apply(raw_3, 1, function(z) {
  z <- z[is.finite(z)]
  if (length(z) >= 2) diff(range(z)) else NA_real_
})
median_valid <- function(z) {
  z <- z[is.finite(z)]
  if (length(z) > 0) median(z) else NA_real_
}
f_3 <- apply(f_individual, 1, median_valid)
good_3 <- n_available == 3L &
  is.finite(range_3) &
  range_3 <= cfg$three_tracer_disagreement_threshold

# Mass-balanced sodium/chloride ratio.
q <- x$Na / x$Cl
dNa <- L[["Na"]] - R[["Na"]]
dCl <- L[["Cl"]] - R[["Cl"]]
denom <- q * dCl - dNa
f_ratio_raw <- (R[["Na"]] - q * R[["Cl"]]) / denom
f_ratio_raw[!is.finite(q) |
              !is.finite(f_ratio_raw) |
              abs(denom) < 1e-10] <- NA_real_
f_ratio <- clip(f_ratio_raw)
ratio_below_end <- is.finite(f_ratio_raw) &
  (f_ratio_raw < 0 | f_ratio_raw > 1)

f_Na <- clip((x$Na-R[["Na"]]) / dNa)
f_Cl <- f_individual[, "Cl"]
f_O18 <- f_individual[, "O18"]
f_air <- f_individual[, "xs_Air"]

fractions <- tibble(
  sample_index = x$sample_index,
  Site_ID = x$Site_ID,
  Depth = x$Depth,
  Na = x$Na, Cl = x$Cl, NaCl_ratio = q,
  f_LSR_Na = f_Na,
  f_LSR_Cl = f_Cl,
  f_LSR_NaCl_ratio_RAW = f_ratio_raw,
  f_LSR_NaCl_ratio_CLIPPED = f_ratio,
  f_LSR_O18 = f_O18,
  f_LSR_excess_air = f_air,
  f_LSR_original_three_tracers = f_3,
  original_three_tracer_QC = good_3,
  NaCl_outside_endmember_range = ratio_below_end,
  NaCl_fraction_minus_Cl_fraction = f_ratio-f_Cl
)
write_csv(fractions, "outputs/NaCl_mixing_fractions.csv")

# Predict a tracer using fractions obtained from the other two.
# This is a fixed-endmember pairwise diagnostic.
pairs <- list(
  list(name="Cl+air", fractions=list(f_Cl,f_air),
       held="O18"),
  list(name="NaCl+air", fractions=list(f_ratio,f_air),
       held="O18"),
  list(name="O18+Cl", fractions=list(f_O18,f_Cl),
       held="xs_Air"),
  list(name="O18+NaCl", fractions=list(f_O18,f_ratio),
       held="xs_Air")
)
pair_predictions <- map_dfr(pairs, function(spec) {
  fraction_pair <- rowMeans(
    do.call(cbind, spec$fractions)
  )
  t <- spec$held
  y <- x[[t]]
  pred <- R[[t]] + fraction_pair*(L[[t]]-R[[t]])
  valid <- is.finite(y) & is.finite(pred)
  tibble(
    sample_index = x$sample_index[valid],
    Site_ID = x$Site_ID[valid],
    Depth = x$Depth[valid],
    mixing_pair = spec$name,
    withheld = t,
    observed = y[valid],
    predicted = pred[valid],
    residual = y[valid]-pred[valid]
  )
})
pair_scores <- pair_predictions %>%
  group_by(mixing_pair, withheld) %>%
  summarise(
    n = n(),
    RMSE = sqrt(mean(residual^2)),
    MAE = mean(abs(residual)),
    bias_predicted_minus_observed = -mean(residual),
    .groups = "drop"
  )
write_csv(pair_predictions, "outputs/NaCl_pairwise_predictions.csv")
write_csv(pair_scores, "outputs/NaCl_pairwise_scores.csv")

# Same baseline three-tracer-QC-selected samples in every conditional
# DO fit, so changes reflect mixing FRACTIONS, not changing inclusion.
base <- x %>%
  select(sample_index, Site_ID, Depth, DO, MRT) %>%
  filter(is.finite(DO), DO > 0, is.finite(MRT), MRT >= 0) %>%
  mutate(qc = good_3[match(sample_index, x$sample_index)]) %>%
  filter(qc)
stopifnot(nrow(base) >= 25L)

fraction_sets <- list(
  "O18+Cl+air" = f_3,
  "O18+NaCl+air" = apply(
    cbind(f_O18,f_ratio,f_air), 1, median_valid
  ),
  "Cl+air" = rowMeans(cbind(f_Cl,f_air)),
  "NaCl+air" = rowMeans(cbind(f_ratio,f_air)),
  "O18+Cl" = rowMeans(cbind(f_O18,f_Cl)),
  "O18+NaCl" = rowMeans(cbind(f_O18,f_ratio))
)

fit_first <- function(dat, f) {
  stopifnot(nrow(dat) == length(f), all(is.finite(f)))
  at_k <- function(k) {
    a <- exp(-k*dat$MRT)
    river_contribution <- (1-f)*cfg$river_DO*a
    w <- f*a
    P <- if (sum(w^2) > 0) {
      sum(w*(dat$DO-river_contribution))/sum(w^2)
    } else 0
    P <- min(cfg$maximum_effective_LSR_DO,
             max(0,P))
    pred <- river_contribution+w*P
    c(k=k,P=P,SSE=sum((pred-dat$DO)^2))
  }
  grid <- seq(0,cfg$maximum_k,length.out=501)
  scored <- lapply(grid, at_k)
  sse <- vapply(scored,function(v) v[["SSE"]],numeric(1))
  im <- which.min(sse)
  lo <- grid[max(1,im-1)]
  hi <- grid[min(length(grid),im+1)]
  best <- if (lo < hi) {
    optimize(function(k) at_k(k)[["SSE"]],
             c(lo,hi),tol=1e-10)$minimum
  } else grid[im]
  cand <- lapply(unique(c(0,grid[im],best,
                         cfg$maximum_k)),at_k)
  cand[[which.min(vapply(
    cand,function(v) v[["SSE"]],numeric(1)
  ))]]
}

DO_results <- imap_dfr(fraction_sets,function(full_f,method) {
  f <- full_f[match(base$sample_index, x$sample_index)]
  if (!all(is.finite(f))) stop("Missing fraction for ",method)
  full <- fit_first(base,f)
  LOO <- vapply(seq_len(nrow(base)),function(i) {
    m <- fit_first(base[-i,],f[-i])
    initial <- (1-f[i])*cfg$river_DO+f[i]*m[["P"]]
    initial*exp(-m[["k"]]*base$MRT[i])
  },numeric(1))
  tibble(
    mixing_method=method,
    n_same_wells=nrow(base),
    apparent_k_per_year=unname(full[["k"]]),
    effective_LSR_DO=unname(full[["P"]]),
    inferred_LSR_DO_loss=
      cfg$maximum_effective_LSR_DO-unname(full[["P"]]),
    fitted_RMSE=sqrt(full[["SSE"]]/nrow(base)),
    leave_one_well_out_RMSE=sqrt(mean(
      (base$DO-LOO)^2
    ))
  )
})
write_csv(DO_results,"outputs/NaCl_effect_on_DO_model.csv")

# Raw-ratio mixing estimates near either endpoint are not trustworthy
# without uncertainty propagation, and ratios share the Cl error.
p1 <- fractions %>%
  filter(is.finite(f_LSR_Cl),
         is.finite(f_LSR_NaCl_ratio_CLIPPED)) %>%
  ggplot(
    aes(f_LSR_Cl,f_LSR_NaCl_ratio_CLIPPED)
  ) +
  geom_abline(slope=1,intercept=0,linetype="dashed") +
  geom_point(aes(shape=NaCl_outside_endmember_range),
             size=2.2) +
  coord_equal(xlim=c(0,1),ylim=c(0,1)) +
  labs(
    title="Mixing fraction: chloride vs Na/Cl",
    subtitle="Physical clipping shown; out-of-range raw ratios flagged",
    x="LSR fraction from Cl alone",
    y="LSR fraction from mass-balanced Na/Cl",
    shape="Outside endpoint ratio"
  ) + theme_bw(base_size=12)
print(p1)
ggsave("plots/NaCl_vs_Cl_source_fractions.png",
       p1,width=7,height=6,dpi=300)

p2 <- pair_predictions %>%
  ggplot(aes(observed,predicted)) +
  geom_abline(slope=1,intercept=0,linetype="dashed") +
  geom_point(alpha=.75) +
  facet_grid(withheld~mixing_pair,scales="free") +
  labs(
    title="Can Na/Cl improve prediction of the omitted tracer?",
    subtitle="Fixed empirical endmembers; tests are not external CV",
    x="Measured tracer",y="Predicted tracer"
  ) + theme_bw(base_size=11)
print(p2)
ggsave("plots/NaCl_predict_heldout_tracers.png",
       p2,width=12,height=7,dpi=300)

cat("\nENDMEMBERS:\n")
print(endmembers,n=Inf)
cat("\nOUT-OF-RANGE Na/Cl FRACTIONS:\n")
print(tibble(
  n_ratio_measured=sum(is.finite(f_ratio_raw)),
  n_ratio_outside_0_1=sum(ratio_below_end,na.rm=TRUE),
  median_abs_difference_cl_vs_ratio=median(
    abs(f_ratio-f_Cl),na.rm=TRUE
  )
))
cat("\nPREDICT O18 OR EXCESS AIR FROM OTHER TWO:\n")
print(pair_scores,n=Inf)
cat("\nSAME-WELL FIRST-ORDER DO FIT SENSITIVITY:\n")
print(DO_results,n=Inf)
cat(
  "\nLIMITATIONS:\n",
  "* Na/Cl mass ratios mix NONLINEARLY; never average endpoint ratios\n",
  "  or apply a linear ratio mixing formula.\n",
  "* Na/Cl has SHARED chloride measurement uncertainty; do not include\n",
  "  ratio and Cl as independent likelihood observations.\n",
  "* Na can undergo ion exchange; Cl can receive extra inputs.\n",
  "* Ratios near the LSR endpoint may imply unstable raw fractions.\n",
  "* 11/42 observed ratios outside endpoint range in inspected\n",
  "  prepared input: this may differ with the local dataset version.\n",
  "* Original QC uses O18, Cl and air; the DO table is a SAME-WELL\n",
  "  sensitivity test, not independent assessment of ratio selection.\n",
  "* Rainfall isotope signature can vary seasonally and spatially;\n",
  "  one fixed O18 endmember is a substantial assumption.\n",
  "* Inferred age-dependent k is not uniquely saturated-aquifer loss.\n",
  sep=""
)
