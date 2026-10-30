# Three-tracer recharge fractions for the JOINT oxygen-loss workflow.
# This is NOT the historical river-rich, 20-percentage-point screen.
library(tidyverse)
setwd(here::here("Oxygen-Reduction-Joint-Model"))

tracers <- c("O18", "Cl", "xs_Air")
maximum_tracer_fraction_range <- 0.35

dat <- read_csv(
  "data/groundwater.csv", show_col_types = FALSE,
  col_types = cols(Site_ID = col_character(), site = col_character(),
                   .default = col_guess())
) %>%
  mutate(sample_id = row_number(), .before = 1)

# Tracer end members are groundwater samples, not independently measured
# pure recharge. In particular, 2801 is set to 100% river BY DEFINITION.
river <- dat %>% filter(Site_ID == "2801")
rain_candidates <- dat %>% filter(
  Site_ID %in% c("1032", "4294") |
    (Site_ID == "5988 HB Prison" & Depth == 108)
)
stopifnot(
  nrow(river) == 1, nrow(rain_candidates) == 3,
  all(complete.cases(rain_candidates[, tracers]))
)

river_endpoint <- vapply(tracers, function(n) river[[n]][1], numeric(1))
rain_endpoint <- vapply(
  tracers, function(n) mean(rain_candidates[[n]]), numeric(1)
)
stopifnot(all(abs(rain_endpoint - river_endpoint) > 1e-9))

estimate_mixing <- function(d, rain_end) {
  raw <- sapply(tracers, function(n) {
    (d[[n]] - river_endpoint[[n]]) /
      (rain_end[[n]] - river_endpoint[[n]])
  })
  clipped <- raw
  clipped[] <- pmax(0, pmin(1, as.vector(raw)))
  n_available <- rowSums(is.finite(raw))

  median_fraction <- apply(clipped, 1, function(v) {
    v <- v[is.finite(v)]
    if (length(v)) median(v) else NA_real_
  })
  disagreement <- apply(raw, 1, function(v) {
    v <- v[is.finite(v)]
    if (length(v) >= 2) diff(range(v)) else NA_real_
  })

  d %>% mutate(
    rain_fraction_18O = clipped[, "O18"],
    rain_fraction_Cl = clipped[, "Cl"],
    rain_fraction_excess_air = clipped[, "xs_Air"],
    raw_rain_fraction_18O = raw[, "O18"],
    raw_rain_fraction_Cl = raw[, "Cl"],
    raw_rain_fraction_excess_air = raw[, "xs_Air"],
    rainfall_contribution = median_fraction,
    river_contribution = 1 - median_fraction,
    tracers_available = n_available,
    mixing_disagreement = disagreement,
    mixing_reliable = n_available == 3 &
      is.finite(disagreement) &
      disagreement <= maximum_tracer_fraction_range
  )
}

mixing <- estimate_mixing(dat, rain_endpoint)
write_csv(mixing, "outputs/three_tracer_mixing.csv")

endmembers <- tibble(
  tracer = tracers,
  river_endpoint = as.numeric(river_endpoint),
  rainfall_endpoint = as.numeric(rain_endpoint),
  river_site = "2801",
  rainfall_sites = "1032; 4294; 5988 HB Prison (108 m)"
)
write_csv(endmembers, "outputs/three_tracer_endmembers.csv")

# Omit one land-surface end member at a time as a selection sensitivity test.
endpoint_sensitivity <- map_dfr(seq_len(nrow(rain_candidates)), function(i) {
  revised <- vapply(tracers, function(n) {
    mean(rain_candidates[[n]][-i])
  }, numeric(1))
  test <- estimate_mixing(dat, revised)
  tibble(
    omitted_endpoint = paste0(rain_candidates$Site_ID[i],
                              " (", rain_candidates$Depth[i], " m)"),
    n_passing = sum(test$mixing_reliable),
    median_absolute_fraction_change = median(
      abs(test$rainfall_contribution - mixing$rainfall_contribution),
      na.rm = TRUE
    ),
    max_absolute_fraction_change = max(
      abs(test$rainfall_contribution - mixing$rainfall_contribution),
      na.rm = TRUE
    )
  )
})
write_csv(endpoint_sensitivity, "outputs/rainfall_endpoint_sensitivity.csv")

production <- mixing %>% filter(
  mixing_reliable, is.finite(MRT), MRT >= 0, is.finite(DO), DO > 0
)
stopifnot(
  nrow(dat) == 42,
  sum(mixing$mixing_reliable) == 29,
  nrow(production) == 29,
  sum(production$rainfall_contribution >= 0.5) == 5
)
cat("Mixing: 42 records; 29 production wells; 5 rainfall-dominated.\n")
print(endpoint_sensitivity)
