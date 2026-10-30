
library(tidyverse)

# ================================================================
# 1. SETTINGS
# ================================================================

# Fonts and text sizes (points)
cfg <- list(
  font_family       = "sans",
  base_font_size    = 14,
  title_size        = 17,
  subtitle_size     = 12,
  axis_title_size   = 13,
  axis_text_size    = 11,
  legend_title_size = 12,
  legend_text_size  = 11,
  facet_title_size  = 12,
  caption_size      = 10,

  # Graphics
  point_size        = 4,
  open_point_size   = 4,
  line_width        = 0.9,
  error_line_width  = 0.7,
  label_size        = 4,
  label_threshold   = 2,       # Label sites with |DO error| > 2
  legend_position   = "bottom",

  # Export
  dpi               = 300,
  figure_width      = 12,
  figure_height     = 7,
  wide_width        = 12,
  tall_height       = 11
)

# All colours
colours <- list(
  # Continuous rainfall fraction: river-rich -> rainfall-rich
  rainfall_low      = "#2166AC",
  rainfall_high     = "#B2182B",

  # Individual plot elements
  observed          = "#2166AC",
  predicted         = "#D6604D",
  reference_line    = "#777777",
  residual_line     = "#555555",
  connector         = "#AAAAAA",
  text              = "#222222",
  grid              = "#E5E5E5",
  facet_background  = "#F2F2F2",
  background        = "#FFFFFF",

  # Components of the oxygen budget
  budget = c(
    "Remaining DO"            = "#4C78A8",
    "Aquifer reduction"       = "#E45756",
    "Recharge-associated loss" = "#F2CF5B"
  ),

  # Individual tracers
  tracers = c(
    "Oxygen-18" = "#4477AA",
    "Chloride"  = "#CC6677",
    "Excess air" = "#228833"
  )
)

# ================================================================
# 2. READ DATA
# ================================================================

required_files <- c(
  "outputs/three_tracer_mixing.csv",
  "outputs/joint_model_fits.csv",
  "outputs/joint_model_LOOWO_predictions.csv"
)

missing <- required_files[!file.exists(required_files)]

if (length(missing) > 0) {
  stop(
    "Missing files: ",
    paste(missing, collapse = ", "),
    "\nRun scripts 01 and 02 first."
  )
}

dir.create(
  "plots",
  showWarnings = FALSE,
  recursive = TRUE
)

mixing <- read_csv(
  "outputs/three_tracer_mixing.csv",
  show_col_types = FALSE,
  col_types = cols(
    Site_ID = col_character(),
    .default = col_guess()
  )
)

fits <- read_csv(
  "outputs/joint_model_fits.csv",
  show_col_types = FALSE
)

predictions <- read_csv(
  "outputs/joint_model_LOOWO_predictions.csv",
  show_col_types = FALSE,
  col_types = cols(
    Site_ID = col_character(),
    .default = col_guess()
  )
)

model_names <- c(
  aquifer_only = "Aquifer only",
  infiltration_only = "Recharge loss only",
  combined = "Combined"
)

# ================================================================
# 3. GLOBAL PLOT THEME
# ================================================================

plot_theme <- theme_bw(
  base_size = cfg$base_font_size,
  base_family = cfg$font_family
) +
  theme(
    plot.title = element_text(
      size = cfg$title_size,
      face = "bold",
      colour = colours$text,
      margin = margin(b = 8)
    ),
    plot.subtitle = element_text(
      size = cfg$subtitle_size,
      colour = colours$text,
      margin = margin(b = 12)
    ),
    axis.title = element_text(
      size = cfg$axis_title_size,
      colour = colours$text
    ),
    axis.text = element_text(
      size = cfg$axis_text_size,
      colour = colours$text
    ),
    legend.title = element_text(
      size = cfg$legend_title_size
    ),
    legend.text = element_text(
      size = cfg$legend_text_size
    ),
    strip.text = element_text(
      size = cfg$facet_title_size,
      face = "bold"
    ),
    strip.background = element_rect(
      fill = colours$facet_background
    ),
    plot.caption = element_text(
      size = cfg$caption_size,
      hjust = 0,
      lineheight = 1.15,
      margin = margin(t = 12)
    ),
    panel.grid.major = element_line(
      colour = colours$grid,
      linewidth = 0.3
    ),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(
      fill = colours$background,
      colour = NA
    ),
    legend.position = cfg$legend_position,
    plot.margin = margin(12, 15, 12, 12)
  )

theme_set(plot_theme)

rainfall_scale <- function() {
  scale_colour_gradient(
    low = colours$rainfall_low,
    high = colours$rainfall_high,
    limits = c(0, 1),
    labels = scales::percent,
    name = "Rainfall contribution"
  )
}

# Turbo is used for the first-order residence-time plot so the two recharge
# end members read directly as blue (river) and red (rainfall).
rainfall_turbo_scale <- function() {
  scale_colour_viridis_c(
    option = "turbo",
    begin = 0.1,
    end = 0.9,
    direction = 1,
    limits = c(0, 1),
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    labels = c("River", "25%", "50%", "75%", "Rainfall"),
    name = "Rainfall contribution",
    guide = guide_colourbar(barwidth = grid::unit(10, "cm"))
  )
}

save_plot <- function(
  p, filename,
  width = cfg$figure_width,
  height = cfg$figure_height
) {
  print(p)

  ggsave(
    filename = file.path("plots", filename),
    plot = p,
    width = width,
    height = height,
    dpi = cfg$dpi,
    bg = colours$background
  )
}

# ================================================================
# 4. OBSERVED VERSUS PREDICTED DO
# ================================================================
# All predictions are leave-one-well-out.

plot_data <- predictions %>%
  mutate(
    model = factor(
      model,
      levels = names(model_names),
      labels = model_names
    ),
    kinetic = factor(
      kinetic,
      levels = c("First order", "Zero order")
    )
  )

max_DO <- ceiling(max(
  plot_data$DO_observed,
  plot_data$DO_predicted,
  na.rm = TRUE
))

p1 <- ggplot(
  plot_data,
  aes(DO_observed, DO_predicted)
) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    colour = colours$reference_line,
    linewidth = cfg$line_width
  ) +
  geom_point(
    aes(colour = rainfall_fraction),
    size = cfg$point_size
  ) +
  facet_grid(kinetic ~ model) +
  coord_equal(
    xlim = c(0, max_DO),
    ylim = c(0, max_DO)
  ) +
  rainfall_scale() +
  labs(
    title = "Observed versus predicted groundwater DO",
    subtitle = "Leave-one-well-out predictions",
    x = "Observed DO (mg/L)",
    y = "Predicted DO (mg/L)"
  )

save_plot(
  p1,
  "01_observed_vs_predicted.png",
  cfg$wide_width,
  8
)

# ================================================================
# 5. DO VERSUS MEAN RESIDENCE TIME
# ================================================================
# Combined first-order model.
# Open circles = independently predicted DO.
# Filled circles = measured DO.
# Connecting lines = prediction error.

age_data <- predictions %>%
  filter(
    kinetic == "First order",
    model == "combined"
  ) %>%
  mutate(
    error = DO_observed - DO_predicted,
    observed_rate = DO_observed / DO_initial,
    predicted_rate = DO_predicted / DO_initial
  )

first_order_fit <- fits %>%
  filter(
    kinetic == "First order",
    model == "combined"
  )

stopifnot(nrow(first_order_fit) == 1)

first_order_rate <- first_order_fit$aquifer_rate[[1]]
best_fit_curve <- tibble(
  MRT = seq(0, max(age_data$MRT), length.out = 300),
  fitted_rate = exp(-first_order_rate * MRT)
)
equation_label <- paste0(
  "Best fit: DOsample / DOinitial = exp(-",
  format(round(first_order_rate, 3), nsmall = 3),
  " × MRT)"
)

p2 <- ggplot(age_data, aes(x = MRT)) +
  geom_segment(
    aes(
      xend = MRT,
      y = observed_rate,
      yend = predicted_rate
    ),
    colour = colours$connector,
    linewidth = cfg$error_line_width
  ) +
  geom_line(
    data = best_fit_curve,
    aes(x = MRT, y = fitted_rate),
    colour = "black",
    linewidth = 1.1,
    inherit.aes = FALSE
  ) +
  geom_point(
    aes(
      y = observed_rate,
      colour = rainfall_fraction
    ),
    size = cfg$point_size
  ) +
  geom_point(
    aes(
      y = predicted_rate,
      colour = rainfall_fraction
    ),
    shape = 1,
    size = cfg$open_point_size,
    stroke = cfg$line_width
  ) +
  geom_text(
    data = ~filter(
      .x,
      abs(error) > cfg$label_threshold
    ),
    aes(
      y = observed_rate,
      label = site
    ),
    size = cfg$label_size,
    colour = colours$text,
    vjust = -0.8,
    check_overlap = TRUE
  ) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = equation_label,
    hjust = 1.05,
    vjust = 1.5,
    size = cfg$label_size,
    colour = colours$text
  ) +
  rainfall_turbo_scale() +
  scale_x_continuous(limits = c(0, NA)) +
  scale_y_log10() +
  labs(
    title = "Groundwater DO rate versus residence time",
    subtitle = "Combined first-order model",
    x = "Mean residence time (years)",
    y = expression(DO[sample] / DO[initial]),
    caption = str_wrap(
      paste(
        "Filled circles: observed DO rate; open circles:",
        "leave-one-well-out predictions.",
        "Vertical lines show prediction errors.",
        "The black curve is the full-data combined first-order best fit.",
        "Initial DO varies with recharge-source proportions;",
        "the y axis is logarithmic."
      ),
      width = 120
    )
  )

save_plot(
  p2,
  "02_DO_vs_residence_time.png"
)

# ================================================================
# 6. RESIDUAL DIAGNOSTICS
# ================================================================
# First-order aquifer-only and combined models.
# Each prediction excludes the observed well from calibration.

residual_base <- predictions %>%
  filter(
    kinetic == "First order",
    model %in% c("aquifer_only", "combined")
  ) %>%
  mutate(
    model = factor(
      model,
      levels = c("aquifer_only", "combined"),
      labels = c("Aquifer only", "Combined")
    )
  )

# Preserve rainfall_fraction for colouring both panels.
residual_data <- bind_rows(
  residual_base %>%
    transmute(
      model,
      rainfall_fraction,
      residual,
      predictor = "Residence time (years)",
      x = MRT
    ),

  residual_base %>%
    transmute(
      model,
      rainfall_fraction,
      residual,
      predictor = "Rainfall contribution",
      x = rainfall_fraction
    )
) %>%
  mutate(
    predictor = factor(
      predictor,
      levels = c(
        "Rainfall contribution",
        "Residence time (years)"
      )
    )
  )

p3 <- ggplot(
  residual_data,
  aes(x, residual)
) +
  geom_hline(
    yintercept = 0,
    colour = colours$reference_line,
    linetype = "dashed",
    linewidth = cfg$error_line_width
  ) +
  geom_point(
    aes(colour = rainfall_fraction),
    size = cfg$point_size
  ) +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = FALSE,
    colour = colours$residual_line,
    linewidth = cfg$line_width
  ) +
  facet_grid(
    model ~ predictor,
    scales = "free_x"
  ) +
  rainfall_scale() +
  labs(
    title = "First-order model residuals",
    subtitle = "Leave-one-well-out prediction errors",
    x = NULL,
    y = "Observed minus predicted DO (mg/L)",
    caption = paste(
      "Linear trends are descriptive diagnostics.",
      "They do not establish the location or mechanism",
      "of oxygen reduction."
    )
  )

save_plot(
  p3,
  "03_residual_diagnostics.png",
  cfg$wide_width,
  8
)

# ================================================================
# 7. OXYGEN LOSS PARTITION
# ================================================================
# Based on the combined first-order model.
#
# IMPORTANT: This is a MODELLED partition.
# It does not establish that recharge-associated consumption
# physically occurred in the unsaturated zone.

fit <- fits %>%
  filter(
    kinetic == "First order",
    model == "combined"
  )

stopifnot(nrow(fit) == 1)

river_DO <- 7.5
reference_DO <- fit$reference_rainfall_DO[[1]]
effective_DO <- fit$effective_rainfall_DO_mg_L[[1]]
k <- fit$aquifer_rate[[1]]

budget <- mixing %>%
  filter(
    mixing_reliable,
    is.finite(DO),
    DO > 0,
    is.finite(MRT),
    MRT >= 0
  ) %>%
  mutate(
    fR = river_contribution,
    fP = rainfall_contribution,

    initial_DO =
      fR * river_DO + fP * reference_DO,

    recharge_loss =
      fP * (reference_DO - effective_DO),

    DO_at_aquifer_entry =
      fR * river_DO + fP * effective_DO,

    aquifer_loss =
      DO_at_aquifer_entry *
      (1 - exp(-k * MRT)),

    remaining_DO =
      DO_at_aquifer_entry * exp(-k * MRT),

    site_label = paste0(
      site, " [", Site_ID, "]"
    )
  )

# Verify that the budget closes.
stopifnot(all(
  abs(
    budget$recharge_loss +
    budget$aquifer_loss +
    budget$remaining_DO -
    budget$initial_DO
  ) < 1e-7
))

budget_long <- budget %>%
  mutate(
    site_label = forcats::fct_reorder(
      site_label, MRT
    )
  ) %>%
  select(
    site_label,
    recharge_loss,
    aquifer_loss,
    remaining_DO
  ) %>%
  pivot_longer(
    cols = c(
      recharge_loss,
      aquifer_loss,
      remaining_DO
    ),
    names_to = "component",
    values_to = "oxygen_mg_L"
  ) %>%
  mutate(
    component = factor(
      component,
      levels = c(
        "remaining_DO",
        "aquifer_loss",
        "recharge_loss"
      ),
      labels = c(
        "Remaining DO",
        "Aquifer reduction",
        "Recharge-associated loss"
      )
    )
  )

p4 <- ggplot(
  budget_long,
  aes(
    site_label,
    oxygen_mg_L,
    fill = component
  )
) +
  geom_col(width = 0.8) +
  coord_flip() +
  scale_fill_manual(
    values = colours$budget
  ) +
  labs(
    title = "Modelled oxygen budget by groundwater site",
    subtitle = paste0(
      "Combined first-order model; ",
      "k = ", round(k, 3), " year^-1"
    ),
    x = NULL,
    y = "Oxygen concentration (mg/L)",
    fill = "Oxygen budget",
    caption = paste(
      "Bars represent the source-weighted initial DO,",
      "partitioned into remaining DO and two inferred losses.",
      "Recharge-associated loss is not independently",
      "identified as infiltration-pathway consumption."
    )
  )

save_plot(
  p4,
  "04_oxygen_loss_partition.png",
  cfg$wide_width,
  cfg$tall_height
)

# ================================================================
# 8. AGREEMENT BETWEEN MIXING TRACERS
# ================================================================
# All wells passing the three-tracer screen.

tracer_data <- mixing %>%
  filter(mixing_reliable) %>%
  select(
    Site_ID,
    rainfall_contribution,
    rain_fraction_18O,
    rain_fraction_Cl,
    rain_fraction_excess_air
  ) %>%
  pivot_longer(
    cols = starts_with("rain_fraction_"),
    names_to = "tracer",
    values_to = "tracer_fraction"
  ) %>%
  mutate(
    tracer = recode(
      tracer,
      rain_fraction_18O = "Oxygen-18",
      rain_fraction_Cl = "Chloride",
      rain_fraction_excess_air = "Excess air"
    )
  )

n_reliable <- n_distinct(tracer_data$Site_ID)

p5 <- ggplot(
  tracer_data,
  aes(
    rainfall_contribution,
    tracer_fraction
  )
) +
  geom_abline(
    slope = 1,
    intercept = 0,
    colour = colours$reference_line,
    linetype = "dashed",
    linewidth = cfg$line_width
  ) +
  geom_point(
    aes(colour = tracer),
    size = cfg$point_size
  ) +
  facet_wrap(~tracer, nrow = 1) +
  scale_colour_manual(
    values = colours$tracers,
    guide = "none"
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    labels = scales::percent
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    labels = scales::percent
  ) +
  coord_equal() +
  labs(
    title = "Agreement between recharge mixing tracers",
    subtitle = paste(
      n_reliable,
      "wells passing the three-tracer screen"
    ),
    x = "Combined rainfall contribution",
    y = "Individual tracer rainfall fraction",
    caption = paste(
      "Dashed line: 1:1 agreement.",
      "The combined fraction is the median of the",
      "three tracer estimates.",
      "Excess-air conservativeness remains an assumption."
    )
  )

save_plot(
  p5,
  "05_three_tracer_agreement.png",
  cfg$wide_width,
  5.5
)

# ================================================================
# 9. FINISHED
# ================================================================

cat(
  "\nFive figures displayed and saved in plots/.\n",
  "Change cfg or colours at the top to alter their appearance.\n"
)
