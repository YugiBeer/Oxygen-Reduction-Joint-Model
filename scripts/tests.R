library(ggplot2)

f <- readRDS(
  "outputs/04_river_DO_prior/2801_masked_no2130/fit_C_prior_sd1p5_order1.rds"
)

k <- as.numeric(f$draws("aquifer_rate_actual"))
q <- quantile(k, c(0.025, 0.5, 0.975))

ggplot(data.frame(k), aes(x = k)) +
  geom_density(fill = "steelblue", alpha = 0.25, linewidth = 1) +
  geom_vline(xintercept = q[2], linewidth = 0.8) +
  geom_vline(xintercept = q[c(1, 3)],
             linetype = "dashed", colour = "grey40") +
  labs(
    x = expression("First-order aquifer reaction rate (yr"^-1*")"),
    y = "Posterior density"
  ) +
  theme_classic(base_size = 14)