library(tidyverse)
library(mgcv)

# Functions ===================================================================

source("data/read_data.R")

# Script ======================================================================

cdc_virus <- read_data("cdc-virus")
cdc_vaccine <- read_data("cdc-vaccine") %>% inner_join(cdc_virus, "virus_full")

cdc_obj1_hi <- read_data("cdc-obj1-hi") %>%
  inner_join(cdc_virus, "virus_full") %>%
  mutate(
    logtitre = log(titre),
    pid_factor = as.factor(pid),
  )

fit <- gam(
  logtitre ~ s(virus_year, by = timepoint) + s(pid_factor, bs = "re") + timepoint,
  data = cdc_obj1_hi,
  method = "REML"
)

predict_data <- expand.grid(
  virus_year = unique(cdc_obj1_hi$virus_year),
  timepoint = unique(cdc_obj1_hi$timepoint)
) %>%
  as_tibble()

predicted_loghi_values <- predict(
  fit,
  predict_data,
  exclude = "s(pid_factor)",
  se = TRUE,
  newdata.guaranteed = TRUE
)

predicted_hi <- predict_data %>%
  mutate(
    predicted_loghi = predicted_loghi_values$fit,
    se = predicted_loghi_values$se.fit,
    predicted_hi = exp(predicted_loghi),
    low = exp(predicted_loghi - qnorm(0.975) * se),
    high = exp(predicted_loghi + qnorm(0.975) * se),
  )

prediction_plot <- predicted_hi %>%
  ggplot(aes(virus_year, predicted_hi, col = timepoint, fill = timepoint)) +
  ggdark::dark_theme_bw(verbose = FALSE) +
  theme(
    strip.background = element_blank(),
    panel.spacing = unit(0, "null"),
    legend.position = "bottom",
    legend.box.spacing = unit(0, "null")
  ) +
  scale_y_log10("Titre", breaks = 5 * 2^(0:15)) +
  scale_x_continuous("Virus year", expand = expansion()) +
  scale_color_manual("Timepoint", values = c("#6e6e6e", "#DF4327", "#EBA85F")) +
  scale_fill_manual("Timepoint", values = c("#6e6e6e", "#DF4327", "#EBA85F")) +
  geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.5, col = NA) +
  geom_line() +
  geom_point() +
  geom_vline(
    aes(xintercept = virus_year),
    data = cdc_vaccine,
    alpha = 0.5,
    lty = "11"
  )

ggdark::ggsave_dark(
  "gam/cdc-obj1-gam.pdf", prediction_plot,
  width = 13, height = 9, units = "cm"
)
