library(tidyverse)
library(mgcv)

# Functions ===================================================================

source("data/read_data.R")

# Script ======================================================================

cdc_virus <- read_data("cdc-virus")
cdc_vaccine <- read_data("cdc-vaccine") %>% inner_join(cdc_virus, "virus_full")

cdc_obj1_hi <- read_data("cdc-obj1-hi") %>%
  inner_join(cdc_virus, "virus_full") %>%
  inner_join(read_data("cdc-obj1-participant"), "pid") %>%
  mutate(
    logtitre = log(titre),
    pid_factor = as.factor(pid),
    prior_vacs_3cat = case_when(
      prior_vacs == 0 ~ "0",
      prior_vacs == 5 ~ "5",
      TRUE ~ "1-4"
    ) %>% as.factor(),
    timepoint_x_prior_vacs_3cat = paste0(timepoint, prior_vacs_3cat) %>%
      as.factor()
  )

fit <- gam(
  logtitre ~ s(virus_year, by = timepoint_x_prior_vacs_3cat) +
    s(pid_factor, bs = "re") +
    timepoint * prior_vacs_3cat,
  data = cdc_obj1_hi,
  method = "REML"
)

predict_data <- expand.grid(
  virus_year = unique(cdc_obj1_hi$virus_year),
  timepoint = unique(cdc_obj1_hi$timepoint),
  prior_vacs_3cat = unique(cdc_obj1_hi$prior_vacs_3cat)
) %>%
  as_tibble() %>%
  mutate(
    timepoint_x_prior_vacs_3cat = paste0(timepoint, prior_vacs_3cat)
  )

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

plot_predictions <- function(data, group_var) {
  # The og colors were c("#6e6e6e", "#DF4327", "#EBA85F")
  group_var_q <- rlang::enquo(group_var)
  data %>%
    ggplot(aes(virus_year, predicted_hi, col = !!group_var_q, fill = !!group_var_q)) +
    ggdark::dark_theme_bw(verbose = FALSE) +
    theme(
      strip.background = element_blank(),
      panel.spacing = unit(0, "null"),
      legend.position = "bottom",
      legend.box.spacing = unit(0, "null")
    ) +
    scale_y_log10("Titre", breaks = 5 * 2^(0:15)) +
    scale_x_continuous("Virus year", expand = expansion()) +
    geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.5, col = NA) +
    geom_line() +
    geom_point() +
    geom_vline(
      aes(xintercept = virus_year),
      data = cdc_vaccine,
      alpha = 0.5,
      lty = "11"
    )
}

prediction_plot_timepoint <- plot_predictions(predicted_hi, timepoint) +
  facet_wrap(
    vars(prior_vacs_3cat),
    ncol = 1,
    strip.position = "right",
    labeller = as_labeller(function(x) paste(x, "prior"))
  ) +
  scale_color_discrete("Timepoint") +
  scale_fill_discrete("Timepoint")

prediction_plot_prior_vacs <- plot_predictions(predicted_hi, prior_vacs_3cat) +
  facet_wrap(vars(timepoint), ncol = 1, strip.position = "right") +
  scale_color_discrete("Prior vaccinations") +
  scale_fill_discrete("Prior vaccinations")

save_plot <- function(plot, name) {
  ggdark::ggsave_dark(
    glue::glue("gam/{name}.pdf"), plot,
    width = 13, height = 13, units = "cm"
  )
}

save_plot(prediction_plot_timepoint, "cdc-obj1-gam-timepoint")
save_plot(prediction_plot_prior_vacs, "cdc-obj1-gam-prior-vacs")
