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

fit_gam <- function(data, key) {
  fit <- gam(
    logtitre ~ s(virus_year, by = timepoint_x_prior_vacs_3cat) +
      s(pid_factor, bs = "re") +
      timepoint * prior_vacs_3cat,
    data = data,
    method = "REML"
  )
  attr(fit, "key") <- key
  fit
}

cdc_obj1_grouped <- cdc_obj1_hi %>% group_by(site, study_year)
future::plan(future::multisession)
fits <- furrr::future_map2(
  group_split(cdc_obj1_grouped),
  group_keys(cdc_obj1_grouped) %>% group_split(row_number(), .keep = FALSE),
  fit_gam
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

predict_values <- function(fit) {
  predicted_values <- predict(
    fit,
    predict_data,
    exclude = "s(pid_factor)",
    se = TRUE,
    newdata.guaranteed = TRUE
  )
  attr(predicted_values, "key") <- attr(fit, "key")
  predicted_values
}

predicted_loghi_values <- map(fits, predict_values)

gen_predicted_hi <- function(predicted_values) {
  predict_data %>%
    mutate(
      predicted_loghi = predicted_values$fit,
      se = predicted_values$se.fit,
      predicted_hi = exp(predicted_loghi),
      low = exp(predicted_loghi - qnorm(0.975) * se),
      high = exp(predicted_loghi + qnorm(0.975) * se),
    ) %>%
    bind_cols(attr(predicted_values, "key"))
}

predicted_hi <- map_dfr(predicted_loghi_values, gen_predicted_hi)

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
    coord_cartesian(ylim = c(5, max(data$predicted_hi))) +
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
  facet_grid(
    study_year + prior_vacs_3cat ~ site,
    labeller = function(x) {
      if ("site" %in% names(x)) {
        return(x)
      }
      x %>%
        mutate(
          study_year = paste0("Year ", study_year),
          prior_vacs_3cat = paste0(prior_vacs_3cat, " prior")
        )
    }
  ) +
  scale_color_discrete("Timepoint") +
  scale_fill_discrete("Timepoint")

prediction_plot_prior_vacs <- plot_predictions(predicted_hi, prior_vacs_3cat) +
  facet_grid(
    study_year + timepoint ~ site,
    labeller = function(x) {
      if ("site" %in% names(x)) {
        return(x)
      }
      x %>%
        mutate(
          study_year = paste0("Year ", study_year),
          timepoint = as.character(timepoint)
        )
    }
  ) +
  scale_color_discrete("Prior vaccinations") +
  scale_fill_discrete("Prior vaccinations")

save_plot <- function(plot, name) {
  ggdark::ggsave_dark(
    glue::glue("gam/{name}.pdf"), plot,
    width = 13, height = 23, units = "cm"
  )
}

save_plot(prediction_plot_timepoint, "cdc-obj1-gam-timepoint")
save_plot(prediction_plot_prior_vacs, "cdc-obj1-gam-prior-vacs")
