cat("drop-off analysis")

library(tidyverse)

drop_off_dir <- "drop-off"
data_dir <- "data"

# Functions ===================================================================

source(file.path(data_dir, "read_data.R"))

calc_auc <- function(hi, ...) {
  hi %>%
    # Remove the viruses that circulated before birth
    filter(virus_year >= year_of_birth, !is.na(logtitre_mid)) %>%
    # Make sure the is one measurement per individual per timepoint
    # per virus year
    group_by(virus_year, ...) %>%
    summarise(logtitre_mid = mean(logtitre_mid), .groups = "drop") %>%
    group_by(...) %>%
    summarise(
      area_under_curve = DescTools::AUC(virus_year, logtitre_mid),
      average_loghimid = area_under_curve / (max(virus_year) - min(virus_year)),
      .groups = "drop"
    )
}

plot_spag <- function(aucs, xdate = TRUE) {
  date_means <- aucs %>%
    group_by(timepoint_global) %>%
    summarise(timepoint_date_mean = mean(timepoint_date), .groups = "drop")
  xvar <- if (xdate) {
    rlang::sym("timepoint_date")
  } else {
    rlang::sym("timepoint_global")
  }
  plot <- aucs %>%
    rename(Timepoint = timepoint) %>%
    ggplot(
      aes(!!xvar, exp(average_loghimid),
        col = Timepoint,
        shape = Timepoint
      )
    ) +
    ggdark::dark_theme_bw(verbose = FALSE) +
    theme(
      legend.position = "bottom",
      legend.box.spacing = unit(0, "null"),
      legend.margin = margin(),
      panel.grid.minor = element_blank()
    ) +
    scale_y_log10("Average titre", breaks = 5 * 2^(0:10)) +
    geom_line(aes(group = pid), alpha = 0.4) +
    geom_point(alpha = 0.4)
  if (xdate) {
    plot <- plot +
      scale_x_date("Date", breaks = date_means$timepoint_date_mean) +
      theme(axis.text.x = element_text(angle = 30, hjust = 1))
  } else {
    plot <- plot +
      scale_x_continuous("Global timepoint", breaks = 1:9)
  }
  plot
}

save_data <- function(data, name, ...) {
  write_csv(
    data,
    file.path(drop_off_dir, paste0(name, ".csv"))
  )
}

save_plot <- function(plot, name, ...) {
  ggdark::ggsave_dark(
    file.path(drop_off_dir, paste0(name, ".png")),
    plot,
    units = "cm", ...
  )
}

# Script ======================================================================

hi_2 <- read_data("hi-obj2")

# Summarise to average titre
aucs <- hi_2 %>% calc_auc(
  pid, year_of_birth, timepoint_global, timepoint, timepoint_date,
  timepoint_date_imputed, study_year, timepoint_num,
  study_year_lab, age, n5y_prior_vacc
)

# Plots -----------------------------------------------------------------------

spag <- plot_spag(aucs)
save_plot(spag, "spag", width = 15, height = 8)

spag_facets <- spag +
  facet_wrap(
    ~study_year_lab,
    nrow = 1, scales = "free_x"
  ) +
  theme(
    strip.background = element_blank(),
    panel.spacing = unit(0, "null")
  )
save_plot(spag_facets, "spag-facets", width = 20, height = 8)

spag_labs <- plot_spag(aucs, xdate = FALSE)
save_plot(spag_labs, "spag-even", width = 15, height = 8)

# Quantify response -----------------------------------------------------------

# Get the vaccine response
vac_resp <- aucs %>%
  group_by(pid, study_year, age, n5y_prior_vacc) %>%
  summarise(
    vac_resp =
      average_loghimid[timepoint_num == 2] -
        average_loghimid[timepoint_num == 1],
    .groups = "drop"
  ) %>%
  ungroup() %>%
  mutate(study_year_cent = study_year - 1)

save_data(vac_resp, "vac-resp")

# Plot vaccine response
vac_resp_plot <- vac_resp %>%
  ggplot(aes(study_year, exp(vac_resp))) +
  ggdark::dark_theme_bw(verbose = FALSE) +
  scale_x_continuous("Study year", breaks = 1:3) +
  scale_y_log10(
    "Vaccine response (av. HI t2 / t1 ratio)",
    breaks = c(0.5, 1, 2, 3, 4, 5)
  ) +
  geom_line(aes(group = pid), alpha = 0.4) +
  geom_point(alpha = 0.4)

save_plot(vac_resp_plot, "vac-resp", width = 10, height = 10)

lme4::lmer(
  vac_resp ~ study_year_cent + n5y_prior_vacc + age + (1 | pid),
  vac_resp
) %>%
  summary()
