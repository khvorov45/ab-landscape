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

# Categorise vaccine response
create_response_categories <- function(resp, tol) {
  diffs <- na.omit(resp - lag(resp))
  diffs[abs(diffs) < tol] <- 0
  if (all(diffs == 0)) {
    return("No change")
  }
  if (any(diffs > 0) & all(diffs >= 0)) {
    return("Increasing")
  }
  if (any(diffs < 0) & all(diffs <= 0)) {
    return("Decreasing")
  }
  "Alternating"
}

plot_response <- function(data) {
  data %>%
    ggplot(aes(study_year, exp(vac_resp), color = pid)) +
    ggdark::dark_theme_bw(verbose = FALSE) +
    theme(legend.position = "none") +
    scale_x_continuous("Study year", breaks = 1:3) +
    scale_y_log10(
      "Vaccine response (av. HI t2 / t1 ratio)",
      breaks = c(0.5, 1, 2, 3, 4, 5)
    ) +
    geom_line(alpha = 0.4) +
    geom_point(alpha = 0.4)
}

plot_categories <- function(data, var_name, var_lab) {
  data %>%
    ggplot(aes(!!rlang::sym(var_name))) +
    ggdark::dark_theme_bw(verbose = FALSE) +
    theme(
      strip.background = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    scale_y_continuous(
      "Count",
      expand = expansion(mult = c(0, .1)), breaks = 1:15
    ) +
    scale_x_continuous(var_lab) +
    facet_wrap(~resp_type, nrow = 2) +
    geom_histogram(binwidth = 1, col = "black")
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
  group_by(pid, study_year, age, n5y_prior_vacc, study_year_lab) %>%
  summarise(
    postvax = average_loghimid[timepoint_num == 2],
    prevax = average_loghimid[timepoint_num == 1],
    vac_resp = postvax - prevax,
    .groups = "drop"
  ) %>%
  ungroup() %>%
  mutate(study_year_cent = study_year - 1)

save_data(vac_resp, "vac-resp")

# Plot vaccine response
vac_resp_plot <- plot_response(vac_resp)

save_plot(vac_resp_plot, "vac-resp", width = 10, height = 10)

# Scatter for pre/post vax titres
vac_resp_scatter <- vac_resp %>%
  ggplot(aes(exp(prevax), exp(postvax))) +
  ggdark::dark_theme_bw(verbose = FALSE) +
  theme(
    strip.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0, "null")
  ) +
  scale_x_log10("Pre-vaccination av. titre", breaks = 5 * 2^(0:10)) +
  scale_y_log10("Post-vaccination av. titre", breaks = 5 * 2^(0:10)) +
  facet_wrap(~study_year_lab) +
  geom_point()

save_plot(vac_resp_scatter, "vac-resp-scatter", width = 15, height = 10)

lme4::lmer(
  postvax ~ prevax + study_year_cent + (1 | pid),
  vac_resp
) %>%
  summary()

# Categorise response

vac_resp_cat <- vac_resp %>%
  group_by(pid, age, n5y_prior_vacc) %>%
  summarise(
    resp_type = create_response_categories(vac_resp, tol = log(1.2)),
    .groups = "drop"
  )

vac_resp_cat %>% count(resp_type)

# Facet response by categories

vac_resp_cat_plot <- vac_resp %>%
  inner_join(select(vac_resp_cat, pid, resp_type), by = "pid") %>%
  plot_response() +
  facet_wrap(~resp_type) +
  theme(strip.background = element_blank())

save_plot(vac_resp_cat_plot, "vac-resp-cat", width = 10, height = 10)

# Plot age for the categories
vac_resp_age_plot <- plot_categories(vac_resp_cat, "age", "Age")

save_plot(vac_resp_age_plot, "vac-resp-age", width = 10, height = 10)

# Plot number of previous vaccinations
vac_resp_vacc_plot <- plot_categories(
  vac_resp_cat, "n5y_prior_vacc", "Vaccinations in past 5 years"
)

save_plot(vac_resp_vacc_plot, "vac-resp-vacc", width = 10, height = 10)
