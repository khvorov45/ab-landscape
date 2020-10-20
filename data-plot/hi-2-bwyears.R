source("data-plot/data-plot.R")

# Data
data <- read_data("hi-obj2") %>%
  filter(timepoint_global %in% c(3, 4, 6, 7)) %>%
  mutate(
    study_year_lab = if_else(
      timepoint_global %in% c(3, 4), "Year 1 2016/17", "Year 2 2017/18"
    ),
    timepoint = recode(
      timepoint,
      "Pre-vax" = "Pre-vax (next year)",
      "Post-season" = "Post-season (this year)"
    ) %>%
      fct_relevel("Post-season (this year)"),
    vaccine_strain = virus == "Hkong/4801/14e"
  )
data_mod <- x_positions_by_year(data)
data_mod_alt <- x_positions_clade_year(data)

# Plots

# Individual plots with a simple year-based x-axis
indiv_plots <- plots_by_pid(data_mod, sex, age_lab, n5y_prior_vacc_lab) %>%
  map(study_year_lab_facets)
# Individual plots with a clade-then-year-based x-axis
indiv_plots_alt <- plots_by_pid(
  data_mod_alt, sex, age_lab, n5y_prior_vacc_lab
) %>%
  map(study_year_lab_facets)
# Contour plots
timepoint_year_grid <- function(plot) {
  plot +
    facet_grid(study_year_lab ~ timepoint)
}
contour_plots <- data %>%
  plots_by_pid(sex, age_lab, n5y_prior_vacc_lab, .plot_fun = plot_contour) %>%
  map(timepoint_year_grid)

# Save
save_plots(indiv_plots, "indiv-hi-2-bwyears", 42, 20)
save_plots(indiv_plots_alt, "indiv-hi-2-bwyears-alt", 42, 20)
save_plots(contour_plots, "indiv-hi-2-bwyears-contour", 18, 14)
