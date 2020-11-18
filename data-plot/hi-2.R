source("data-plot/data-plot.R")

# Data
data <- read_data("hi-obj2")
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
save_plots(indiv_plots, "indiv-hi-2", width = 42, height = 42)
save_plots(indiv_plots_alt, "indiv-hi-2-alt", width = 42, height = 42)
save_plots(contour_plots, "indiv-hi-2-contour", width = 25.5, height = 23.5)
