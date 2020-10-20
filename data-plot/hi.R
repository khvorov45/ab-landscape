source("data-plot/data-plot.R")

# Data
data <- read_data("hi")
data_mod <- x_positions_by_year(data)
data_mod_alt <- x_positions_clade_year(data)

# Plots

# Individual plots with a simple year-based x-axis
indiv_plots <- plots_by_pid(data_mod, group, sex, age_lab)
# Individual plots with a clade-then-year-based x-axis
indiv_plots_alt <- plots_by_pid(data_mod_alt, group, sex, age_lab)
# Contour plots
contour_plots <- data %>%
  plots_by_pid(group, sex, age_lab, .plot_fun = plot_contour)

# Save
save_plots(indiv_plots, "indiv-hi", 42, 15)
save_plots(indiv_plots_alt, "indiv-hi-alt", 42, 15)
save_plots(contour_plots, "indiv-hi-countour", 21, 9.5)
