source("data-plot/data-plot.R")

# Data
data <- read_data("hi-rmh-hcw")
data_mod <- x_positions_by_year(data)
data_mod_alt <- x_positions_clade_year(data)

# Plots

# Individual plots with a simple year-based x-axis
indiv_plots <- plots_by_pid(data_mod, group, sex, age_lab)
# Individual plots with a clade-then-year-based x-axis
indiv_plots_alt <- plots_by_pid(data_mod_alt, group, sex, age_lab)

# Save
save_plots(indiv_plots, "indiv-hi-rmh-hcw", 42, 15)
save_plots(indiv_plots_alt, "indiv-hi-rmh-hcw-alt", 42, 15)
