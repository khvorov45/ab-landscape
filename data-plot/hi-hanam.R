source("data-plot/data-plot.R")

# Data
data <- read_data("hi-hanam")
data_mod <- x_positions_by_year(data)
data_mod_alt <- x_positions_clade_year(data)

# Plots

# Individual plots with a simple year-based x-axis
indiv_plots <- plots_by_pid(data_mod, prior_h3_lab)
# Individual plots with a clade-then-year-based x-axis
indiv_plots_alt <- plots_by_pid(data_mod_alt, prior_h3_lab)
# Not enough viruses with known coordinates for the contour plots

# Save
save_plots(indiv_plots, "indiv-hi-hanam", width = 42, height = 15)
save_plots(indiv_plots_alt, "indiv-hi-hanam-alt", width = 42, height = 15)
