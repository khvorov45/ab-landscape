library(tidyverse)

# Functions ===================================================================

read_data <- function(name) {
  read_csv(glue::glue("cluster/{name}.csv"), col_types = cols())
}

plot_trace <- function(data) {
  data %>%
    ggplot(aes(iteration, estimate, col = as.factor(chain))) +
    ggdark::dark_theme_bw(verbose = FALSE) +
    theme(
      legend.position = "bottom",
      legend.box.spacing = unit(0, "null"),
      panel.spacing = unit(0, "null")
    ) +
    scale_color_discrete("Chain") +
    facet_grid(n_clusters + parameter ~ virus, scales = "free_y") +
    geom_line(alpha = 0.5)
}

save_plot <- function(plot, name, ...) {
  ggdark::ggsave_dark(glue::glue("cluster/{name}.pdf"), units = "cm", ...)
}

# Script ======================================================================

parameters <- read_data("parameters")

trace_plots <- parameters %>%
  # Each virus has one mean and sd estimate but some have numbers appended
  # to those for no reason
  mutate(parameter = str_replace(parameter, "\\.\\d{2}$", "")) %>%
  plot_trace()

save_plot(trace_plots, "trace", width = 80, height = 50)
