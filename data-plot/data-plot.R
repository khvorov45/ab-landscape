cat("Plot the data")

library(tidyverse)

data_dir <- "data"
data_plot_dir <- "data-plot"

# Functions ===================================================================

source(file.path(data_dir, "read_data.R"))

plot_one_pid <- function(data, key) {
  virus_names <- data %>%
    group_by(virus) %>%
    summarise(year_shifted = unique(year_shifted), .groups = "drop")
  plot <- data %>%
    ggplot(
      aes(year_shifted, titre,
        col = as.factor(timepoint),
        shape = as.factor(timepoint),
        linetype = as.factor(timepoint)
      )
    ) +
    ggdark::dark_theme_bw(verbose = FALSE) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      legend.position = "bottom",
      legend.box.spacing = unit(0, "null"),
      panel.grid.minor = element_blank()
    ) +
    scale_color_discrete("Timepoint") +
    scale_shape_discrete("Timepoint") +
    scale_linetype_discrete("Timepoint") +
    scale_y_log10("Titre", breaks = 5 * 2^(0:10)) +
    scale_x_continuous(
      "Virus",
      breaks = virus_names$year_shifted, labels = virus_names$virus
    ) +
    geom_line() +
    geom_point()
  attr(plot, "name") <- paste(key$pid)
  plot
}

save_pdf <- function(plot, name, dir = ".", width = 20, height = 15) {
  plotdir <- file.path(data_plot_dir, dir)
  if (!dir.exists(plotdir)) dir.create(plotdir)
  ggdark::ggsave_dark(
    file.path(data_plot_dir, dir, paste0(name, ".pdf")), plot,
    width = width,
    height = height,
    units = "cm"
  )
}

# Script ======================================================================

hi <- read_data("hi") %>%
  # Each pid at each timepoint should have one virus label
  group_by(pid, timepoint, virus_year) %>%
  mutate(year_shifted = virus_year + ((row_number() - 1) / n())) %>%
  ungroup()

indiv_hi_plots <- hi %>%
  group_by(pid) %>%
  group_map(plot_one_pid)

walk(
  indiv_hi_plots, ~ save_pdf(.x, attr(.x, "name"), "indiv-hi", 45, 15)
)
