cat("Plot the data")

library(tidyverse)

data_dir <- "data"
data_plot_dir <- "data-plot"

# Functions ===================================================================

source(file.path(data_dir, "read_data.R"))

plot_one_pid <- function(data, key) {
  virus_names <- data %>%
    group_by(virus) %>%
    summarise(x_position = unique(x_position), .groups = "drop")
  plot <- data %>%
    ggplot(
      aes(x_position, titre,
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
      breaks = virus_names$x_position, labels = virus_names$virus,
      expand = c(0, 0.5)
    ) +
    coord_cartesian(
      clip = "off", ylim = c(5, max(data$titre)),
      xlim = c(min(data$x_position), max(data$x_position))
    ) +
    geom_line() +
    geom_point() +
    labs(caption = paste(key$pid, key$group, key$sex, key$age, "years"))
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
  # Each pid should have one virus label per year
  group_by(pid, virus_year) %>%
  mutate(
    x_position =
      virus_year + (map_int(virus, ~ which(.x == unique(virus))) - 1) /
        length(unique(virus))
  ) %>%
  ungroup()

indiv_hi_plots <- hi %>%
  filter(pid == "HIA15611") %>%
  group_by(pid, group, sex, age) %>%
  group_map(plot_one_pid)

walk(
  indiv_hi_plots, ~ save_pdf(.x, attr(.x, "name"), "indiv-hi", 45, 15)
)
