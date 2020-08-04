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
    labs(caption = pmap(key, ~ paste(...))) +
    geom_text(
      aes(
        label = if_else(clade == "(Missing)", "", clade),
        y = 5, x = x_position
      ),
      angle = 90, hjust = 0, col = "gray40", alpha = 0.5, size = 3,
      inherit.aes = FALSE
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

remap_range <- function(range, output_start, output_end) {
  if (length(unique(range)) == 1L) {
    return(output_start)
  }
  output_start +
    ((output_end - output_start) /
      (max(range) - min(range))) * (range - min(range))
}

x_positions_by_year <- function(hi) {
  hi %>%
    # Arrange by year
    group_by(virus_year) %>%
    mutate(
      x_position =
        virus_year + (map_int(virus, ~ which(.x == unique(virus))) - 1) /
          length(unique(virus)),
    ) %>%
    ungroup()
}

x_positions_clade_year <- function(hi) {
  hi_mod <- x_positions_by_year(hi)
  clade_postions <- gen_clade_positions(hi)
  # Work out the position of viruses within clades
  # Use the temporal x_position because it's unique for each virus
  hi_mod %>%
    group_by(clade) %>%
    mutate(
      rel_position = remap_range(x_position, 0.12, 0.88)
    ) %>%
    ungroup() %>%
    inner_join(clade_postions, by = "clade") %>%
    mutate(
      x_position = clade_start + rel_position * clade_len
    )
}

gen_clade_positions <- function(hi) {
  hi %>%
    group_by(clade) %>%
    summarise(
      start_year = min(virus_year),
      n_viruses = length(unique(virus)),
      .groups = "drop"
    ) %>%
    group_by(start_year) %>%
    # Arrange by year, make sure the starting positions are unique
    mutate(
      x_position =
        start_year + (map_int(clade, ~ which(.x == unique(clade))) - 1) /
          length(unique(clade)),
    ) %>%
    ungroup() %>%
    # Stretch the length of clades with more viruses than time to next clade
    arrange(x_position) %>%
    mutate(
      next_x = lead(x_position, default = Inf),
      diff = pmax(0, n_viruses - (next_x - x_position)),
      x_position_mod = x_position + lag(cumsum(diff), default = 0),
      next_x_mod = lead(x_position_mod, default = Inf),
      clade_len = next_x_mod - x_position_mod,
      clade_len = ifelse(clade_len == Inf, n_viruses, clade_len)
    ) %>%
    select(clade, clade_start = x_position_mod, clade_len)
}

# Script ======================================================================

# HI data
hi <- read_data("hi")
hi_annette_extra <- read_data("hi-annette-extra")

# Each pid should have one virus label per x_position
hi_mod <- x_positions_by_year(hi)
hi_annette_extra_mod <- x_positions_by_year(hi_annette_extra)

# Individual plots with a simple year-based x-axis
indiv_hi_plots <- hi_mod %>%
  # filter(pid == "HIA15611") %>%
  group_by(pid, group, sex, age_lab) %>%
  group_map(plot_one_pid)

indiv_hi_plots_annette_extra <- hi_annette_extra_mod %>%
  # filter(pid == first(pid)) %>%
  group_by(pid, prior_h3_lab) %>%
  group_map(plot_one_pid)

walk(
  indiv_hi_plots, ~ save_pdf(.x, attr(.x, "name"), "indiv-hi", 45, 15)
)
walk(
  indiv_hi_plots_annette_extra,
  ~ save_pdf(.x, attr(.x, "name"), "indiv-hi-annette-extra", 45, 15)
)

# A diffrent x-axis

hi_mod_alt <- x_positions_clade_year(hi)
hi_mod_alt_annette_extra <- x_positions_clade_year(hi_annette_extra)

indiv_hi_plots_alt <- hi_mod_alt %>%
  # filter(pid == "HIA15611") %>%
  group_by(pid, group, sex, age_lab) %>%
  group_map(plot_one_pid)

indiv_hi_plots_alt_annette_extra <- hi_mod_alt_annette_extra %>%
  # filter(pid == first(pid)) %>%
  group_by(pid, dob, prior_h3_lab) %>%
  group_map(plot_one_pid)

walk(
  indiv_hi_plots_alt, ~ save_pdf(.x, attr(.x, "name"), "indiv-hi-alt", 45, 15)
)
walk(
  indiv_hi_plots_alt_annette_extra,
  ~ save_pdf(.x, attr(.x, "name"), "indiv-hi-alt-annette-extra", 45, 15)
)
