cat("plot the data\n")

library(tidyverse)
library(furrr)

plan(multiprocess)

data_dir <- "data"
data_plot_dir <- "data-plot"

# Functions ===================================================================

source(file.path(data_dir, "read_data.R"))

minmax <- function(vals, length.out = 50) {
  seq(min(vals), max(vals), length.out = length.out)
}

fit_loess <- function(data, to_predict, span = 2) {
  loess(logtitre_mid ~ ag_x_coord + ag_y_coord, data, span = span) %>%
    predict(to_predict, se = TRUE) %>%
    as_tibble() %>%
    bind_cols(to_predict)
}

# Expects to be used as a function object in group_map
plot_one_pid <- function(data, key, name_gen = function(key) paste(key$pid)) {
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
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 9),
      text = element_text(size = 20),
      legend.position = "bottom",
      legend.box.spacing = unit(0, "null"),
      panel.grid.minor = element_blank(),
      strip.background = element_blank()
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
    )
  if ("vaccine_strain" %in% names(data)) {
    plot <- plot +
      geom_vline(
        aes(xintercept = x_position),
        data = filter(data, vaccine_strain),
        alpha = 0.8, col = "blue", lty = "11"
      )
  }
  plot <- plot +
    geom_line() +
    geom_point()
  attr(plot, "name") <- name_gen(key)
  plot
}

plot_contour <- function(data, key,
                         name_gen = function(key) paste(key$pid),
                         span = 2) {
  data <- data %>%
    filter(!is.na(ag_x_coord), !is.na(ag_y_coord))
  to_predict <- expand.grid(
    ag_x_coord = minmax(data$ag_x_coord),
    ag_y_coord = minmax(data$ag_y_coord),
    KEEP.OUT.ATTRS = FALSE
  ) %>%
    as_tibble()
  plot <- data %>%
    select(timepoint, ag_x_coord, ag_y_coord, logtitre_mid) %>%
    group_by(timepoint) %>%
    group_modify(~ fit_loess(.x, to_predict, span)) %>%
    ggplot(aes(ag_x_coord, ag_y_coord)) +
    ggdark::dark_theme_bw(verbose = FALSE) +
    theme(strip.background = element_blank()) +
    scale_fill_continuous(
      "Titre",
      type = "viridis", breaks = log(5 * 2^(0:10)), labels = 5 * 2^(0:10)
    ) +
    facet_wrap(~timepoint, nrow = 1) +
    scale_x_continuous("X", expand = expansion()) +
    scale_y_continuous("Y", expand = expansion()) +
    labs(caption = pmap(key, ~ paste(...))) +
    geom_tile(aes(fill = fit)) +
    geom_point(data = data) +
    ggrepel::geom_label_repel(
      data = data,
      mapping = aes(label = virus), fill = "#ffffff56",
      size = 2.5
    )
  if ("vaccine_strain" %in% names(data)) {
    plot <- plot + ggrepel::geom_label_repel(
      data = data %>% filter(vaccine_strain),
      mapping = aes(label = virus), fill = "#ffffff56", color = "blue"
    )
  }
  attr(plot, "name") <- name_gen(key)
  plot
}

# Returns a list with each entry being a plot by pid
# Add grouping variables for the caption
plots_by_pid <- function(data, ..., .plot_fun = plot_one_pid) {
  data %>%
    # filter(pid == first(pid)) %>%
    group_by(pid, ...) %>%
    group_map(.plot_fun)
}

save_plot <- function(plot,
                      name,
                      width = 45, height = 15,
                      device = "pdf",
                      ...) {
  ggdark::ggsave_dark(
    paste0(name, ".", device),
    plot,
    dark = FALSE,
    width = width,
    height = height,
    units = "cm",
    ...
  )
}

save_plots <- function(plots, dir,
                       width = 45, height = 15,
                       device = "png",
                       ...) {
  plotdir <- file.path(data_plot_dir, dir)
  if (!dir.exists(plotdir)) dir.create(plotdir)
  future_map(
    plots,
    ~ save_plot(.x, file.path(plotdir, attr(.x, "name")), width, height, device)
  )
  invisible(NULL)
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

study_year_lab_facets <- function(x) {
  x +
    facet_wrap(~study_year_lab, ncol = 1, strip.position = "right") +
    theme(panel.spacing = unit(0, "null"))
}

# Script ======================================================================

# Split into multiple files, otherwise takes too long
