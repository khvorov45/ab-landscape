# Average the landscapes

source("data-plot/data-plot.R")

data <- map(
  c("hi" = "hi", "rmh" = "hi-rmh-hcw"),
  ~ read_data(.x) %>%
    select(
      pid, group, timepoint, virus,
      virus_year, clade,
      ag_x_coord, ag_y_coord, titre,
    ),
)

data_summ <- map(data, function(data) {
  data %>%
    group_by(
      group, timepoint, virus, virus_year, clade,
      ag_x_coord, ag_y_coord
    ) %>%
    summarise(titre = mean(log(titre)) %>% exp(), .groups = "drop")
})

# X positions for the linear landscape
data_mod <- map(data_summ, x_positions_by_year)
data_mod_alt <- map(data_summ, x_positions_clade_year)

# Plot
add_facets <- function(plot) {
  plot +
    facet_wrap(~group, ncol = 1)
}
lin <- data_mod %>%
  map(plot_lin_landscape) %>%
  map(add_facets)

lin_alt <- data_mod_alt %>%
  map(plot_lin_landscape) %>%
  map(add_facets)

add_facets_2 <- function(plot) {
  plot +
    facet_wrap(~timepoint, ncol = 1)
}
lin_2 <- data_mod %>%
  map(plot_lin_landscape, group_var = group, group_var_lab = "Group") %>%
  map(add_facets_2)
lin_2_alt <- data_mod_alt %>%
  map(plot_lin_landscape, group_var = group, group_var_lab = "Group") %>%
  map(add_facets_2)

save_plot_averaged <- function(plot,
                               name,
                               suff = "",
                               heights = c("rmh" = 30, "hi" = 20),
                               ...) {
  save_plot(
    plot,
    glue::glue("data-plot/averaged-{name}{suff}"),
    width = 37,
    height = heights[[name]],
    ...
  )
}
iwalk(lin, save_plot_averaged)
iwalk(lin_alt, save_plot_averaged, suff = "-alt")
iwalk(
  lin_2, save_plot_averaged,
  suff = "-2", heights = c("rmh" = 30, "hi" = 30)
)
iwalk(
  lin_2_alt, save_plot_averaged,
  suff = "-2-alt",
  heights = c("rmh" = 30, "hi" = 30)
)

# Contours
add_cont_facets <- function(plot) {
  plot +
    facet_grid(group ~ timepoint)
}
cont <- data_summ %>%
  map(plot_contour, span = 4) %>%
  map(add_cont_facets)
iwalk(cont, save_plot_averaged, "-contour")
