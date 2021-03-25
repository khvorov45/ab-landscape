# Individual HI plots

library(tidyverse)

future::plan(future::multisession)

# Functions ===================================================================

source("data/read_data.R")

plot_one <- function(data) {
  pid_info <- data %>%
    select(pid, age_first_bleed, gender, prior_vacs, study_year) %>%
    distinct()

  xgrid <- data %>%
    select(virus_short) %>%
    distinct()
  ygrid <- tibble(titre = 5 * 2^(0:10)) %>% filter(titre <= max(data$titre))

  clades <- data %>%
    select(virus_short, clade) %>%
    distinct()

  plot <- data %>%
    ggplot(
      aes(
        virus_short,
        titre,
        color = timepoint, shape = timepoint, linetype = timepoint,
        fill = timepoint,
        group = paste0(pid, timepoint)
      )
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom",
      legend.box.spacing = unit(0, "null"),
      panel.grid = element_blank()
    ) +
    scale_y_log10("Titre", breaks = 5 * 2^(0:10), expand = expansion(c(0, 0.1))) +
    scale_x_discrete("Virus", expand = expansion(0.01)) +
    scale_fill_discrete("Timepoint") +
    scale_color_discrete("Timepoint") +
    scale_linetype_discrete("Timepoint") +
    scale_shape_manual("Timepoint", values = c(22, 23, 24)) +
    guides(fill = guide_legend(override.aes = list(alpha = 0.5))) +
    labs(
      caption = paste0(
        pid_info$pid,
        " Prior vacs ", pid_info$prior_vacs,
        " Age ", round(pid_info$age_first_bleed),
        " ", pid_info$gender,
        " Year ", pid_info$study_year
      )
    ) +
    # Ribbon should be before gridlines
    geom_ribbon(aes(ymin = 5, ymax = titre)) +
    # Grid lines
    geom_vline(
      aes(xintercept = virus_short),
      data = xgrid,
      alpha = 0.3,
      color = "gray50"
    ) +
    geom_hline(
      aes(yintercept = titre),
      data = ygrid,
      alpha = 0.3,
      color = "gray50"
    ) +
    # Vaccine strain
    geom_vline(
      aes(xintercept = virus_short),
      data = data %>% filter(vaccine_strain),
      alpha = 0.6
    ) +
    # Main geoms
    geom_line(color = "gray50", lty = "solid", size = 1.1, alpha = 1) +
    geom_line() +
    geom_point(fill = "black", color = "black", alpha = 1, size = 2) +
    geom_point() +
    geom_text(
      aes(virus_short, 5.5, label = ifelse(clade == "(missing)", "", clade)),
      vjust = 0.5, hjust = 0,
      data = clades,
      angle = 90,
      inherit.aes = FALSE,
      color = "gray25"
    )

  attr(plot, "pid") <- pid_info$pid
  plot
}

save_plot <- function(plot, name, ext = "pdf", ...) {
  ggsave(
    paste0("indiv-hi/", name, ".", ext),
    plot,
    units = "cm",
    ...
  )
}

# Script ======================================================================

cdc_obj1_hi <- read_data("cdc-obj1-hi") %>%
  inner_join(read_data("cdc-obj1-participant"), "pid") %>%
  inner_join(read_data("cdc-virus"), "virus_full") %>%
  left_join(
    read_data("cdc-vaccine") %>% mutate(vaccine_strain = TRUE),
    c("virus_full", "study_year")
  ) %>%
  # Filter for if I want to change things - re-doing everything takes a while
  # filter(pid == first(pid)) %>%
  mutate(
    virus_full = fct_reorder(virus_full, virus_year),
    virus_short = fct_reorder(virus_short, virus_year),
    vaccine_strain = replace_na(vaccine_strain, FALSE)
  )

plots <- cdc_obj1_hi %>%
  group_split(pid) %>%
  map(~ plot_one(.x), .keep = TRUE)
plots
if (!dir.exists("indiv-hi/cdc-obj1")) dir.create("indiv-hi/cdc-obj1")
furrr::future_walk(
  plots,
  ~ save_plot(
    .x, paste0("cdc-obj1/", attr(.x, "pid")),
    ext = "png", width = 20, height = 15
  )
)
