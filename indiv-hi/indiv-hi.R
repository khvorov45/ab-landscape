# Individual HI plots

library(tidyverse)

future::plan(future::multisession)

# Functions ===================================================================

source("data/read_data.R")

#' `pid_info` is a string that will go into the caption
plot_one <- function(data, pid_info = "",
                     group_var = timepoint, group_var_lab = "Timepoint") {
  group_var_q <- enquo(group_var)

  xgrid <- data %>%
    select(virus_short) %>%
    distinct()
  ygrid <- tibble(titre = 5 * 2^(0:10)) %>% filter(titre <= max(data$titre))

  clades <- data %>%
    select(virus_short, clade) %>%
    distinct()

  data %>%
    ggplot(
      aes(
        virus_short,
        titre,
        color = !!group_var_q,
        shape = !!group_var_q,
        linetype = !!group_var_q,
        fill = !!group_var_q,
        group = !!group_var_q
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
    scale_linetype_discrete(group_var_lab) +
    scale_fill_discrete(group_var_lab) +
    scale_color_discrete(group_var_lab) +
    scale_shape_manual(group_var_lab, values = c(22, 23, 24)) +
    guides(fill = guide_legend(override.aes = list(alpha = 0.5))) +
    ifelse(pid_info == "", list(), list(labs(caption = pid_info))) +
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
}

plot_longitidinal <- function(data, pid_info) {
  add_theme <- theme(
    panel.spacing = unit(0, "null"),
    strip.background = element_blank(),
    legend.position = "right"
  )
  rm_x_axis <- theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.length.x = unit(0, "null")
  )
  pl1 <- data %>%
    plot_one() +
    facet_wrap(
      ~study_year,
      ncol = 1,
      strip.position = "right",
      labeller = as_labeller(function(x) paste0("Year ", x))
    ) +
    add_theme

  pl2 <- data %>%
    plot_one("", as.factor(study_year) %>% paste0("Year ", .), "Study year") +
    facet_wrap(
      ~timepoint,
      ncol = 1,
      strip.position = "right",
    ) +
    add_theme

  # This can be empty because some participants in longitudinal objectives
  # don't have enough data
  bwyear_subset <- data %>%
    filter(
      timepoint %in% c("Pre-vax", "Post-season"),
      !(study_year == 3 & timepoint == "Post-season"),
      !(study_year == 1 & timepoint == "Pre-vax")
    ) %>%
    mutate(
      bwyear_period = if_else(
        study_year == 1 | (study_year == 2 & timepoint == "Pre-vax"), 1, 2
      ),
      vaccine_strain = if_else(study_year == 3, FALSE, vaccine_strain)
    ) %>%
    group_by(bwyear_period) %>%
    filter(length(unique(timepoint)) > 1) %>%
    ungroup()

  pl3 <- NULL
  if (nrow(bwyear_subset) == 0) {
    pl2 <- pl2 + labs(caption = pid_info)
  } else {
    pl3 <- bwyear_subset %>%
      plot_one(
        pid_info,
        group_var = timepoint %>%
          fct_rev() %>%
          recode(
            "Pre-vax" = "Pre-vax\n(next year)",
            "Post-season" = "Post-season\n(this year)"
          )
      ) +
      facet_wrap(
        ~bwyear_period,
        ncol = 1,
        labeller = as_labeller(
          function(x) recode(x, "1" = "Years 1/2", "2" = "Years 2/3")
        ),
        strip.position = "right"
      ) +
      add_theme
  }

  if (is.null(pl3)) {
    pl_all <- list(pl1 + rm_x_axis, pl2)
  } else {
    pl_all <- list(pl1 + rm_x_axis, pl2 + rm_x_axis, pl3)
  }

  cowplot::plot_grid(
    plotlist = pl_all,
    ncol = 1, align = "v"
  )
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

# CDC -------------------------------------------------------------------------

cdc_join_vaccine <- function(data) {
  data %>%
    left_join(
      read_data("cdc-vaccine") %>% mutate(vaccine_strain = TRUE),
      c("virus_full", "study_year")
    ) %>%
    mutate(
      vaccine_strain = replace_na(vaccine_strain, FALSE)
    )
}

cdc_read_collapsed_vax_hist <- function(obj) {
  read_data(glue::glue("cdc-obj{obj}-vax-hist")) %>%
    filter(status == 1) %>%
    group_by(pid) %>%
    summarise(vax_years = list(year), .groups = "drop")
}

cdc_reorder_viruses <- function(data) {
  data %>%
    mutate(
      virus_full = fct_reorder(virus_full, virus_year),
      virus_short = fct_reorder(virus_short, virus_year),
    )
}

cdc_virus <- read_data("cdc-virus")

cdc_read_hi <- function(obj) {
  read_data(glue::glue("cdc-obj{obj}-hi")) %>%
    inner_join(read_data(glue::glue("cdc-obj{obj}-participant")), "pid") %>%
    inner_join(cdc_virus, "virus_full") %>%
    cdc_join_vaccine() %>%
    inner_join(cdc_read_collapsed_vax_hist(obj), "pid") %>%
    # Filter for if I want to change things - re-doing everything takes a while
    # filter(pid == first(pid)) %>%
    cdc_reorder_viruses()
}

# Objective 1 -----------------------------------------------------------------

cdc_obj1_hi <- cdc_read_hi("1")

plots <- cdc_obj1_hi %>%
  group_by(
    pid, prior_vacs, age_first_bleed, gender, recruitment_year,
    pre_post_vax_days, post_vax_post_season_days, vax_years
  ) %>%
  group_map(function(data, key) {
    pid_info <- paste0(
      key$pid,
      " Prior vacs ", key$prior_vacs,
      " Age ", round(key$age_first_bleed),
      " ", key$gender,
      " Recrutment year ", key$recruitment_year,
      "\n",
      "Pre-Post vax days ", key$pre_post_vax_days,
      "; Post vax Post season days ", key$post_vax_post_season_days,
      "\n",
      "All vax years: ", paste(key$vax_years[[1]], collapse = " ")
    )
    plot <- plot_one(data %>% filter(timepoint != "Post-vax"), pid_info)
    attr(plot, "pid") <- key$pid
    plot
  })

if (!dir.exists("indiv-hi/cdc-obj1")) dir.create("indiv-hi/cdc-obj1")
furrr::future_walk(
  plots,
  ~ save_plot(
    .x, paste0("cdc-obj1/", attr(.x, "pid")),
    ext = "png", width = 20, height = 15
  )
)

# Objective 2 -----------------------------------------------------------------

cdc_obj2_hi <- cdc_read_hi("2")

plots_obj2 <- cdc_obj2_hi %>%
  group_by(
    pid, prior_vacs, age_first_bleed, gender, vax_years, recruitment_year
  ) %>%
  group_map(function(data, key) {
    pid_info <- paste0(
      key$pid,
      " Prior vacs ", key$prior_vacs,
      " Age ", round(key$age_first_bleed),
      " ", key$gender,
      " Recrutment year ", key$recruitment_year,
      "\n",
      "All vax years: ", paste(key$vax_years[[1]], collapse = " ")
    )

    pl_all <- plot_longitidinal(data, pid_info)

    attr(pl_all, "pid") <- key$pid
    pl_all
  })

if (!dir.exists("indiv-hi/cdc-obj2")) dir.create("indiv-hi/cdc-obj2")
furrr::future_walk(
  plots_obj2,
  ~ save_plot(
    .x, paste0("cdc-obj2/", attr(.x, "pid")),
    ext = "png", width = 20, height = 30
  )
)

# Objective 3 -----------------------------------------------------------------

cdc_obj3_hi <- cdc_read_hi("3")

# This is longitudinal sometimes
plots_obj3 <- cdc_obj3_hi %>%
  group_by(
    pid, prior_vacs, age_first_bleed, gender, vax_years, recruitment_year
  ) %>%
  group_map(function(data, key) {
    pid_info <- paste0(
      key$pid,
      " Prior vacs ", key$prior_vacs,
      " Age ", round(key$age_first_bleed),
      " ", key$gender,
      " Recruitment year ", key$recruitment_year,
      "\n",
      "All vax years: ", paste(key$vax_years[[1]], collapse = " ")
    )

    if (length(unique(data$study_year)) > 1) {
      pl <- plot_longitidinal(data, pid_info)
    } else {
      pl <- plot_one(data, pid_info)
    }

    attr(pl, "pid") <- key$pid
    pl
  })

if (!dir.exists("indiv-hi/cdc-obj3")) dir.create("indiv-hi/cdc-obj3")
furrr::future_walk(
  plots_obj3,
  ~ save_plot(
    .x, paste0("cdc-obj3/", attr(.x, "pid")),
    ext = "png", width = 20, height = 30
  )
)

# Objective 4 -----------------------------------------------------------------

cdc_obj4_hi <- cdc_read_hi("4")

# One year's data for each partcipant here
plots_obj4 <- cdc_obj4_hi %>%
  group_by(
    pid, prior_vacs, age_first_bleed, gender, recruitment_year, vax_years
  ) %>%
  group_map(function(data, key) {
    pid_info <- paste0(
      key$pid,
      " Prior vacs ", key$prior_vacs,
      " Age ", round(key$age_first_bleed),
      " ", key$gender,
      " Recrutment year ", key$recruitment_year,
      "\n",
      "All vax years: ", paste(key$vax_years[[1]], collapse = " ")
    )
    plot <- plot_one(data %>% filter(timepoint != "Post-vax"), pid_info)
    attr(plot, "pid") <- key$pid
    plot
  })

if (!dir.exists("indiv-hi/cdc-obj4")) dir.create("indiv-hi/cdc-obj4")
furrr::future_walk(
  plots_obj4,
  ~ save_plot(
    .x, paste0("cdc-obj4/", attr(.x, "pid")),
    ext = "png", width = 20, height = 15
  )
)
