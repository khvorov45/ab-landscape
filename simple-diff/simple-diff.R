cat("simple summary-measures analysis")

library(tidyverse)

data_dir <- "data"
simple_diff_dir <- "simple-diff"

# Functions ===================================================================

source(file.path(data_dir, "read_data.R"))

calc_auc <- function(hi) {
  hi %>%
    # Remove the viruses that circulated before birth
    filter(virus_year >= year_of_birth, !is.na(logtitre_mid)) %>%
    # Make sure the is one measurement per individual per timepoint
    # per virus year
    group_by(pid, group, year_of_birth, timepoint_num, virus_year) %>%
    summarise(logtitre_mid = mean(logtitre_mid), .groups = "drop") %>%
    group_by(pid, group, year_of_birth, timepoint_num) %>%
    summarise(
      area_under_curve = DescTools::AUC(virus_year, logtitre_mid),
      average_loghimid = area_under_curve / (max(virus_year) - min(virus_year)),
      .groups = "drop"
    )
}

calc_diffs <- function(aucs) {
  aucs %>%
    filter(timepoint_num %in% c(1, 2)) %>%
    select(-area_under_curve) %>%
    pivot_wider(names_from = timepoint_num, values_from = average_loghimid) %>%
    mutate(logtitre_mid_diff = `2` - `1`) %>%
    select(pid, group, logtitre_mid_diff)
}

plot_diffs <- function(diffs, ylab = "Normalised difference b/w T1 & T2") {
  diffs %>%
    ggplot(aes(group, logtitre_mid_diff)) +
    ggdark::dark_theme_bw(verbose = FALSE) +
    theme(
      strip.background = element_blank(),
      panel.spacing = unit(0, "null")
    ) +
    xlab("Group") +
    ylab(ylab) +
    geom_jitter(width = 0.05, height = 0) +
    geom_boxplot(fill = NA, col = "#0e0e2e", outlier.alpha = 0)
}

save_plot <- function(plot, name, ext = "png", ...) {
  write_csv(
    select(plot$data),
    file.path(simple_diff_dir, paste0(name, ".csv"))
  )
  ggdark::ggsave_dark(
    file.path(simple_diff_dir, paste0(name, ".", ext)), plot,
    units = "cm", ...
  )
}

process_hanam <- function(hanam, t2 = "d14", virus_threshold = 2014) {
  hanam %>%
    filter(virus_year >= virus_threshold) %>%
    mutate(group = prior_h3_lab) %>%
    filter(timepoint %in% c("BL", t2)) %>%
    mutate(
      timepoint_num = if_else(timepoint == "BL", 1, 2) %>%
        as.character() %>%
        as.numeric()
    )
}

egg_cell <- function(data, .f) {
  aucs_both <- .f(data) %>% mutate(egg = "both")
  aucs_egg_cell <- data %>%
    group_by(egg) %>%
    group_modify(~ .f(.x)) %>%
    mutate(egg = if_else(egg, "egg", "cell"))
  bind_rows(aucs_both, aucs_egg_cell)
}

# Script ======================================================================

hanam_virus_threshold <- 2014

hi <- read_data("hi")
rmh_hcw <- read_data("hi-rmh-hcw") %>% filter(group != "Moderate")
hanam <- read_data("hi-hanam") %>% process_hanam("d14", hanam_virus_threshold)
hanam_d280 <- read_data("hi-hanam") %>%
  process_hanam("d280", hanam_virus_threshold)

# Work out area under curve with midpoints
hi_aucs_all_v <- egg_cell(hi, calc_auc) %>% mutate(landscape = "after birth")
hi_aucs_same_time <- hi %>%
  filter(virus_year >= min(virus_year[egg])) %>%
  egg_cell(calc_auc) %>%
  mutate(landscape = paste0("after ", min(hi$virus_year[hi$egg])))
hi_aucs <- bind_rows(hi_aucs_all_v, hi_aucs_same_time)
rmh_hcw_aucs <- calc_auc(rmh_hcw)
hanam_aucs <- calc_auc(hanam)
hanam_d280_aucs <- calc_auc(hanam_d280)

# Differences between timepoints 1 and 2
hi_diffs <- hi_aucs %>%
  group_by(egg, landscape) %>%
  group_modify(~ calc_diffs(.x))
rmh_hcw_diffs <- calc_diffs(rmh_hcw_aucs)
hanam_diffs <- calc_diffs(hanam_aucs)
hanam_d280_diffs <- calc_diffs(hanam_d280_aucs)

# Differences b/w egg and cell
hi_cell_egg_diffs <- hi_diffs %>%
  group_by(pid, group, landscape) %>%
  summarise(
    logtitre_mid_diff =
      logtitre_mid_diff[egg == "egg"] - logtitre_mid_diff[egg == "cell"],
    .groups = "drop"
  )

# Plot the differences
hi_diffs_plot_egg_cell <- plot_diffs(hi_diffs) +
  facet_grid(landscape ~ egg) +
  labs(caption = element_blank())
hi_diffs_plot <- plot_diffs(
  filter(hi_diffs, egg == "both", landscape == "after birth")
)
rmh_hcw_diffs_plot <- plot_diffs(rmh_hcw_diffs)
hanam_diffs_plot <- plot_diffs(
  hanam_diffs,
  paste("Difference b/w BL & d14 for viruses after", hanam_virus_threshold)
)
hanam_d280_diffs_plot <- plot_diffs(
  hanam_d280_diffs,
  paste("Difference b/w BL & d280 for viruses after", hanam_virus_threshold)
)
hi_cell_egg_diffs_plot <- plot_diffs(
  hi_cell_egg_diffs,
  "Egg response minus cell response"
) +
  facet_wrap(~landscape, nrow = 1) +
  labs(caption = element_blank()) +
  geom_hline(yintercept = 0)

save_plot(hi_diffs_plot, "simple-diff", width = 10, height = 10)
save_plot(
  hi_diffs_plot_egg_cell, "simple-diff-egg-cell",
  width = 25, height = 15
)
save_plot(rmh_hcw_diffs_plot, "simple-diff-rmh-hcw", width = 10, height = 10)
save_plot(hanam_diffs_plot, "simple-diff-hanam", width = 10, height = 10)
save_plot(
  hanam_d280_diffs_plot, "simple-diff-hanam-d280",
  width = 10, height = 10
)
save_plot(
  hi_cell_egg_diffs_plot, "simple-diff-egg-minus-cell",
  width = 15, height = 10
)
