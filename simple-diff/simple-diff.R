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
    group_by(pid, group, year_of_birth, timepoint, virus_year) %>%
    summarise(logtitre_mid = mean(logtitre_mid), .groups = "drop") %>%
    group_by(pid, group, year_of_birth, timepoint) %>%
    summarise(
      area_under_curve = DescTools::AUC(virus_year, logtitre_mid),
      average_loghimid = area_under_curve / (max(virus_year) - min(virus_year)),
      .groups = "drop"
    )
}

calc_diffs <- function(aucs) {
  aucs %>%
    filter(timepoint %in% c(1, 2)) %>%
    select(-area_under_curve) %>%
    pivot_wider(names_from = timepoint, values_from = average_loghimid) %>%
    mutate(logtitre_mid_diff = `2` - `1`) %>%
    select(pid, group, logtitre_mid_diff)
}

plot_diffs <- function(diffs, ylab = "Normalised difference b/w T1 & T2") {
  unique_groups <- unique(diffs$group)
  test_res <- with(diffs, {
    wilcox.test(
      logtitre_mid_diff[group == unique_groups[[1]]],
      logtitre_mid_diff[group == unique_groups[[2]]],
      exact = TRUE, conf.int = TRUE
    )
  })
  diffs %>%
    group_by(group) %>%
    mutate(
      pid_outlier = if_else(
        logtitre_mid_diff %in% DescTools::Outlier(
          logtitre_mid_diff,
          na.rm = TRUE
        ),
        pid, ""
      )
    ) %>%
    ungroup() %>%
    ggplot(aes(group, logtitre_mid_diff)) +
    ggdark::dark_theme_bw(verbose = FALSE) +
    xlab("Group") +
    ylab(ylab) +
    geom_jitter(width = 0.05, height = 0) +
    geom_text(aes(label = pid_outlier), hjust = 1, nudge_x = -0.05) +
    labs(
      caption = paste0("Mann-Whitney test p=", signif(test_res$p.value, 3))
    )
}

save_plot <- function(plot, name, ext = "pdf", ...) {
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
      timepoint = if_else(timepoint == "BL", 1, 2) %>%
        as.character() %>%
        as.numeric()
    )
}

# Script ======================================================================

hanam_virus_threshold <- 2014

hi <- read_data("hi")
rmh_hcw <- read_data("hi-rmh-hcw") %>% filter(group != "Moderate")
hanam <- read_data("hi-hanam") %>% process_hanam("d14", hanam_virus_threshold)
hanam_d280 <- read_data("hi-hanam") %>%
  process_hanam("d280", hanam_virus_threshold)

# Work out area under curve with midpoints
hi_aucs <- calc_auc(hi)
rmh_hcw_aucs <- calc_auc(rmh_hcw)
hanam_aucs <- calc_auc(hanam)
hanam_d280_aucs <- calc_auc(hanam_d280)

# Differences between timepoints 1 and 2
hi_diffs <- calc_diffs(hi_aucs)
rmh_hcw_diffs <- calc_diffs(rmh_hcw_aucs)
hanam_diffs <- calc_diffs(hanam_aucs)
hanam_d280_diffs <- calc_diffs(hanam_d280_aucs)

# Plot the differences
hi_diffs_plot <- plot_diffs(hi_diffs)
rmh_hcw_diffs_plot <- plot_diffs(rmh_hcw_diffs)
hanam_diffs_plot <- plot_diffs(
  hanam_diffs,
  paste("Difference b/w BL & d14 for viruses after", hanam_virus_threshold)
)
hanam_d280_diffs_plot <- plot_diffs(
  hanam_d280_diffs,
  paste("Difference b/w BL & d280 for viruses after", hanam_virus_threshold)
)

save_plot(hi_diffs_plot, "simple-diff", width = 10, height = 10)
save_plot(hi_diffs_plot, "simple-diff", ext = "png", width = 10, height = 10)
save_plot(rmh_hcw_diffs_plot, "simple-diff-rmh-hcw", width = 10, height = 10)
save_plot(hanam_diffs_plot, "simple-diff-hanam", width = 10, height = 10)
save_plot(
  hanam_d280_diffs_plot, "simple-diff-hanam-d280",
  width = 10, height = 10
)
