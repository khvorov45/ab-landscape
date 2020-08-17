cat("simple summary-measures analysis")

library(tidyverse)

data_dir <- "data"
simple_diff_dir <- "simple-diff"

# Functions ===================================================================

source(file.path(data_dir, "read_data.R"))

calc_auc <- function(hi) {
  hi %>%
    # Remove the viruses that circulated before birth
    filter(virus_year >= year_of_birth) %>%
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

plot_diffs <- function(diffs) {
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
        logtitre_mid_diff %in% DescTools::Outlier(logtitre_mid_diff), pid, ""
      )
    ) %>%
    ungroup() %>%
    ggplot(aes(group, logtitre_mid_diff)) +
    ggdark::dark_theme_bw(verbose = FALSE) +
    xlab("Group") +
    ylab("Normalised difference b/w T1 & T2") +
    geom_jitter(width = 0.05, height = 0) +
    geom_text(aes(label = pid_outlier), hjust = 1, nudge_x = -0.05) +
    labs(
      caption = paste0("Mann-Whitney test p=", signif(test_res$p.value, 3))
    )
}

save_pdf <- function(plot, name, ...) {
  ggdark::ggsave_dark(
    file.path(simple_diff_dir, paste0(name, ".pdf")), plot,
    units = "cm", ...
  )
}

# Script ======================================================================

hi <- read_data("hi")
rmh_hcw <- read_data("hi-rmh-hcw") %>% filter(group != "Moderate")

# Work out area under curve with midpoints
hi_aucs <- calc_auc(hi)
rmh_hcw_aucs <- calc_auc(rmh_hcw)

# Differences between timepoints 1 and 2
hi_diffs <- calc_diffs(hi_aucs)
rmh_hcw_diffs <- calc_diffs(rmh_hcw_aucs)

# Plot the differences
hi_diffs_plot <- plot_diffs(hi_diffs)
rmh_hcw_diffs_plot <- plot_diffs(rmh_hcw_diffs)

save_pdf(hi_diffs_plot, "simple-diff", width = 10, height = 10)
save_pdf(rmh_hcw_diffs_plot, "simple-diff-rmh-hcw", width = 10, height = 10)
