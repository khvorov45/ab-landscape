cat("simple summary-measures analysis")

library(tidyverse)

data_dir <- "data"
simple_diff_dir <- "simple-diff"

# Functions ===================================================================

source(file.path(data_dir, "read_data.R"))

save_pdf <- function(plot, name, ...) {
  ggdark::ggsave_dark(
    file.path(simple_diff_dir, paste0(name, ".pdf")),
    units = "cm", ...
  )
}

# Script ======================================================================

hi <- read_data("hi")

# Work out area under curve with midpoints
aucs <- hi %>%
  # Make sure the is one measurement per individual per timepoint per year
  group_by(pid, group, timepoint, virus_year) %>%
  summarise(logtitre_mid = mean(logtitre_mid), .groups = "drop") %>%
  group_by(pid, group, timepoint) %>%
  summarise(
    area_under_curve = DescTools::AUC(virus_year, logtitre_mid),
    .groups = "drop"
  )

# Differences between timepoints 1 and 2
diffs <- aucs %>%
  filter(timepoint %in% c(1, 2)) %>%
  pivot_wider(names_from = timepoint, values_from = area_under_curve) %>%
  mutate(logtitre_mid_diff = `2` - `1`) %>%
  select(pid, group, logtitre_mid_diff)

# Compare the differences
test_res <- with(diffs, {
  wilcox.test(
    logtitre_mid_diff[group == "Case"],
    logtitre_mid_diff[group == "Control"],
    exact = TRUE, conf.int = TRUE
  )
})

# Plot the differences
diffs_plot <- diffs %>%
  mutate(
    pid_outlier = if_else(
      logtitre_mid_diff >= (25 + 15 * (group == "Case")), pid, ""
    )
  ) %>%
  ggplot(aes(group, logtitre_mid_diff)) +
  ggdark::dark_theme_bw(verbose = FALSE) +
  xlab("Group") +
  ylab("AUC difference b/w T1 & T2") +
  geom_jitter(width = 0.05, height = 0) +
  geom_text(aes(label = pid_outlier), hjust = 1, nudge_x = -0.05) +
  labs(
    caption = paste0("Mann-Whitney test p=", signif(test_res$p.value, 3))
  )

save_pdf(diffs_plot, "simple-diff", width = 10, height = 10)
