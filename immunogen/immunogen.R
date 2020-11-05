cat("immunogenicity")

library(tidyverse)

data_dir <- "data"
immunogen_dir <- "immunogen"

# Functions ===================================================================

source(file.path(data_dir, "read_data.R"))

calc_measures <- function(t2, t1) {
  calc_mean <- function(titres, measure) {
    f <- function(n) signif(exp(n), 2)
    n <- length(titres)
    se <- sd(titres) / sqrt(n)
    tibble(
      measure = measure,
      n = n,
      point = mean(titres),
      low = point - qnorm(0.975) * se,
      high = point + qnorm(0.975) * se,
      summary = glue::glue("{f(point)} ({f(low)}, {f(high)}) [{n}]"),
    )
  }

  diff <- t2 - t1
  baseline_mean <- mean(t1)

  baseline <- calc_mean(t1, "baseline")
  gmr_uncorrected <- calc_mean(diff, "gmr_uncorrected")
  gmt_uncorrected <- calc_mean(t2, "gmt_uncorrected")

  calc_prop <- function(bin_success, measure) {
    f <- function(n) signif(n, 2)
    n <- length(bin_success)
    success <- sum(bin_success)
    failure <- n - success
    tibble(
      measure = measure,
      n = n,
      point = success / n,
      low = qbeta(0.025, success + 1, failure + 1),
      high = qbeta(0.975, success + 1, failure + 1),
      summary = glue::glue(
        "{f(point)} ({f(low)}, {f(high)}) [{success} / {n}]"
      ),
    )
  }
  seroprot_uncorrected <- calc_prop(t2 >= log(40), "seroprot_uncorrected")
  seroprot_only_below40 <- calc_prop(
    t2[t1 < log(40)] >= log(40), "seroprot_only_below40"
  )
  seroconv <- calc_prop(
    if_else(t1 < log(10), t2 >= log(40), diff >= log(4)), "seroconv"
  )

  fit <- lm(t2 ~ t1)
  b1 <- fit$coef[["t1"]]
  correction_baseline <- log(10)
  t2_corrected <- t2 - b1 * (t1 - correction_baseline)

  gmt_corrected <- calc_mean(t2_corrected, "gmt_corrected")
  gmr_corrected <- calc_mean(
    t2_corrected - correction_baseline, "gmr_corrected"
  )

  bind_rows(
    baseline,
    gmr_uncorrected,
    gmt_uncorrected,
    gmr_corrected,
    gmt_corrected,
    seroprot_uncorrected,
    seroprot_only_below40,
    seroconv,
  )
}

save_data <- function(data, name) {
  write_csv(data, glue::glue("immunogen/{name}.csv"))
  data
}

save_plot <- function(plot, name, ...) {
  ggdark::ggsave_dark(
    glue::glue("immunogen/{name}.png"),
    plot,
    units = "cm",
    ...
  )
}

arrange_plots <- function(...) {
  ggdark::lighten_geoms()
  arr <- ggpubr::ggarrange(...)
  ggdark::darken_geoms()
  arr
}

# Script ======================================================================

rmh <- read_data("hi-rmh-hcw")

rmh_egg_cell_pairs <- list(
  c("Ncast/30/16", "Sing/16-0019/16e"),
  c("Michigan/15/14", "Hkong/4801/14e"),
  c("Perth/16/09", "Perth/16/09e"),
  c("Vic/361/11", "Vic/361/11e"),
  c("Texas/50/12", "Texas/50/12e"),
  c("Switz/9715293/13", "Switz/9715293/13e")
)

rmh_immun <- rmh %>%
  filter(virus %in% flatten_chr(rmh_egg_cell_pairs)) %>%
  group_by(virus, egg) %>%
  summarise(
    calc_measures(
      logtitre[timepoint == "Post-vax"],
      logtitre[timepoint == "Pre-vax"]
    ),
    .groups = "drop"
  ) %>%
  mutate(
    virus = fct_relevel(
      virus,
      "Michigan/15/14", "Hkong/4801/14e",
      "Ncast/30/16", "Sing/16-0019/16e"
    )
  ) %>%
  arrange(virus)

rmh_immun %>%
  select(virus, measure, summary) %>%
  pivot_wider(names_from = "measure", values_from = "summary") %>%
  save_data("rmh")

# Plot egg-cell correlation ---------------------------------------------------

plot_egg_cell <- function(pair, data) {
  cell_q <- rlang::sym(pair[[1]])
  egg_q <- rlang::sym(pair[[2]])
  data %>%
    filter(virus %in% pair) %>%
    select(pid, timepoint, virus, titre) %>%
    # So that 20 -> tile from 20 to 40
    mutate(titre = exp(log(titre) + log(2) / 2)) %>%
    pivot_wider(names_from = "virus", values_from = "titre") %>%
    count(timepoint, !!cell_q, !!egg_q) %>%
    ggplot(aes(!!egg_q, !!cell_q)) +
    theme_bw() +
    theme(
      strip.background = element_blank(),
      panel.spacing = unit(0, "null"),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 30, hjust = 1),
    ) +
    scale_x_log10(breaks = 5 * 2^(0:10)) +
    scale_y_log10(breaks = 5 * 2^(0:10)) +
    scale_fill_gradient("Count", low = "lightgrey", high = "black") +
    facet_wrap(~timepoint, nrow = 1) +
    geom_tile(aes(fill = n)) +
    geom_abline(intercept = 0, slope = 1)
}

rmh_egg_cell_plots <- map(rmh_egg_cell_pairs, plot_egg_cell, rmh)

rmh_egg_cell_plots_arr <- arrange_plots(
  plotlist = rmh_egg_cell_plots, ncol = 1, common.legend = TRUE
)

save_plot(rmh_egg_cell_plots_arr, "rmh-egg-cell-corr", width = 20, height = 50)

# Look at Singapore titres for the infected

rmh_infected_rel_vir <- rmh %>%
  filter(
    virus %in% c("Ncast/30/16", "Sing/16-0019/16e")
  )

rmh_inf_plot <- rmh_infected_rel_vir %>%
  ggplot(aes(timepoint, titre, group = pid)) +
  ggdark::dark_theme_bw(verbose = FALSE) +
  theme(strip.background = element_blank()) +
  scale_y_log10("Titre", breaks = 5 * 2^(0:10)) +
  scale_x_discrete("Timepoint", expand = expansion(mult = 0.2)) +
  facet_wrap(~virus) +
  geom_line(alpha = 0.4) +
  geom_point(alpha = 0.4, shape = 16) +
  geom_line(
    data = filter(rmh_infected_rel_vir, infected),
    color = "red"
  ) +
  geom_point(
    data = filter(rmh_infected_rel_vir, infected),
    color = "red"
  )

save_plot(rmh_inf_plot, "rmh-inf", width = 12, height = 7)
