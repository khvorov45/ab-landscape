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
  seroporot_only_below40 <- calc_prop(
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
    gmt_corrected,
    gmr_corrected,
    seroprot_uncorrected,
    seroporot_only_below40,
    seroconv,
  )
}

save_data <- function(data, name) {
  write_csv(data, glue::glue("immunogen/{name}.csv"))
  data
}

# Script ======================================================================

rmh <- read_data("hi-rmh-hcw")

rmh_immun <- rmh %>%
  mutate(virus_temp = str_replace(virus, "e$", "")) %>%
  # Keep only the viruses that have both egg and cell data
  group_by(virus_temp) %>%
  filter(
    (any(egg) & any(!egg)) |
      virus %in% c(
        "Michigan/15/14", "Hkong/4801/14e",
        "Ncast/30/16", "Sing/16-0019/16e"
      )
  ) %>%
  group_by(virus, egg) %>%
  summarise(
    calc_measures(
      logtitre_mid[timepoint == "Post-vax"],
      logtitre_mid[timepoint == "Pre-vax"]
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
