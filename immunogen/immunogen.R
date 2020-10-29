cat("immunogenicity")

library(tidyverse)

data_dir <- "data"
immunogen_dir <- "immunogen"

# Functions ===================================================================

source(file.path(data_dir, "read_data.R"))

calc_measures <- function(t2, t1) {
  f <- function(n) signif(n, 2)
  calc_mean <- function(titres, measure) {
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

  gmr_uncorrected <- calc_mean(diff, "gmr_uncorrected")
  gmt_uncorrected <- calc_mean(t2, "gmt_uncorrected")

  calc_prop <- function(bin_success, measure) {
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
  seropos_uncorrected <- calc_prop(t2 >= log(40), "seropos_uncorrected")
  seropos_only_below40 <- calc_prop(
    t2[t1 < log(40)] >= log(40), "seropos_only_below40"
  )
  seroconv <- calc_prop(
    if_else(t1 < log(10), t2 >= log(40), diff >= log(4)), "seroconv"
  )

  fit <- lm(t2 ~ t1)
  b1 <- fit$coef[["t1"]]
  t2_corrected <- t2 - b1 * t1

  # This can be interpreted as both the corrected gmt and the corrected gmr
  gmtr_corrected <- calc_mean(t2_corrected, "gmtr_corrected")

  bind_rows(
    gmr_uncorrected,
    gmt_uncorrected,
    gmtr_corrected,
    seropos_uncorrected,
    seropos_only_below40,
    seroconv,
  )
}

# Script ======================================================================

rmh <- read_data("hi-rmh-hcw")

rmh_immun <- rmh %>%
  mutate(virus_temp = str_replace(virus, "e$", "")) %>%
  # Keep only the viruses that have both egg and cell data
  group_by(virus_temp) %>%
  filter(any(egg) & any(!egg)) %>%
  group_by(virus, egg) %>%
  summarise(
    calc_measures(
      logtitre_mid[timepoint == "Post-vax"],
      logtitre_mid[timepoint == "Pre-vax"]
    ),
    .groups = "drop"
  )

rmh_immun %>%
  select(virus, measure, summary) %>%
  pivot_wider(names_from = "measure", values_from = "summary")
