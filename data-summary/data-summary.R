# Data summaries

library(tidyverse)
library(kableExtra)

# Functions ===================================================================

source("data/read_data.R")

format_decimal <- function(d) {
  signif(d, 2)
}

format_percent <- function(d) {
  paste0(format_decimal(d) * 100, "%")
}

summarise_numeric <- function(nums) {
  nums <- na.omit(nums)
  glue::glue("{format_decimal(mean(nums))} ({format_decimal(sd(nums))})")
}

summarise_factor <- function(fac) {
  fac <- na.omit(fac)
  n <- length(fac)
  one_level <- function(lvl) {
    n_level <- sum(fac == lvl)
    glue::glue("{n_level} ({format_percent(n_level / n)}) {lvl}")
  }
  map(sort(unique(fac)), one_level) %>% paste(collapse = "; ")
}

summarise_binary <- function(bin_vec, return_everything = FALSE) {
  bin_vec <- na.omit(bin_vec)
  total <- length(bin_vec)
  success <- sum(bin_vec)
  failure <- total - success
  point <- success / total
  ci <- PropCIs::exactci(success, total, 0.95)$conf.int
  low <- ci[[1]]
  high <- ci[[2]]
  summ <- glue::glue(
    "{success} / {total} ",
    "{format_percent(point)} ({format_percent(low)}, {format_percent(high)})"
  )
  if (return_everything) {
    return(tibble(total, success, failure, point, low, high, summ))
  }
  summ
}

summarise_logmean <- function(titres) {
  titres <- na.omit(titres)
  logtitres <- log(titres)
  logmn <- mean(logtitres)
  logse <- sd(logtitres) / sqrt(length(titres))
  q <- qnorm(0.975)
  glue::glue(
    "{format_decimal(exp(logmn))} ",
    "({format_decimal(exp(logmn - q * logse))}, ",
    "{format_decimal(exp(logmn + q * logse))})",
  )
}

summarise_count <- function(count) {
  count <- na.omit(count)
  glue::glue(
    "{median(count)} [{quantile(count, 0.25)}, {quantile(count, 0.75)}] ",
    "{length(count)}"
  )
}

save_data <- function(data, name) {
  write_csv(data, glue::glue("data-summary/{name}.csv"))
  data
}

save_plot <- function(plot, name, ...) {
  ggdark::ggsave_dark(
    glue::glue("data-summary/{name}.pdf"), plot,
    units = "cm",
    ...
  )
}

save_table <- function(table, name) {
  write(table, file.path("data-summary", paste0(name, ".tex")))
}

# Script ======================================================================

cdc_participant_obj1 <- read_data("cdc-participant-obj1")

cdc_participant_obj1 %>%
  mutate(age_years = 2016 - yob) %>%
  group_by(group) %>%
  summarise(
    sex_prop = summarise_factor(sex),
    age_summ = summarise_numeric(age_years),
    .groups = "drop"
  ) %>%
  pivot_longer(-group, names_to = "summary_type", values_to = "summary") %>%
  pivot_wider(names_from = "group", values_from = "summary") %>%
  mutate(
    summary_type = recode(
      summary_type,
      "sex_prop" = "Sex",
      "age_summ" = "Age (years) in 2016"
    )
  ) %>%
  save_data("cdc-participant-summary") %>%
  kbl(
    format = "latex",
    caption = "Summaries of CDC objective 1 participants.
      Format: count (percentage); mean (sd)",
    booktabs = TRUE,
    col.names = c("", colnames(.)[-1]),
    label = "cdc-participant-summary"
  ) %>%
  save_table("cdc-participant-summary")
