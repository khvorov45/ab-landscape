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

summarise_logmean <- function(titres, out = "string") {
  titres <- na.omit(titres)
  logtitres <- log(titres)
  logmn <- mean(logtitres)
  logse <- sd(logtitres) / sqrt(length(titres))
  q <- qnorm(0.975)
  loglow <- logmn - q * logse
  loghigh <- logmn + q * logse
  mn <- exp(logmn)
  low <- exp(loglow)
  high <- exp(loghigh)
  if (out == "string") {
    glue::glue(
      "{format_decimal(mn)} ",
      "({format_decimal(low)}, ",
      "{format_decimal(high)})",
    )
  }
  else {
    tibble(logmn, logse, loglow, loghigh, mn, low, high)
  }
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

# CDC Objective 1

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
  save_data("cdc-participant-summary-obj1") %>%
  kbl(
    format = "latex",
    caption = "Summaries of CDC objective 1 participants.
      Format: count (percentage); mean (sd)",
    booktabs = TRUE,
    col.names = c("", colnames(.)[-1]),
    label = "cdc-participant-summary-obj1"
  ) %>%
  save_table("cdc-participant-summary-obj1")

label_cdc_timepoints <- as_labeller(
  c("prevax" = "Pre-Vax", "postvax" = "Post-Vax", "postseas" = "Post-Season")
)

cdc_hi_obj1 <- read_data("cdc-hi-obj1") %>%
  inner_join(cdc_participant_obj1, by = "pid")

# Titre plot - the usual stuff, variable at prevax, rises at postvax, stays
# the same/drops slightly to postseas.
cdc_hi_obj1 %>%
  ggplot(aes(timepoint, titre, group = pid, color = pid)) +
  ggdark::dark_theme_bw(verbose = FALSE) +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    panel.spacing = unit(0, "null"),
    axis.text.x = element_text(angle = -35, hjust = 0),
    plot.margin = margin(10, 25, 10, 10)
  ) +
  facet_wrap(~virus, ncol = 7) +
  scale_y_log10(
    "Titre",
    breaks = 5 * 2^(0:15),
  ) +
  scale_x_discrete("Timepoint", expand = expansion(0.1)) +
  geom_point(alpha = 0.5, shape = 18) +
  geom_line(alpha = 0.5)

# Timepoint GMT's
cdc_obj1_timepoint_gmts <- cdc_hi_obj1 %>%
  group_by(timepoint, virus, group) %>%
  summarise(summarise_logmean(titre, out = "tibble"), .groups = "drop") %>%
  ggplot(aes(virus, mn, color = group, shape = group)) +
  ggdark::dark_theme_bw(verbose = FALSE) +
  theme(
    legend.position = "bottom",
    legend.box.spacing = unit(0, "null"),
    strip.placement = "right",
    panel.spacing = unit(0, "null"),
    strip.background = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = -45, hjust = 0),
    plot.margin = margin(10, 40, 10, 10)
  ) +
  facet_wrap(
    ~timepoint,
    ncol = 1, strip.position = "right", labeller = label_cdc_timepoints
  ) +
  scale_y_log10("GMT (95% CI)", breaks = 5 * 2^(0:15)) +
  scale_x_discrete("Virus") +
  scale_color_discrete("Group") +
  scale_shape_discrete("Group") +
  geom_pointrange(
    aes(ymin = low, ymax = high),
    position = position_dodge(width = 0.5)
  )

save_plot(
  cdc_obj1_timepoint_gmts, "cdc-obj1-timepoint-gmts",
  width = 20, height = 20
)

# Differences between timepoints
timepoint_diffs <- cdc_hi_obj1 %>%
  filter(timepoint %in% c("prevax", "postvax")) %>%
  pivot_wider(names_from = "timepoint", values_from = "titre") %>%
  mutate(ratio = postvax / prevax) %>%
  group_by(virus, group) %>%
  summarise(summarise_logmean(ratio, out = "tibble"), .groups = "drop") %>%
  ggplot(aes(virus, mn, color = group, shape = group)) +
  ggdark::dark_theme_bw(verbose = FALSE) +
  theme(
    legend.position = "bottom",
    legend.box.spacing = unit(0, "null"),
    strip.placement = "right",
    panel.spacing = unit(0, "null"),
    strip.background = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = -45, hjust = 0),
    plot.margin = margin(10, 60, 10, 10)
  ) +
  geom_pointrange(
    aes(ymin = low, ymax = high),
    position = position_dodge(width = 0.5)
  ) +
  scale_y_log10("Post/Pre vax ratio (95% CI)", breaks = 1:10) +
  scale_x_discrete("Virus") +
  scale_color_discrete("Group") +
  scale_shape_discrete("Group")

save_plot(
  timepoint_diffs, "cdc-obj1-timepoint-diffs",
  width = 20, height = 15
)

# CDC Objective 2

cdc_participant_obj2 <- read_data("cdc-participant-obj2")
cdc_vacc_hist_obj2 <- read_data("cdc-vacc-hist-obj2")

cdc_vacc_hist_obj2_summ <- cdc_vacc_hist_obj2 %>%
  group_by(pid) %>%
  summarise(vacc_number = sum(vacc_status), .groups = "drop")

cdc_participant_obj2_extra <- cdc_participant_obj2 %>%
  inner_join(cdc_vacc_hist_obj2_summ, by = "pid") %>%
  mutate(age_years = 2016 - yob)

cdc_participant_obj2_extra %>%
  summarise(
    sex_prop = summarise_factor(sex),
    age_summ = summarise_numeric(age_years),
    vax_summ = summarise_factor(vacc_number),
    .groups = "drop"
  ) %>%
  pivot_longer(everything(), names_to = "summary_type", values_to = "summary") %>%
  mutate(
    summary_type = recode(
      summary_type,
      "sex_prop" = "Sex",
      "age_summ" = "Age (years) in 2016",
      "vax_summ" = "Vaccinations in 2011-2015"
    )
  ) %>%
  save_data("cdc-participant-summary-obj2") %>%
  kbl(
    format = "latex",
    caption = "Summaries of CDC objective 2 participants.
      Format: count (percentage); mean (sd)",
    booktabs = TRUE,
    col.names = c("", "Summary"),
    label = "cdc-participant-summary-obj2"
  ) %>%
  save_table("cdc-participant-summary-obj2")


# Titres
cdc_hi_obj2 <- read_data("cdc-hi-obj2")
cdc_viruses_obj2 <- read_data("cdc-virus-obj2")

cdc_hi_obj2_extra <- cdc_hi_obj2 %>%
  inner_join(cdc_viruses_obj2, "virus_n") %>%
  inner_join(cdc_participant_obj2, "pid") %>%
  mutate(
    study_year_lbl = as.factor(study_year),
    egg_lbl = if_else(egg, "Egg", "Cell")
  )

# Titre summaries
cdc_obj2_gmts <- cdc_hi_obj2_extra %>%
  group_by(timepoint, virus_full, study_year_lbl, site) %>%
  summarise(summarise_logmean(titre, out = "tibble"), .groups = "drop")

plot1_obj2_gmts <- function(data) {
  data %>%
    ggplot(aes(virus_full, mn, color = timepoint, shape = timepoint)) +
    ggdark::dark_theme_bw(verbose = FALSE) +
    theme(
      legend.position = "bottom",
      legend.box.spacing = unit(0, "null"),
      strip.placement = "right",
      panel.spacing = unit(0, "null"),
      strip.background = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = -45, hjust = 0),
      plot.margin = margin(10, 40, 10, 10)
    ) +
    facet_wrap(
      ~study_year_lbl,
      ncol = 1, strip.position = "right",
      labeller = as_labeller(function(x) paste("Year", x))
    ) +
    scale_y_log10("GMT (95% CI)", breaks = 5 * 2^(0:15)) +
    scale_x_discrete("Virus") +
    scale_color_discrete("Timepoint", labels = label_cdc_timepoints) +
    scale_shape_discrete("Timepoint", labels = label_cdc_timepoints) +
    geom_pointrange(
      aes(ymin = low, ymax = high),
      position = position_dodge(width = 0.5)
    )
}

cdc_obj2_gmts_plot1_israel <- cdc_obj2_gmts %>%
  filter(site == "Israel") %>%
  plot1_obj2_gmts()

cdc_obj2_gmts_plot1_peru <- cdc_obj2_gmts %>%
  filter(site == "Peru") %>%
  plot1_obj2_gmts()

save_plot(
  cdc_obj2_gmts_plot1_israel, "cdc-obj2-gmts-1-israel",
  width = 20, height = 20
)

save_plot(
  cdc_obj2_gmts_plot1_peru, "cdc-obj2-gmts-1-peru",
  width = 20, height = 20
)

cdc_obj2_gmts_plot2 <- cdc_obj2_gmts %>%
  ggplot(aes(virus_full, mn, color = study_year_lbl, shape = study_year_lbl)) +
  ggdark::dark_theme_bw(verbose = FALSE) +
  theme(
    legend.position = "bottom",
    legend.box.spacing = unit(0, "null"),
    strip.placement = "right",
    panel.spacing = unit(0, "null"),
    strip.background = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = -45, hjust = 0),
    plot.margin = margin(10, 40, 10, 10)
  ) +
  facet_wrap(
    ~timepoint,
    ncol = 1, strip.position = "right",
    labeller = label_cdc_timepoints
  ) +
  scale_y_log10("GMT (95% CI)", breaks = 5 * 2^(0:15)) +
  scale_x_discrete("Virus") +
  scale_color_discrete("Study year") +
  scale_shape_discrete("Study year") +
  geom_pointrange(
    aes(ymin = low, ymax = high),
    position = position_dodge(width = 0.5)
  )

save_plot(
  cdc_obj2_gmts_plot2, "cdc-obj2-gmts-2",
  width = 20, height = 20
)

# Titre changes
cdc_hi_obj2_wide <- cdc_hi_obj2_extra %>%
  pivot_wider(names_from = "timepoint", values_from = "titre") %>%
  mutate(vax_resp = postvax / prevax)

cdc_obj2_vax_resp_virus_plot <- cdc_hi_obj2_wide %>%
  group_by(virus_full, study_year) %>%
  summarise(summarise_logmean(vax_resp, out = "tibble"), .groups = "drop") %>%
  mutate(study_year_lbl = as.factor(study_year)) %>%
  ggplot(aes(virus_full, mn, color = study_year_lbl, shape = study_year_lbl)) +
  ggdark::dark_theme_bw(verbose = FALSE) +
  theme(
    legend.position = "bottom",
    legend.box.spacing = unit(0, "null"),
    strip.placement = "right",
    panel.spacing = unit(0, "null"),
    strip.background = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = -45, hjust = 0),
    plot.margin = margin(10, 60, 10, 10)
  ) +
  geom_pointrange(
    aes(ymin = low, ymax = high),
    position = position_dodge(width = 0.5)
  ) +
  scale_y_log10("Post/Pre vax ratio (95% CI)", breaks = 1:10) +
  scale_x_discrete("Virus") +
  scale_color_discrete("Study year") +
  scale_shape_discrete("Study year")

save_plot(
  cdc_obj2_vax_resp_virus_plot, "cdc-obj2-vax-resp-virus",
  width = 25, height = 15
)

cdc_obj2_average_response <- cdc_hi_obj2_wide %>%
  filter(!is.na(vax_resp)) %>%
  group_by(pid, study_year, egg_lbl) %>%
  summarise(average_response = exp(mean(log(vax_resp))), .groups = "drop") %>%
  ggplot(aes(study_year, average_response, color = pid)) +
  ggdark::dark_theme_bw(verbose = FALSE) +
  theme(
    legend.position = "none",
    legend.box.spacing = unit(0, "null"),
    strip.placement = "right",
    panel.spacing = unit(0, "null"),
    strip.background = element_blank(),
    panel.grid.minor = element_blank(),
  ) +
  facet_wrap(~egg_lbl) +
  geom_line(alpha = 0.7) +
  geom_point(shape = 18, alpha = 0.7) +
  scale_y_log10(
    "Average post/pre vax ratio",
    breaks = c(1:5, 7, 10, 15, 20, 25)
  ) +
  scale_x_continuous("Study year", breaks = 1:3)

save_plot(
  cdc_obj2_average_response, "cdc-obj2-average-response",
  width = 15, height = 9
)
