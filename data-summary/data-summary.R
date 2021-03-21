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

#' Takes non-logged titres
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

cdc_viruses <- read_data("cdc-virus")
cdc_vaccine_obj1 <- read_data("cdc-vaccine-obj1")

# Viruses in each clade
cdc_viruses %>%
  filter(!clade %in% c("1", "2", "(missing)")) %>%
  mutate(egg_lbl = if_else(egg, "Egg", "Cell")) %>%
  select(Clade = clade, Virus = virus_full, Type = egg_lbl) %>%
  kbl(
    format = "latex",
    caption = "Viruses in each clade.",
    booktabs = TRUE,
    label = "cdc-clade-viruses"
  ) %>%
  collapse_rows(columns = 1, latex_hline = "major") %>%
  save_table("cdc-clade-viruses")

cdc_hi_obj1 <- read_data("cdc-hi-obj1") %>%
  inner_join(cdc_participant_obj1, by = "pid") %>%
  inner_join(cdc_viruses, "virus_full") %>%
  mutate(
    vaccine_strain = virus_full %in% cdc_vaccine_obj1$virus_full,
    egg_lbl = if_else(egg, "Egg", "Cell")
  )

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
  facet_wrap(~virus_full, ncol = 7) +
  scale_y_log10(
    "Titre",
    breaks = 5 * 2^(0:15),
  ) +
  scale_x_discrete("Timepoint", expand = expansion(0.1)) +
  geom_point(alpha = 0.5, shape = 18) +
  geom_line(alpha = 0.5)

# Timepoint GMT's
cdc_obj1_timepoint_gmts <- cdc_hi_obj1 %>%
  group_by(timepoint, virus_full, group) %>%
  summarise(summarise_logmean(titre, out = "tibble"), .groups = "drop") %>%
  ggplot(aes(virus_full, mn, color = group, shape = group)) +
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
  geom_vline(
    xintercept = cdc_hi_obj1$virus_full[cdc_hi_obj1$vaccine_strain],
    size = 4, alpha = 0.2
  ) +
  geom_pointrange(
    aes(ymin = low, ymax = high),
    position = position_dodge(width = 0.5)
  )

save_plot(
  cdc_obj1_timepoint_gmts, "cdc-obj1-timepoint-gmts",
  width = 20, height = 20
)

# Differences between timepoints
# Highlight vaccine strain here as well
timepoint_diffs <- cdc_hi_obj1 %>%
  select(-bleed_date) %>%
  filter(timepoint %in% c("prevax", "postvax")) %>%
  pivot_wider(names_from = "timepoint", values_from = "titre") %>%
  mutate(ratio = postvax / prevax) %>%
  group_by(virus_full, group) %>%
  summarise(summarise_logmean(ratio, out = "tibble"), .groups = "drop") %>%
  ggplot(aes(virus_full, mn, color = group, shape = group)) +
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
  geom_vline(
    xintercept = cdc_hi_obj1$virus_full[cdc_hi_obj1$vaccine_strain],
    size = 4, alpha = 0.2
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
  group_by(site) %>%
  summarise(
    sex_prop = summarise_factor(sex),
    age_summ = summarise_numeric(age_years),
    vax_summ = summarise_factor(vacc_number),
    .groups = "drop"
  ) %>%
  pivot_longer(-site, names_to = "summary_type", values_to = "summary") %>%
  pivot_wider(names_from = "site", values_from = "summary") %>%
  mutate(
    summary_type = recode(
      summary_type,
      "sex_prop" = "Sex",
      "age_summ" = "Age (years) in 2016",
      "vax_summ" = "Vaccinations in 2011-15"
    )
  ) %>%
  save_data("cdc-participant-summary-obj2") %>%
  kbl(
    format = "latex",
    caption = "Summaries of CDC objective 2 participants.
      Format: count (percentage); mean (sd)",
    booktabs = TRUE,
    col.names = c("", colnames(.)[-1]),
    label = "cdc-participant-summary-obj2"
  ) %>%
  kable_styling(latex_options = "scale_down") %>%
  save_table("cdc-participant-summary-obj2")


# Titres
cdc_hi_obj2 <- read_data("cdc-hi-obj2")
cdc_viruses <- read_data("cdc-virus")
cdc_vaccine_obj2 <- read_data("cdc-vaccine-obj2")

cdc_hi_obj2_extra <- cdc_hi_obj2 %>%
  inner_join(cdc_viruses, "virus_n") %>%
  inner_join(cdc_participant_obj2, "pid") %>%
  mutate(
    study_year_lbl = as.factor(study_year) %>% paste0("Year ", .),
    egg_lbl = if_else(egg, "Egg", "Cell"),
  ) %>%
  left_join(
    mutate(
      cdc_vaccine_obj2,
      vaccine_strain = TRUE
    ),
    c("study_year", "virus_full")
  ) %>%
  mutate(vaccine_strain = replace_na(vaccine_strain, FALSE))

# Titre summaries
cdc_obj2_gmts <- cdc_hi_obj2_extra %>%
  group_by(timepoint, virus_full, vaccine_strain, study_year_lbl, site) %>%
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
    geom_vline(
      aes(xintercept = virus_full),
      data = data %>% filter(vaccine_strain),
      size = 6, alpha = 0.2
    ) +
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

plot2_obj2_gmts <- function(data) {
  data %>%
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
    geom_vline(
      aes(xintercept = virus_full),
      data = data %>% filter(vaccine_strain),
      size = 6, alpha = 0.2
    ) +
    geom_pointrange(
      aes(ymin = low, ymax = high),
      position = position_dodge(width = 0.5)
    )
}

cdc_obj2_gmts_plot2_israel <- cdc_obj2_gmts %>%
  filter(site == "Israel") %>%
  plot2_obj2_gmts()

cdc_obj2_gmts_plot2_peru <- cdc_obj2_gmts %>%
  filter(site == "Peru") %>%
  plot2_obj2_gmts()

save_plot(
  cdc_obj2_gmts_plot2_israel, "cdc-obj2-gmts-2-israel",
  width = 20, height = 20
)

save_plot(
  cdc_obj2_gmts_plot2_peru, "cdc-obj2-gmts-2-peru",
  width = 20, height = 20
)

# Titre changes
cdc_hi_obj2_wide <- cdc_hi_obj2_extra %>%
  pivot_wider(names_from = "timepoint", values_from = "titre") %>%
  mutate(vax_resp = postvax / prevax)

cdc_obj2_vax_resp_virus_plot <- cdc_hi_obj2_wide %>%
  group_by(virus_full, study_year, site, vaccine_strain) %>%
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
  facet_wrap(~site, ncol = 1, strip.position = "right") +
  geom_vline(
    aes(xintercept = virus_full),
    data = . %>% filter(vaccine_strain),
    size = 6, alpha = 0.2
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
  width = 20, height = 20
)

cdc_obj2_average_response <- cdc_hi_obj2_wide %>%
  filter(!is.na(vax_resp)) %>%
  group_by(pid, study_year, egg_lbl, site) %>%
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
  facet_grid(site ~ egg_lbl) +
  geom_line(alpha = 0.7) +
  geom_point(shape = 18, alpha = 0.7) +
  scale_y_log10(
    "Average post/pre vax ratio",
    breaks = c(1:5, 7, 10, 15, 20, 25)
  ) +
  scale_x_continuous("Study year", breaks = 1:3)

save_plot(
  cdc_obj2_average_response, "cdc-obj2-average-response",
  width = 15, height = 15
)

# Circulating clades
cdc_clade_frequencies <- read_data("cdc-clade-frequencies")

label_years <- function(year) {
  year <- as.numeric(year)
  glue::glue(
    "{floor(year)} ({ifelse(year %% 1 == 0.25, '1st', '2nd')} half)"
  )
}

cdc_clade_freq_plot <- cdc_clade_frequencies %>%
  filter(clade %in% cdc_viruses$clade) %>%
  ggplot(aes(year, freq)) +
  ggdark::dark_theme_bw(verbose = FALSE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_continuous(
    "Year",
    breaks = unique(cdc_clade_frequencies$year),
    labels = as_labeller(label_years),
  ) +
  scale_fill_discrete("Clade") +
  scale_y_continuous(
    "Frequency",
    labels = scales::percent_format(1),
    expand = expansion(c(0, 0.05))
  ) +
  geom_bar(aes(fill = clade), stat = "identity", color = "black")

save_plot(cdc_clade_freq_plot, "cdc-clade-freq", width = 15, height = 10)

# The 'season' was the first half of 2019
cdc_obj1_bleed_dates <- cdc_hi_obj1 %>%
  ggplot(aes(bleed_date, pid, color = timepoint)) +
  ggdark::dark_theme_bw(verbose = FALSE) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    legend.box.spacing = unit(0, "null")
  ) +
  scale_x_date("Bleed date", breaks = "month") +
  scale_color_discrete("Timepoint", labels = label_cdc_timepoints) +
  geom_point()

save_plot(cdc_obj1_bleed_dates, "cdc-obj1-bleed-dates", width = 15, height = 15)

cdc_obj1_clades_summ <- cdc_hi_obj1 %>%
  filter(!clade %in% c("(missing)", "1", "2")) %>%
  select(clade, egg_lbl, group) %>%
  distinct() %>%
  inner_join(filter(cdc_clade_frequencies, year == 2019.25), "clade") %>%
  filter(freq > 0) %>%
  group_by(egg_lbl, group) %>%
  summarise(clades = list(clade), freq_total = sum(freq), .groups = "drop")

# This is supposed to show titres against what was circulating at the time
cdc_obj1_hi_against_circulating <- cdc_hi_obj1 %>%
  filter(clade %in% cdc_clade_frequencies$clade) %>%
  group_by(pid, timepoint, clade, group, egg_lbl) %>%
  summarise(
    mean_log_titre = mean(log(titre)),
    .groups = "drop"
  ) %>%
  inner_join(filter(cdc_clade_frequencies, year == 2019.25), "clade") %>%
  group_by(pid, timepoint, group, egg_lbl) %>%
  summarise(
    wmean_log_titre = sum(mean_log_titre * freq) / sum(freq),
    wmean_titre = exp(wmean_log_titre),
    .groups = "drop"
  )

cdc_obj1_ind_av_circulating <- cdc_obj1_hi_against_circulating %>%
  ggplot(aes(timepoint, wmean_titre, group = pid, color = pid)) +
  ggdark::dark_theme_bw(verbose = FALSE) +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    panel.spacing = unit(0, "null"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_y_log10(
    "Average titre against cirulating strains",
    breaks = 5 * 2^(0:10)
  ) +
  scale_x_discrete(
    "Timepoint",
    labels = label_cdc_timepoints,
    expand = expansion(0.1)
  ) +
  facet_grid(egg_lbl ~ group) +
  geom_line(alpha = 0.5) +
  geom_point() +
  geom_text(
    aes(1, 640, label = paste0(signif(freq_total * 100, 2), "%")),
    data = cdc_obj1_clades_summ,
    inherit.aes = FALSE
  )

save_plot(
  cdc_obj1_ind_av_circulating, "cdc-obj1-ind-av-circulating",
  width = 15, height = 15
)

cdc_obj1_gmt_circulating <- cdc_obj1_hi_against_circulating %>%
  group_by(group, timepoint, egg_lbl) %>%
  summarise(
    summarise_logmean(wmean_titre, out = "tibble"),
    .groups = "drop"
  ) %>%
  ggplot(aes(timepoint, mn, color = group)) +
  ggdark::dark_theme_bw(verbose = FALSE) +
  theme(
    legend.position = "bottom",
    legend.box.spacing = unit(0, "null"),
    strip.background = element_blank(),
    panel.spacing = unit(0, "null")
  ) +
  facet_wrap(~egg_lbl) +
  scale_y_log10(
    "GMT against circulating strains",
    breaks = 5 * 2^(0:10),
  ) +
  scale_x_discrete("Timepoint", labels = label_cdc_timepoints) +
  scale_color_discrete("Group") +
  geom_pointrange(
    aes(ymin = low, ymax = high),
    position = position_dodge(width = 0.4)
  )

save_plot(
  cdc_obj1_gmt_circulating, "cdc-obj1-gmt-circulating",
  width = 15, height = 10
)

# See what the seasons were in objective 2
cdc_obj2_bleed_dates_plot <- cdc_hi_obj2_extra %>%
  ggplot(aes(bleed_date, pid, color = timepoint)) +
  ggdark::dark_theme_bw(verbose = FALSE) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_blank(),
    strip.placement = "right",
    legend.box.spacing = unit(0, "null")
  ) +
  scale_x_date("Bleed date", breaks = "month") +
  scale_y_discrete("PID") +
  scale_color_discrete("Timepoint", labels = label_cdc_timepoints) +
  geom_point() +
  geom_hline(yintercept = 15.5) +
  annotate(
    geom = "text",
    x = lubridate::ymd("2019-04-01"),
    y = 30, label = "Peru"
  ) +
  annotate(
    geom = "text",
    x = lubridate::ymd("2016-06-01"),
    y = 9, label = "Israel"
  )

save_plot(
  cdc_obj2_bleed_dates_plot, "cdc-obj2-bleed-dates",
  width = 20, height = 20
)

cdc_hi_obj2_extra_circulating <- cdc_hi_obj2_extra %>%
  mutate(
    year = case_when(
      study_year == 1 & site == "Peru" ~ 2016.75,
      study_year == 1 & site == "Israel" ~ 2017.25,
      study_year == 2 & site == "Peru" ~ 2017.75,
      study_year == 2 & site == "Israel" ~ 2018.25,
      study_year == 3 & site == "Peru" ~ 2018.75,
      study_year == 3 & site == "Israel" ~ 2019.25,
      TRUE ~ NA_real_
    )
  )

cdc_obj2_clades_summ <- cdc_hi_obj2_extra_circulating %>%
  filter(!clade %in% c("(missing)", "1", "2")) %>%
  select(clade, egg_lbl, site, year, study_year_lbl) %>%
  distinct() %>%
  inner_join(cdc_clade_frequencies, c("clade", "year")) %>%
  filter(freq > 0) %>%
  group_by(egg_lbl, site, year, study_year_lbl) %>%
  summarise(clades = list(clade), freq_total = sum(freq), .groups = "drop")

cdc_obj2_ind_circulating <- cdc_hi_obj2_extra_circulating %>%
  group_by(
    pid, site, study_year, study_year_lbl, timepoint, clade, egg, egg_lbl
  ) %>%
  summarise(logtitre_mean = mean(log(titre)), .groups = "drop") %>%
  inner_join(cdc_clade_frequencies, c("clade", "year")) %>%
  group_by(pid, study_year, study_year_lbl, site, timepoint, egg, egg_lbl) %>%
  summarise(
    wmean_titre = exp(sum(logtitre_mean * freq) / sum(freq)),
    .groups = "drop"
  )

cdc_obj2_ind_circulating_plot <- cdc_obj2_ind_circulating %>%
  ggplot(aes(timepoint, wmean_titre, group = pid, color = pid)) +
  ggdark::dark_theme_bw(verbose = FALSE) +
  theme(
    strip.background = element_blank(),
    panel.spacing = unit(0, "null"),
    legend.position = "none"
  ) +
  facet_grid(site + egg_lbl ~ study_year_lbl) +
  scale_x_discrete("Timepoint", labels = label_cdc_timepoints) +
  scale_y_log10("Average titre against circulating", breaks = 5 * 2^(0:10)) +
  geom_line(alpha = 0.4) +
  geom_point() +
  geom_text(
    aes(1, 500, label = signif(freq_total * 100, 2) %>% paste0("%")),
    data = cdc_obj2_clades_summ,
    inherit.aes = FALSE
  )

save_plot(
  cdc_obj2_ind_circulating_plot, "cdc-obj2-ind-circulating",
  width = 20, height = 25
)

cdc_obj2_gmt_circulating <- cdc_obj2_ind_circulating %>%
  group_by(study_year_lbl, site, timepoint, egg_lbl) %>%
  summarise(
    summarise_logmean(wmean_titre, out = "tibble"),
    .groups = "drop"
  ) %>%
  ggplot(aes(timepoint, mn)) +
  ggdark::dark_theme_bw(verbose = FALSE) +
  theme(
    panel.spacing = unit(0, "null"),
    strip.background = element_blank(),
    legend.position = "bottom",
    legend.box.spacing = unit(0, "null")
  ) +
  facet_grid(site ~ egg_lbl) +
  scale_y_log10("GMT against circulating", breaks = 5 * 2^(0:10)) +
  scale_color_discrete(
    "Study year",
    labels = as_labeller(function(x) str_replace(x, "Year ", ""))
  ) +
  scale_x_discrete("Timepoint", labels = label_cdc_timepoints) +
  geom_pointrange(
    aes(ymin = low, ymax = high, color = study_year_lbl),
    position = position_dodge(0.4)
  )

save_plot(
  cdc_obj2_gmt_circulating, "cdc-obj2-gmt-circulating",
  width = 15, height = 15
)
