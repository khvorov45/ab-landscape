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
  if (length(nums) == 1) {
    return(glue::glue("{format_decimal(nums)}"))
  }
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

summarise_baseline <- function(data, ...) {
  data %>%
    group_by(...) %>%
    summarise(
      Count = n() %>% as.character(),
      gender_prop = summarise_factor(recode(gender, "Female" = "F", "Male" = "M")),
      age_summ = summarise_numeric(age_first_bleed),
      .groups = "drop"
    ) %>%
    pivot_longer(
      -c(...),
      names_to = "summary_type", values_to = "summary"
    ) %>%
    pivot_wider(names_from = "prior_vacs", values_from = "summary") %>%
    mutate(
      summary_type = recode(
        summary_type,
        "gender_prop" = "Gender",
        "age_summ" = "Age (years) at first bleed"
      )
    )
}

# Script ======================================================================

# CDC Objective 1

cdc_obj1_participant <- read_data("cdc-obj1-participant") %>%
  mutate(
    recruitment_year_lbl = recode(
      recruitment_year,
      "1" = "1 (2016)", "2" = "2 (2017)", "3" = "3 (2018)"
    ),
    prior_vacs2 = recode(
      prior_vacs,
      "0" = "0", "1" = "1-3", "2" = "1-3", "3" = "1-3", "5" = "5"
    ),
  )

summarise_baseline(
  cdc_obj1_participant, prior_vacs, site, recruitment_year_lbl
) %>%
  bind_rows(summarise_baseline(cdc_obj1_participant, prior_vacs, site) %>%
    mutate(recruitment_year_lbl = "Overall")) %>%
  bind_rows(summarise_baseline(cdc_obj1_participant, prior_vacs) %>%
    mutate(recruitment_year_lbl = "Overalll", site = "Both")) %>%
  mutate(site = fct_relevel(site, "Both", after = Inf)) %>%
  arrange(site, recruitment_year_lbl) %>%
  save_data("cdc-participant-summary-obj1") %>%
  kbl(
    format = "latex",
    caption = "Summaries of CDC objective 1 participants.
      Format: count (percentage); mean (sd)",
    booktabs = TRUE,
    col.names = c("Site", "Year", "", colnames(.)[-(1:3)]),
    label = "cdc-participant-summary-obj1"
  ) %>%
  add_header_above(
    c(" " = 3, "Vaccinations in 5 years before recruitment" = 5)
  ) %>%
  column_spec(4:8, width = "2.2cm", latex_valign = "m") %>%
  column_spec(1, width = "1cm") %>%
  column_spec(2, width = "1.5cm") %>%
  column_spec(3, width = "2.5cm", latex_valign = "m") %>%
  collapse_rows(
    columns = 1:2, latex_hline = "custom", custom_latex_hline = 2
  ) %>%
  kable_styling(latex_options = "scale_down") %>%
  save_table("cdc-participant-summary-obj1")

cdc_viruses <- read_data("cdc-virus")
cdc_vaccine <- read_data("cdc-vaccine")

# Viruses in each clade
cdc_viruses %>%
  filter(!clade %in% c("1", "2", "(missing)")) %>%
  mutate(egg_lbl = if_else(egg, "Egg", "Cell")) %>%
  select(Clade = clade, Virus = virus_full, Type = egg_lbl) %>%
  arrange(Clade) %>%
  kbl(
    format = "latex",
    caption = "Viruses in each clade.",
    booktabs = TRUE,
    label = "cdc-clade-viruses"
  ) %>%
  collapse_rows(columns = 1, latex_hline = "major") %>%
  save_table("cdc-clade-viruses")

cdc_circulating_year <- function(study_year, site) {
  case_when(
    study_year == 1 & site == "Peru" ~ 2016.75,
    study_year == 1 & site == "Israel" ~ 2017.25,
    study_year == 2 & site == "Peru" ~ 2017.75,
    study_year == 2 & site == "Israel" ~ 2018.25,
    study_year == 3 & site == "Peru" ~ 2018.75,
    study_year == 3 & site == "Israel" ~ 2019.25,
    TRUE ~ NA_real_
  )
}

cdc_obj1_hi <- read_data("cdc-obj1-hi") %>%
  inner_join(cdc_obj1_participant, by = "pid") %>%
  inner_join(cdc_viruses, "virus_full") %>%
  mutate(
    vaccine_strain = virus_full %in% cdc_vaccine$virus_full,
    egg_lbl = if_else(egg, "Egg", "Cell"),
    circulating_year = cdc_circulating_year(study_year, site)
  )

# Titre plot - the usual stuff, variable at prevax, rises at postvax, stays
# the same/drops slightly to postseas.
cdc_obj1_hi %>%
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
cdc_obj1_timepoint_gmts <- cdc_obj1_hi %>%
  group_by(timepoint, virus_full, prior_vacs2, site) %>%
  summarise(summarise_logmean(titre, out = "tibble"), .groups = "drop") %>%
  ggplot(
    aes(
      virus_full, mn,
      color = prior_vacs2, shape = prior_vacs2
    )
  ) +
  ggdark::dark_theme_bw(verbose = FALSE) +
  theme(
    legend.position = "bottom",
    legend.box.spacing = unit(0, "null"),
    strip.placement = "right",
    panel.spacing = unit(0.5, "lines"),
    strip.background = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = -45, hjust = 0),
    plot.margin = margin(10, 40, 10, 10)
  ) +
  facet_wrap(
    ~ site + timepoint,
    ncol = 1, strip.position = "right"
  ) +
  scale_y_log10("GMT (95% CI)", breaks = 5 * 2^(0:15)) +
  scale_x_discrete("Virus") +
  scale_color_discrete("Prior vaccinations") +
  scale_shape_discrete("Prior vaccinations") +
  geom_vline(
    xintercept = cdc_obj1_hi$virus_full[cdc_obj1_hi$vaccine_strain],
    size = 7, alpha = 0.2
  ) +
  geom_pointrange(
    aes(ymin = low, ymax = high),
    position = position_dodge(width = 0.8)
  )

save_plot(
  cdc_obj1_timepoint_gmts, "cdc-obj1-timepoint-gmts",
  width = 20, height = 25
)

# Differences between timepoints
# Highlight vaccine strain here as well
timepoint_diffs <- cdc_obj1_hi %>%
  select(-bleed_date) %>%
  filter(timepoint %in% c("Pre-vax", "Post-vax")) %>%
  pivot_wider(names_from = "timepoint", values_from = "titre") %>%
  mutate(ratio = `Post-vax` / `Pre-vax`) %>%
  group_by(virus_full, prior_vacs2) %>%
  summarise(summarise_logmean(ratio, out = "tibble"), .groups = "drop") %>%
  ggplot(aes(virus_full, mn, color = prior_vacs2, shape = prior_vacs2)) +
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
    xintercept = cdc_obj1_hi$virus_full[cdc_obj1_hi$vaccine_strain],
    size = 4, alpha = 0.2
  ) +
  geom_pointrange(
    aes(ymin = low, ymax = high),
    position = position_dodge(width = 0.5)
  ) +
  scale_y_log10("Post/Pre vax ratio (95% CI)", breaks = 1:10) +
  scale_x_discrete("Virus") +
  scale_color_discrete("Prior vaccinations") +
  scale_shape_discrete("Prior vaccinations")

save_plot(
  timepoint_diffs, "cdc-obj1-timepoint-diffs",
  width = 20, height = 15
)

# CDC Objective 2

cdc_obj2_participants <- read_data("cdc-obj2-participant")

summarise_baseline(
  cdc_obj2_participants, prior_vacs, site
) %>%
  bind_rows(summarise_baseline(cdc_obj2_participants, prior_vacs) %>%
    mutate(site = "Both")) %>%
  save_data("cdc-participant-summary-obj2") %>%
  kbl(
    format = "latex",
    caption = "Summaries of CDC objective 2 participants.
      Format: count (percentage); mean (sd)",
    booktabs = TRUE,
    col.names = c("Site", "", colnames(.)[-(1:2)]),
    label = "cdc-participant-summary-obj2"
  ) %>%
  add_header_above(
    c(" " = 2, "Vaccinations in 5 years before recruitment" = 3)
  ) %>%
  column_spec(1, width = "1cm", latex_valign = "m") %>%
  column_spec(2, width = "2.5cm", latex_valign = "m") %>%
  column_spec(3:5, width = "2.2cm", latex_valign = "m") %>%
  collapse_rows(
    columns = 1, latex_hline = "custom", custom_latex_hline = 1
  ) %>%
  save_table("cdc-participant-summary-obj2")

# Titres
cdc_obj2_hi <- read_data("cdc-obj2-hi") %>%
  inner_join(cdc_obj2_participants, by = "pid") %>%
  inner_join(cdc_viruses, "virus_full") %>%
  mutate(
    vaccine_strain = virus_full %in% cdc_vaccine$virus_full,
    egg_lbl = if_else(egg, "Egg", "Cell"),
    study_year_lbl = recode(
      study_year,
      "1" = "1 (2016)", "2" = "2 (2017)", "3" = "3 (2018)"
    ),
    circulating_year = cdc_circulating_year(study_year, site)
  )

# Titre summaries
cdc_obj2_gmts <- cdc_obj2_hi %>%
  group_by(timepoint, virus_full, vaccine_strain, study_year_lbl, site) %>%
  summarise(summarise_logmean(titre, out = "tibble"), .groups = "drop")

plot_obj2_gmts <- function(data, group_var, group_var_lab) {
  group_var_q <- enquo(group_var)
  data %>%
    ggplot(aes(virus_full, mn, color = !!group_var_q, shape = !!group_var_q)) +
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
    scale_y_log10("GMT (95% CI)", breaks = 5 * 2^(0:15)) +
    scale_x_discrete("Virus") +
    geom_vline(
      xintercept = data$virus_full[data$vaccine_strain],
      size = 6, alpha = 0.2
    ) +
    geom_pointrange(
      aes(ymin = low, ymax = high),
      position = position_dodge(width = 0.5)
    ) +
    scale_color_discrete(group_var_lab) +
    scale_shape_discrete(group_var_lab)
}

cdc_obj2_gmts_plot1 <- cdc_obj2_gmts %>%
  plot_obj2_gmts(timepoint, "Timepoint") +
  facet_wrap(
    ~ site + study_year_lbl,
    ncol = 1, strip.position = "right",
    labeller = function(labels) {
      mutate(labels, study_year_lbl = paste0("Year ", study_year_lbl))
    }
  )

save_plot(
  cdc_obj2_gmts_plot1, "cdc-obj2-gmts-1",
  width = 20, height = 25
)

cdc_obj2_gmts_plot2 <- cdc_obj2_gmts %>%
  plot_obj2_gmts(study_year_lbl, "Study year") +
  facet_wrap(
    ~ site + timepoint,
    ncol = 1, strip.position = "right",
  )

save_plot(
  cdc_obj2_gmts_plot2, "cdc-obj2-gmts-2",
  width = 20, height = 25
)

# Titre changes
cdc_obj2_hi_wide <- cdc_obj2_hi %>%
  select(-bleed_date) %>%
  pivot_wider(names_from = "timepoint", values_from = "titre") %>%
  mutate(vax_resp = `Post-vax` / `Pre-vax`)

cdc_obj2_vax_resp_virus_plot <- cdc_obj2_hi_wide %>%
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

cdc_obj1_bleed_dates <- cdc_obj1_hi %>%
  ggplot(aes(bleed_date, pid, color = timepoint)) +
  ggdark::dark_theme_bw(verbose = FALSE) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "bottom",
    legend.box.spacing = unit(0, "null")
  ) +
  scale_x_date("Bleed date", breaks = "month") +
  scale_color_discrete("Timepoint") +
  geom_point() +
  geom_hline(yintercept = 177.5) +
  annotate(
    geom = "text",
    x = lubridate::ymd("2019-04-01"),
    y = 230, label = "Peru"
  ) +
  annotate(
    geom = "text",
    x = lubridate::ymd("2016-06-01"),
    y = 50, label = "Israel"
  )

save_plot(cdc_obj1_bleed_dates, "cdc-obj1-bleed-dates", width = 15, height = 20)

cdc_obj1_clades_summ <- cdc_obj1_hi %>%
  filter(!clade %in% c("(missing)", "1", "2")) %>%
  select(site, egg_lbl, circulating_year, study_year, prior_vacs2, clade) %>%
  distinct() %>%
  inner_join(cdc_clade_frequencies, c("clade", "circulating_year" = "year")) %>%
  filter(freq > 0) %>%
  group_by(site, egg_lbl, circulating_year, study_year, prior_vacs2) %>%
  summarise(clades = list(clade), freq_total = sum(freq), .groups = "drop")

# This is supposed to show titres against what was circulating at the time
cdc_obj1_hi_against_circulating <- cdc_obj1_hi %>%
  filter(clade %in% cdc_clade_frequencies$clade) %>%
  group_by(
    pid, timepoint, egg_lbl, site, circulating_year, study_year, prior_vacs2,
    clade
  ) %>%
  summarise(
    mean_log_titre = mean(log(titre)),
    .groups = "drop"
  ) %>%
  inner_join(cdc_clade_frequencies, c("clade", "circulating_year" = "year")) %>%
  group_by(
    pid, timepoint, egg_lbl, site, circulating_year, study_year, prior_vacs2
  ) %>%
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
    expand = expansion(0.1)
  ) +
  facet_grid(
    site + egg_lbl + prior_vacs2 ~ study_year,
    labeller = function(labels) {
      if ("study_year" %in% colnames(labels)) {
        mutate(labels, study_year = paste0("Year ", study_year))
      } else if ("prior_vacs2" %in% colnames(labels)) {
        mutate(labels, prior_vacs2 = paste0(prior_vacs2, " prior"))
      } else {
        labels
      }
    }
  ) +
  geom_line(alpha = 0.3) +
  geom_point() +
  geom_text(
    aes(1, 640, label = paste0(signif(freq_total * 100, 2), "%")),
    data = cdc_obj1_clades_summ,
    inherit.aes = FALSE
  )

save_plot(
  cdc_obj1_ind_av_circulating, "cdc-obj1-ind-av-circulating",
  width = 18, height = 25
)

cdc_obj1_gmt_circulating <- cdc_obj1_hi_against_circulating %>%
  group_by(prior_vacs2, timepoint, site, egg_lbl) %>%
  summarise(
    summarise_logmean(wmean_titre, out = "tibble"),
    .groups = "drop"
  ) %>%
  ggplot(aes(timepoint, mn, color = prior_vacs2)) +
  ggdark::dark_theme_bw(verbose = FALSE) +
  theme(
    legend.position = "bottom",
    legend.box.spacing = unit(0, "null"),
    strip.background = element_blank(),
    panel.spacing = unit(0, "null")
  ) +
  facet_grid(site ~ egg_lbl) +
  scale_y_log10(
    "GMT against circulating strains",
    breaks = 5 * 2^(0:10),
  ) +
  scale_x_discrete("Timepoint") +
  scale_color_discrete("Prior vaccinations") +
  geom_pointrange(
    aes(ymin = low, ymax = high),
    position = position_dodge(width = 0.4)
  )

save_plot(
  cdc_obj1_gmt_circulating, "cdc-obj1-gmt-circulating",
  width = 15, height = 15
)

# See what the seasons were in objective 2
cdc_obj2_bleed_dates_plot <- cdc_obj2_hi %>%
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
  scale_color_discrete("Timepoint") +
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

cdc_obj2_clades_summ <- cdc_obj2_hi %>%
  filter(!clade %in% c("(missing)", "1", "2")) %>%
  select(egg_lbl, site, circulating_year, study_year, clade) %>%
  distinct() %>%
  inner_join(cdc_clade_frequencies, c("clade", "circulating_year" = "year")) %>%
  filter(freq > 0) %>%
  group_by(egg_lbl, site, circulating_year, study_year) %>%
  summarise(clades = list(clade), freq_total = sum(freq), .groups = "drop")

cdc_obj2_ind_circulating <- cdc_obj2_hi %>%
  group_by(
    pid, site, study_year, timepoint, egg, egg_lbl,
    circulating_year, clade
  ) %>%
  summarise(logtitre_mean = mean(log(titre)), .groups = "drop") %>%
  inner_join(cdc_clade_frequencies, c("clade", "circulating_year" = "year")) %>%
  group_by(
    pid, site, study_year, timepoint, egg, egg_lbl,
    circulating_year
  ) %>%
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
  facet_grid(site + egg_lbl ~ study_year, labeller = function(labels) {
    if ("study_year" %in% names(labels)) {
      mutate(labels, study_year = paste0("Year ", study_year))
    } else {
      labels
    }
  }) +
  scale_x_discrete("Timepoint") +
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
  width = 20, height = 20
)

cdc_obj2_gmt_circulating <- cdc_obj2_ind_circulating %>%
  group_by(study_year, site, timepoint, egg_lbl) %>%
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
  ) +
  scale_x_discrete("Timepoint") +
  geom_pointrange(
    aes(ymin = low, ymax = high, color = as.factor(study_year)),
    position = position_dodge(0.4)
  )

save_plot(
  cdc_obj2_gmt_circulating, "cdc-obj2-gmt-circulating",
  width = 15, height = 15
)

# Look into Bilthoven viruses more closely ====================================

yob_cutoffs <- c(-Inf, 1968, 1980, Inf)
yob_labels <- c("1951-1967", "1968-1979", "1980-1991")

cdc_obj1_yobs_plot <- cdc_obj1_participant %>%
  ggplot(aes(lubridate::year(dob))) +
  ggdark::dark_theme_bw(verbose = FALSE) +
  theme(
    legend.position = "bottom",
    legend.box.spacing = unit(0, "null")
  ) +
  scale_y_continuous("Count", expand = expansion(c(0, 0.1))) +
  scale_x_continuous("Year of birth") +
  scale_fill_discrete("Prior vaccinations") +
  geom_histogram(aes(fill = prior_vacs2), color = "black", binwidth = 1)

save_plot(cdc_obj1_yobs_plot, "cdc-obj1-yobs", width = 14, height = 10)

cdc_obj1_ind_bilthoven <- cdc_obj1_hi %>%
  filter(str_detect(virus_full, "bilthoven")) %>%
  mutate(yob = lubridate::year(dob)) %>%
  group_by(pid, timepoint, prior_vacs2, yob) %>%
  summarise(bilthoven_mean = exp(mean(log(titre))), .groups = "drop") %>%
  mutate(
    yob_cat = cut(yob, yob_cutoffs, right = FALSE),
    x_position = as.integer(yob_cat) + 0.2 * (as.integer(timepoint) - 2),
  )

cdc_obj1_ind_bilthoven_plot <- cdc_obj1_ind_bilthoven %>%
  ggplot(aes(x_position, bilthoven_mean, color = timepoint, group = pid)) +
  ggdark::dark_theme_bw(verbose = FALSE) +
  theme(
    legend.position = "bottom", legend.box.spacing = unit(0, "null"),
    panel.spacing = unit(0, "null"),
    strip.background = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  facet_wrap(
    ~prior_vacs2,
    labeller = as_labeller(function(x) paste0(x, " prior"))
  ) +
  scale_y_log10("Average Bilthoven titre", breaks = 5 * 2^(0:10)) +
  scale_x_continuous(
    "Year of birth",
    breaks = 1:3,
    labels = yob_labels
  ) +
  scale_color_discrete("Timepoint") +
  geom_line() +
  geom_point() +
  geom_text(
    aes(x_position, 5, label = n),
    data = cdc_obj1_ind_bilthoven %>%
      group_by(yob_cat, prior_vacs2) %>%
      summarise(n = length(unique(pid)), .groups = "drop") %>%
      mutate(x_position = as.integer(yob_cat)),
    inherit.aes = FALSE
  )

save_plot(
  cdc_obj1_ind_bilthoven_plot, "cdc-obj1-ind-bilthoven",
  width = 15, height = 12
)

cdc_obj1_gmt_bilthoven <- cdc_obj1_ind_bilthoven %>%
  group_by(timepoint, prior_vacs2, yob_cat) %>%
  summarise(summarise_logmean(bilthoven_mean, out = "tibble"), .groups = "drop") %>%
  ggplot(aes(yob_cat, mn, color = timepoint)) +
  ggdark::dark_theme_bw(verbose = FALSE) +
  theme(
    panel.spacing = unit(0, "null"),
    strip.background = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_y_log10("Mean of averaged Bilthoven titres", breaks = 5 * 2^(0:10)) +
  scale_x_discrete("Year of birth", labels = yob_labels) +
  scale_color_discrete("Timepoint") +
  facet_wrap(
    ~prior_vacs2,
    nrow = 1,
    labeller = as_labeller(function(x) paste0(x, " prior"))
  ) +
  geom_pointrange(
    aes(ymin = low, ymax = high),
    position = position_dodge(width = 0.5)
  )

save_plot(
  cdc_obj1_gmt_bilthoven, "cdc-obj1-gmt-bilthoven",
  width = 15, height = 8
)
