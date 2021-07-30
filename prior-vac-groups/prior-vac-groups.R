# Annette's original email:

# We have talked about further exploring whether there is evidence that some
# sera may be focused on particular epitopes using reverse genetics viruses that
# each have one of the  five antigenic sites messed up by substituting with
# amino acids from H10 HA. I’m finding it almost impossible to choose sera that
# may be of interest by just looking at them and there vaccine histories.

# We know that people responded better with greater rises & somewhat higher post
# vaccine titres if they had fewer prior vaccinations, but can we investigate if
# there are any trends in terms of the regions/clades against which people
# responded best and the sequences of vaccines they received (not just number).
# For example, there is a small group of people who only received vaccine in
# 2015 and then again in 2018, and others who only received vaccine in 2014 and
# then 2018 (Y3 vaccine group). There is also a group in Y 1 who only received
# vaccine in 2013. In the first instance, could we at least see if the area
# under the curve for geometric titre ratio against 2007-2019 strains clusters
# in any way by sequence of vaccines. Not sure how many sequences there are, or
# how many there are with a reasonable n, but maybe not too many?? Objective 1
# has 357 people of whom 120 had some prior vaccination, but were not vaccinated
# every year. Other than area under the curve, are you able to extract
# information about where the peak of titre rise occurs in terms of virus or
# clade (perhaps doing this with an without excluding the vaccine strain)?

library(tidyverse)
library(kableExtra)

source("data/read_data.R")

cdc_obj1_hi <- read_data("cdc-obj1-hi") %>%
  inner_join(read_data("cdc-virus"), "virus_full") %>%
  mutate(virus_full = fct_reorder(virus_full, virus_year))
cdc_obj1_participants <- read_data("cdc-obj1-participant")
cdc_obj1_vax_hist <- read_data("cdc-obj1-vax-hist")

# NOTE(sen) What are the different vaccination history patterns? We are probably
# more interested in the specific year of vaccination rather than the offset
# from bleed. BUt then for any given prior vaccination patten, the year of the
# bleed is also important.

write_and_return <- function(data, name) {
  write_csv(data, glue::glue("prior-vac-groups/{name}.csv"))
  data
}

cdc_obj1_vax_hist_patterns <- cdc_obj1_vax_hist %>%
  inner_join(cdc_obj1_participants, "pid") %>%
  # NOTE(sen) Remove the multi-bleed year complexity for now, filter the data to
  # what we would have had if only the 2016 data collection year happened
  filter(recruitment_year == 1, year < 2016) %>%
  group_by(pid) %>%
  summarise(
    .groups = "drop",
    prior_vax_pattern = paste(year[status == 1], collapse = " "),
    prior_vax_pattern =
      if_else(prior_vax_pattern == "", "None", prior_vax_pattern),
    # prior_vax_pattern = fct_relevel(factor(prior_vax_pattern), "None")
  ) %>%
  mutate(
    prior_vax_pattern = factor(prior_vax_pattern) %>%
      fct_relevel("2011 2012 2013 2014 2015", "2013")
  )

cdc_obj1_vax_hist_patterns_counts <- cdc_obj1_vax_hist_patterns %>%
  count(prior_vax_pattern, name = "count") %>%
  arrange(desc(count))

cdc_obj1_vax_hist_patterns_counts %>%
  rename(Years = prior_vax_pattern, Count = count) %>%
  write_and_return("cdc-obj1-2016-prior-vax-patterns") %>%
  kbl(
    format = "latex",
    caption = "CDC objective 1 year 2016 prior vaccination patterns",
    booktabs = TRUE,
    label = "cdc-obj1-prior-vax-patterns"
  ) %>%
  write("prior-vac-groups/cdc-obj1-prior-vax-patterns.tex")

# NOTE(sen) Look at the landscapes for each pattern. Right now this is only the
# participants who were recruited in year 1

summarise_logmean <- function(titres) {
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
  tibble(logmn, logse, loglow, loghigh, mn, low, high)
}

cdc_obj1_prior_vac_hi <- cdc_obj1_hi %>%
  inner_join(cdc_obj1_vax_hist_patterns, "pid") %>%
  inner_join(cdc_obj1_vax_hist_patterns_counts, "prior_vax_pattern") %>%
  filter(count >= 9) %>%
  group_by(study_year, timepoint, virus_full, prior_vax_pattern) %>%
  summarise(.groups = "drop", summarise_logmean(titre))

# NOTE(sen) There will be missing values here if there is only one observation
# in a subset in the summary above
cdc_obj1_prior_vac_hi %>% filter(!complete.cases(.))

virus_labeller <- as_labeller(function(virus_names) {
  virus_names %>%
    tools::toTitleCase() %>%
    str_replace("^a/", "A/")
})

cdc_obj1_prior_vac_hi_plot <- cdc_obj1_prior_vac_hi %>%
  ggplot(aes(virus_full, mn, ymin = low, ymax = high, col = prior_vax_pattern)) +
  ggdark::dark_theme_bw(verbose = FALSE) +
  theme(
    legend.position = "bottom",
    panel.spacing = unit(0, "null"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(10, 10, 10, 20)
  ) +
  facet_wrap(
    ~ study_year + timepoint,
    ncol = 1, strip.position = "right",
    labeller = function(fac) {
      if ("study_year" %in% names(fac)) {
        fac$study_year <- paste("Study Year", fac$study_year)
      }
      if ("timepoint" %in% names(fac)) {
        fac$timepoint <- as.character(fac$timepoint)
      }
      fac
    }
  ) +
  scale_y_log10("Titre", breaks = 5 * 2^(0:10)) +
  scale_x_discrete("Virus", labels = virus_labeller) +
  scale_color_discrete("Prior vaccination pattern") +
  geom_pointrange(position = position_dodge(width = 1))

ggdark::ggsave_dark(
  "prior-vac-groups/cdc-obj1-prior-vac-hi.pdf",
  cdc_obj1_prior_vac_hi_plot,
  width = 30, height = 30, units = "cm"
)
