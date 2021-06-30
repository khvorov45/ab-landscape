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

# Script ======================================================================

cluster_probabilities <- read_csv("cluster/clusters.csv", col_types = cols())

cdc_obj1_cluster2_participants <- cluster_probabilities %>%
  filter(n_clusters == 2, chain == 1, cluster == 1) %>%
  mutate(
    cluster = if_else(prob_q50 > 0.5, "1", "2"),
  ) %>%
  select(pid, cluster)

cdc_obj1_participants <- cdc_obj1_cluster2_participants %>%
  inner_join(read_data("cdc-obj1-participant"), "pid") %>%
  mutate(
    prior_vacs_cat3 = case_when(
      prior_vacs == 0 ~ "0",
      prior_vacs == 5 ~ "5",
      TRUE ~ "1-4"
    )
  )

cdc_obj1_serology <- cdc_obj1_cluster2_participants %>%
  inner_join(read_data("cdc-obj1-hi"), "pid") %>%
  inner_join(read_data("cdc-virus"), "virus_full") %>%
  mutate(virus_full = fct_reorder(virus_full, virus_year))

cdc_obj1_participants %>%
  group_by(cluster) %>%
  summarise(
    count = n(),
    prior_vacs = summarise_factor(prior_vacs_cat3),
    age = summarise_numeric(age_first_bleed),
    gender = summarise_factor(gender),
    .groups = "drop"
  ) %>%
  kbl(
    format = "latex",
    caption = "Summaries of CDC objective 1 participants.
      Format: count (percentage); mean (sd)",
    booktabs = TRUE,
    col.names = c("Cluster", "Count", "Prior vaccinations", "Age", "Gender"),
    label = "cdc-obj1-participants-cluster2"
  ) %>%
  kable_styling(latex_options = "scale_down") %>%
  write("cluster/cdc-obj1-participants-cluster2.tex")

virus_labeller <- as_labeller(function(breaks) {
  breaks %>%
    tools::toTitleCase() %>%
    str_replace("^a/", "A/")
})

cdc_obj1_cluster_gmts <- cdc_obj1_serology %>%
  group_by(virus_full, cluster, timepoint) %>%
  summarise(summarise_logmean(titre, out = "tibble"), .groups = "drop") %>%
  ggplot(aes(virus_full, mn, color = cluster, shape = cluster)) +
  ggdark::dark_theme_bw(verbose = FALSE) +
  theme(
    legend.position = "bottom",
    legend.box.spacing = unit(0, "null"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.spacing = unit(0, "null"),
    strip.background = element_blank(),
    plot.margin = margin(5, 5, 5, 15)
  ) +
  facet_wrap(~timepoint, ncol = 1, strip.position = "right") +
  scale_y_log10("Titre", breaks = 5 * 2^(0:10)) +
  scale_x_discrete("Virus", labels = virus_labeller) +
  scale_color_discrete("Cluster") +
  scale_shape_discrete("Cluster") +
  geom_pointrange(aes(ymin = low, ymax = high), position = position_dodge(1))

ggdark::ggsave_dark(
  "cluster/cdc-obj1-cluster2-gmts.pdf",
  cdc_obj1_cluster_gmts,
  width = 20, height = 15, units = "cm"
)
