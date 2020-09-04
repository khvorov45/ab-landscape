cat("extract data for analysis")

library(tidyverse)

data_dir <- "data"
data_raw_dir <- "data-raw"

# Functions ===================================================================

read_raw <- function(name, ...) {
  readxl::read_excel(file.path(data_raw_dir, paste0(name, ".xlsx")), ...)
}

read_raw_csv <- function(name, ...) {
  read_csv(file.path(data_raw_dir, paste0(name, ".csv")), ...)
}

equally_unique <- function(name1, name2, data) {
  # Alignment
  cond1 <- data %>%
    group_by(!!rlang::sym(name1)) %>%
    summarise(
      unique2 = length(unique(!!rlang::sym(name2))), .groups = "drop"
    ) %>%
    filter(unique2 != 1L)
  cond1 <- nrow(cond1) == 0L
  # Unique amount
  cond2 <- length(unique(data[[name1]])) == length(unique(data[[name2]]))
  cond1 & cond2
}

save_csv <- function(data, name) {
  write_csv(data, file.path(data_dir, paste0(name, ".csv")))
}

recode_3_timepoints <- function(timepoint_num) {
  recode(
    timepoint_num,
    "1" = "Pre-vax", "2" = "Post-vax", "3" = "Post-season",
  )
}

# Script ======================================================================

# The raw Israel data ---------------------------------------------------------

viruses_raw <- read_raw("Viruses")
dilutions_raw <- read_raw("VirusDiln_BT")
sera_raw <- read_raw("Sera List")
samples_raw <- read_raw("Obj1_Sample_info")
hi_raw <- read_raw("HI")

# Select the required columns

sera <- sera_raw %>%
  select(sample = Sample_ID, pid = PID, timepoint_num = Timepoint) %>%
  filter(timepoint_num %in% 1:3) %>%
  mutate(timepoint = recode_3_timepoints(timepoint_num))

hi <- select(
  hi_raw,
  sample = Sample_ID, virus = Virus, virus_n = VirusN, titre = Titer
)
viruses <- select(
  viruses_raw,
  virus = Virus_Name, virus_n = VirusN, virus_year = Year, clade = Clade
) %>%
  mutate(
    virus_year = as.integer(virus_year),
    clade = replace_na(clade, "(Missing)")
  )
participants <- select(
  samples_raw,
  pid = PID, group = `Case/Control`, sex = Sex, age = Age
) %>%
  distinct(pid, .keep_all = TRUE) %>%
  mutate(
    age_lab = paste("Age:", age),
    # It'd be good to know actual dates of birth
    year_of_birth = floor(2019.5 - age),
    group = recode(group, "Case" = "Infrequent", "Control" = "Frequent")
  )

# HI results

# See if the sample id's match
setdiff(hi$sample, sera$sample)
setdiff(sera$sample, hi$sample)

# See if virus names match
setdiff(viruses$virus, hi$virus) # %>% write_lines("in-viruses-not-in-hi.txt")
setdiff(hi$virus, viruses$virus) # %>% write_lines("in-hi-not-in-viruses.txt")

# See if virus numbers are as unique as names
equally_unique("virus", "virus_n", hi)
equally_unique("virus", "virus_n", viruses)

# Extra variables for HI results
hi_extra <- hi %>%
  mutate(
    logtitre = log(titre),
    logtitre_mid = if_else(titre == 5L, logtitre, logtitre + log(2) / 2)
  )

# Trust viruses for names
hi_no_virname <- select(hi_extra, -virus)

# See if virus numbers match
setdiff(viruses$virus_n, hi$virus_n)
setdiff(hi$virus_n, viruses$virus_n)

# See if pid's match
setdiff(sera$pid, participants$pid) # KK38 is Annette, exclude that
setdiff(participants$pid, sera$pid)

# Join the relevant bits of HI data
hi_full <- inner_join(hi_no_virname, sera, by = "sample") %>%
  inner_join(viruses, by = "virus_n") %>%
  inner_join(participants, by = "pid") %>%
  filter(!str_detect(sample, "KK38")) # Because Annette's

save_csv(hi_full, "hi")

# Israel objective 2 ----------------------------------------------------------

hi_2 <- read_raw("Obj2_timecourse_200904_complete")

hi_2_final <- hi_2 %>%
  select(
    pid = PID, sex = Sex, age = Age, virus_year = Virus_Year,
    n5y_prior_vacc = n5Y_prior_vacc,
    virus = Short_name,
    cluster,
    contains("L2HI"),
  ) %>%
  pivot_longer(
    contains("L2HI"),
    names_to = "timepoint_global", values_to = "titre"
  ) %>%
  mutate(
    titre = as.integer(round(2^titre, 0)),
    logtitre = log(titre),
    logtitre_mid = if_else(titre == 5L, logtitre, logtitre + log(2) / 2),
    timepoint_global = str_replace(
      timepoint_global, "^time\\s?(\\d)\\.L2HI$", "\\1"
    ) %>% as.integer(),
    timepoint_num = ((timepoint_global - 1) %% 3) + 1,
    timepoint = recode_3_timepoints(timepoint_num),
    age_lab = paste("Age:", age),
    clade = paste("cluster", cluster),
    study_year = ceiling(timepoint_global / 3),
    study_year_lab = glue::glue(
      "Year {study_year} 20{15 + study_year}/{16 + study_year}"
    ),
    n5y_prior_vacc_lab = paste0("Vac in past 5 years: ", n5y_prior_vacc)
  )

save_csv(hi_2_final, "hi-obj2")

# The Hanam dataset Annette gave me -------------------------------------------

hi_hanam <- read_csv(
  file.path(data_raw_dir, "HI_hanam.csv"),
  col_types = cols_only(
    Subject_ID = col_character(),
    time = col_integer(),
    Titer = col_integer(),
    Year = col_integer(),
    Clade = col_character(),
    prior_H3 = col_integer(),
    DoBS = col_date(),
    Virus_Abbrv = col_character()
  )
) %>%
  rename(
    pid = Subject_ID, timepoint = time, titre = Titer,
    virus_year = Year, virus = Virus_Abbrv, clade = Clade, prior_h3 = prior_H3,
    dob = DoBS
  ) %>%
  mutate(
    pid = str_replace(pid, "/", "-"),
    clade = replace_na(clade, "(Missing)"),
    prior_h3_lab = recode(
      prior_h3,
      "0" = "Prior H3: No", "1" = "Prior H3: Yes"
    ),
    year_of_birth = lubridate::year(dob),
    logtitre = log(titre),
    logtitre_mid = if_else(titre == 5, logtitre, logtitre + log(2) / 2)
  ) %>%
  filter(
    timepoint %in% 1:6,
    # Exclude all egg-grown except HK14e
    str_detect(virus, "e$", negate = TRUE) | virus == "HK14e"
  ) %>%
  mutate(
    # Recode timepoint
    timepoint_num = timepoint,
    timepoint = recode(
      timepoint_num,
      "1" = "BL", "2" = "d4", "3" = "d7", "4" = "d14", "5" = "d21", "6" = "d280"
    )
  )

save_csv(hi_hanam, "hi-hanam")

# RMH HCW study ---------------------------------------------------------------

hi_rmh_hcw <- read_raw_csv("HI_long", col_types = cols())

hi_rmh_hcw_reduced <- hi_rmh_hcw %>%
  mutate(
    virus_year = lubridate::year(IsolDate),
    group = case_when(
      Vacc5Yn >= 3 ~ "Frequent",
      Vacc5Yn == 2 ~ "Moderate",
      Vacc5Yn <= 1 ~ "Infrequent"
    ),
    year_of_birth = lubridate::year(dob)
  ) %>%
  select(
    pid = PID, timepoint_num = TimeN, virus = Short_Name, clade = Clade,
    year_of_birth, virus_year,
    titre = Titer, group,
    age = Age, sex
  ) %>%
  mutate(
    timepoint = recode_3_timepoints(timepoint_num),
    age_lab = paste("Age:", age),
    logtitre = log(titre),
    logtitre_mid = if_else(titre == 5, logtitre, logtitre + log(2) / 2)
  )

save_csv(hi_rmh_hcw_reduced, "hi-rmh-hcw")
