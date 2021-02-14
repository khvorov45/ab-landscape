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

year_keep_last_2 <- function(virus_name) {
  str_replace(virus_name, "/\\d{0,2}(\\d{2}[[:alpha:]]?)$", "/\\1")
}

standardise_full_virus_name <- function(virus_name) {
  virus_name %>%
    tolower() %>%
    year_keep_last_2() %>%
    str_replace("hongkong", "hong kong") %>%
    str_replace("infimh16/0019", "16-0019") %>%
    str_replace("/switz/", "/switzerland/") %>%
    str_replace("new castle", "newcastle")
}

# Convert short names to what's in the viruses table to look up long names to
# then look up coordinates
standardise_short_names <- function(names) {
  str_replace(names, "_", "") %>%
    str_replace("HKong", "Hkong") %>%
    str_replace("Philippine", "Phil") %>%
    str_replace("NethLand", "Neth") %>%
    str_replace("Tasmania", "Tas") %>%
    str_replace("NCastle", "Ncast") %>%
    str_replace("NCaled", "Ncal") %>%
    str_replace("Townsville", "T'ville") %>%
    str_replace("Brisbane", "Bris")
}

save_data <- function(data, name) {
  write_csv(data, glue::glue("data/{name}.csv"))
}

compare_vectors <- function(vec1, vec2, vec1_lbl = "in1", vec2_lbl = "in2") {
  in1 <- setdiff(vec1, vec2)
  in2 <- setdiff(vec2, vec1)
  if (length(in1) > length(in2)) {
    longer <- in1
    shorter <- in2
  } else {
    longer <- in2
    shorter <- in1
  }
  shorter_padded <- c(sort(shorter), rep(NA, length(longer) - length(shorter)))
  cmp <- tibble(sort(longer), shorter_padded)
  if (length(in1) > length(in2)) {
    colnames(cmp) <- c(vec1_lbl, vec2_lbl)
  } else {
    colnames(cmp) <- c(vec2_lbl, vec1_lbl)
  }
  cmp %>%
    select(!!rlang::sym(vec1_lbl), !!rlang::sym(vec2_lbl)) %>%
    filter(!is.na(!!rlang::sym(vec1_lbl)) | !is.na(!!rlang::sym(vec2_lbl)))
}

# Script ======================================================================

# The map ---------------------------------------------------------------------

agmap <- read_raw_csv("miniH3.coords", col_types = cols()) %>%
  rename(virus_full = Virus, virus_x_coord = X, virus_y_coord = Y) %>%
  mutate(
    virus_full = virus_full %>%
      str_replace_all("_", "/") %>%
      standardise_full_virus_name()
  )

save_data(agmap, "map")

# CDC -------------------------------------------------------------------------

# Viruses

cdc_viruses_raw_obj1 <- read_raw("cdc-obj1/Viruses")
cdc_viruses_raw_obj2 <- read_raw("cdc-obj2/Viruses")

fmt_cdc_viruses <- function(data) {
  data %>%
    select(
      virus_full = Virus_Name,
      virus_short = Short_name,
      virus_n = VirusN,
      virus_year = Year,
      clade = Clade,
      egg = Egg_Cell,
    ) %>%
    mutate(
      virus_full = standardise_full_virus_name(virus_full),
      virus_year = as.integer(virus_year),
      clade = replace_na(clade, "(Missing)"),
      egg = egg == 0,
      # Use the Egg_Cell column as the source of truth for egg status
      virus_short = str_replace(virus_short, "e$", "") %>%
        if_else(egg, paste0(., "e"), .),
      virus_full = str_replace(virus_full, "e$", "") %>%
        if_else(egg, paste0(., "e"), .),
    )
}

cdc_viruses_obj1 <- fmt_cdc_viruses(cdc_viruses_raw_obj1)
cdc_viruses_obj2 <- cdc_viruses_raw_obj2 %>%
  rename(Year = Virus_Year) %>%
  fmt_cdc_viruses()

# The differences are:
# - 'p' on the end
# - hanam vs hanoi
compare_vectors(cdc_viruses_obj1$virus_full, cdc_viruses_obj2$virus_full)
# The differences are:
# - 'p' on the end
compare_vectors(cdc_viruses_obj1$virus_short, cdc_viruses_obj2$virus_short)

# See if they would match with the map
# - 'p' on the end
compare_vectors(
  cdc_viruses_obj1$virus_full, agmap$virus_full, "ob1", "map"
) %>% print(n = 99)

# - 'p' on the end
# - hanam vs hanoi
compare_vectors(
  cdc_viruses_obj2$virus_full, agmap$virus_full, "ob2", "map"
) %>% print(n = 99)

cdc_viruses_obj1 %>% filter(!complete.cases(.))
cdc_viruses_obj2 %>% filter(!complete.cases(.))

# They may be united if the names are sorted out
save_data(cdc_viruses_obj1, "cdc-virus-obj1")
save_data(cdc_viruses_obj2, "cdc-virus-obj2")

# Participants for objective 1

cdc_hi_time_obj1 <- read_raw("cdc-obj1/HI_timecourse")

cdc_participants_obj1 <- select(
  cdc_hi_time_obj1,
  pid = PID, group = `Case/Control`, sex = Sex, yob = YoB
) %>%
  distinct(pid, .keep_all = TRUE) %>%
  mutate(
    group = str_replace(group, "Vax", "") %>% str_trim()
  ) %>%
  filter(pid != "KK38") # Annette

cdc_participants_obj1 %>% filter(!complete.cases(.))

save_data(cdc_participants_obj1, "cdc-participant-obj1")

# Participants for objective 2

cdc_hi_prior_vacc_obj2 <- read_raw("cdc-obj2/prior_vacc")

cdc_participants_obj2 <- cdc_hi_prior_vacc_obj2 %>%
  select(pid = study_id, sex = Sex, yob = yob) %>%
  distinct(pid, .keep_all = TRUE)

cdc_participants_obj2 %>% filter(!complete.cases(.))

save_data(cdc_participants_obj2, "cdc-participant-obj2")

# Prior vaccinations for objective 2

cdc_vacc_hist_obj2 <- cdc_hi_prior_vacc_obj2 %>%
  select(pid = study_id, contains("Vacc_")) %>%
  select(pid, matches("c$")) %>%
  pivot_longer(
    contains("Vacc_"),
    names_to = "vacc_year", values_to = "vacc_status"
  ) %>%
  mutate(
    vacc_year = str_replace(vacc_year, "Vacc_(\\d{4})c", "\\1") %>%
      as.integer()
  )

cdc_vacc_hist_obj2 %>% filter(!complete.cases(.))
unique(cdc_vacc_hist_obj2$vacc_year)
unique(cdc_vacc_hist_obj2$vacc_status)

compare_vectors(cdc_vacc_hist_obj2$pid, cdc_participants_obj2$pid)

save_data(cdc_vacc_hist_obj2, "cdc-vacc-hist-obj2")

# HI for objective 1

cdc_hi_obj1 <- cdc_hi_time_obj1 %>%
  select(pid = PID, contains("Titer"), virus = Virus_Name) %>%
  pivot_longer(
    contains("Titer"),
    names_to = "timepoint", values_to = "titre"
  ) %>%
  mutate(
    timepoint = recode(
      timepoint,
      "Pre.Titer" = "prevax",
      "Post.Titer" = "postvax",
      "PostS.Titer" = "postseas"
    )
  ) %>%
  filter(pid != "KK38")

cdc_hi_obj1 %>% filter(!complete.cases(.))

compare_vectors(cdc_hi_obj1$pid, cdc_participants_obj1$pid)

cdc_hi_obj1 %>%
  group_by(pid, virus) %>%
  filter(n() != 3)

save_data(cdc_hi_obj1, "cdc-hi-obj1")

cdc_hi_change <- cdc_hi_obj1 %>%
  filter(timepoint %in% c("prevax", "postvax")) %>%
  pivot_wider(names_from = timepoint, values_from = titre) %>%
  mutate(ratio = postvax / prevax)

save_data(cdc_hi_change, "cdc-hi-change-obj1")

# HI for objective 2

cdc_hi_obj2_raw <- read_raw("cdc-obj2/HIresult", guess_max = 1e5)

cdc_hi_obj2 <- cdc_hi_obj2_raw %>%
  select(
    pid = PID, study_year = SpecYearN, timepoint = DrawN, titre = Titer
  ) %>%
  mutate(
    timepoint = recode(
      timepoint,
      "1" = "prevax", "2" = "postvax", "3" = "postseas"
    )
  ) %>%
  filter(!is.na(titre))

cdc_hi_obj2 %>% filter(!complete.cases(.))

compare_vectors(cdc_hi_obj2$pid, cdc_participants_obj2$pid)

save_data(cdc_hi_obj2, "cdc-hi-obj2")

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
    Virus_Abbrv = col_character(),
    Short_Name = col_character()
  )
) %>%
  rename(
    pid = Subject_ID, timepoint = time, titre = Titer,
    virus_year = Year, virus = Virus_Abbrv, clade = Clade, prior_h3 = prior_H3,
    dob = DoBS, virus_short = Short_Name
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
    logtitre_mid = if_else(titre == 5, logtitre, logtitre + log(2) / 2),
    virus_short = standardise_short_names(virus_short)
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

# Attach the virus full name and coordinates
setdiff(unique(hi_hanam$virus_short), cdc_viruses_obj1$virus_short)

hi_hanam_virusmeta <- hi_hanam %>%
  left_join(
    select(cdc_viruses_obj1, virus_short, virus_full),
    by = "virus_short"
  ) %>%
  left_join(agmap, by = "virus_full")

save_csv(hi_hanam_virusmeta, "hi-hanam")

# RMH HCW study ---------------------------------------------------------------

hi_rmh_hcw <- read_raw_csv("HI_long", col_types = cols())

hi_rmh_hcw_reduced <- hi_rmh_hcw %>%
  mutate(
    virus_year = lubridate::year(IsolDate),
    infected = PID %in%
      c("RMH0052", "RMH0072", "RMH0044", "RMH0099", "RMH0117", "RMH0154"),
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
    titre = Titer, group, infected,
    age = Age, sex
  ) %>%
  mutate(
    timepoint = recode_3_timepoints(timepoint_num),
    age_lab = paste("Age:", age),
    logtitre = log(titre),
    logtitre_mid = if_else(titre == 5, logtitre, logtitre + log(2) / 2),
    virus = standardise_short_names(virus),
    egg = str_detect(virus, "e$"),
  )

setdiff(unique(hi_rmh_hcw_reduced$virus), cdc_viruses_obj1$virus_short)

hi_rmh_hcw_virusmeta <- hi_rmh_hcw_reduced %>%
  left_join(
    select(cdc_viruses_obj1, virus = virus_short, virus_full),
    by = "virus"
  ) %>%
  left_join(agmap, by = "virus_full")

save_csv(hi_rmh_hcw_virusmeta, "hi-rmh-hcw")
