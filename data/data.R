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
    str_replace("new castle", "newcastle") %>%
    str_replace("p$", "") %>%
    str_replace("hanam", "hanoi")
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
    str_replace("Brisbane", "Bris") %>%
    str_replace("p$", "")
}

standardise_clades <- function(clades) {
  clades %>%
    tolower() %>%
    str_replace_all("\\.", "") %>%
    str_replace_all(" ", "")
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

parse_access_date <- function(date_strings, century = "19") {
  date_strings %>%
    str_replace("^(\\d\\d/\\d\\d)/(\\d\\d) ", paste0("\\1/", century, "\\2 ")) %>%
    lubridate::mdy_hms() %>%
    lubridate::as_date()
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

# Viruses ---------------------------

cdc_viruses_raw_obj1 <- read_raw_csv("cdc-obj1/Viruses", col_types = cols())
cdc_viruses_raw_obj2 <- read_raw_csv("cdc-obj2/Viruses", col_types = cols())
cdc_viruses_raw_obj3 <- read_raw_csv("cdc-obj3/Viruses", col_types = cols())
cdc_viruses_raw_obj4 <- read_raw_csv("cdc-obj4/Viruses", col_types = cols())

fmt_cdc_viruses <- function(data) {
  data %>%
    select(
      virus_full = Virus_Name,
      virus_short = Short_name,
      virus_n = VirusN,
      virus_year = Virus_Year,
      clade = Clade,
      egg = Egg_Cell,
    ) %>%
    mutate(
      virus_full = standardise_full_virus_name(virus_full),
      virus_short = standardise_short_names(virus_short),
      virus_year = as.integer(virus_year),
      clade = replace_na(clade, "(Missing)"),
      egg = egg == 0,
      # Use the Egg_Cell column as the source of truth for egg status
      virus_short = str_replace(virus_short, "e$", "") %>%
        if_else(egg, paste0(., "e"), .),
      virus_full = str_replace(virus_full, "e$", "") %>%
        if_else(egg, paste0(., "e"), .),
      clade = standardise_clades(clade)
    )
}

cdc_viruses_obj1 <- fmt_cdc_viruses(cdc_viruses_raw_obj1)
cdc_viruses_obj2 <- fmt_cdc_viruses(cdc_viruses_raw_obj2)
cdc_viruses_obj3 <- fmt_cdc_viruses(cdc_viruses_raw_obj3)
cdc_viruses_obj4 <- fmt_cdc_viruses(cdc_viruses_raw_obj4)

compare_vectors(cdc_viruses_obj1$virus_full, cdc_viruses_obj2$virus_full)
compare_vectors(cdc_viruses_obj1$virus_full, cdc_viruses_obj3$virus_full)
compare_vectors(cdc_viruses_obj1$virus_full, cdc_viruses_obj4$virus_full)

compare_vectors(cdc_viruses_obj1$virus_short, cdc_viruses_obj2$virus_short)
compare_vectors(cdc_viruses_obj1$virus_short, cdc_viruses_obj3$virus_short)
compare_vectors(cdc_viruses_obj1$virus_short, cdc_viruses_obj4$virus_short)

# See if clades match
compare_vectors(cdc_viruses_obj1$clade, cdc_viruses_obj2$clade)
compare_vectors(cdc_viruses_obj1$clade, cdc_viruses_obj3$clade)
compare_vectors(cdc_viruses_obj1$clade, cdc_viruses_obj4$clade)

cdc_viruses_obj1 %>% filter(!complete.cases(.))
cdc_viruses_obj2 %>% filter(!complete.cases(.))
cdc_viruses_obj3 %>% filter(!complete.cases(.))
cdc_viruses_obj4 %>% filter(!complete.cases(.))

# They are the same table, so save one
cdc_viruses <- cdc_viruses_obj1

# See how the viruses match to the map
compare_vectors(
  cdc_viruses$virus_full, agmap$virus_full, "ob1", "map"
) %>% print(n = 99)

save_data(cdc_viruses, "cdc-virus")

# Vaccine viruses ---------------------------

cdc_vaccine <- tibble(
  virus_full = c(
    "a/hong kong/4801/14e", "a/hong kong/4801/14e", "a/singapore/16-0019/16e"
  ),
  study_year = c(1, 2, 3)
)
cdc_vaccine$virus_full %in% cdc_viruses_obj2$virus_full

save_data(cdc_vaccine, "cdc-vaccine")

# Look for virus/clade frequencies ------------------------

nextrain_freqs_raw <- read_raw_csv(
  "nexstrain-virus-frequencies",
  col_types = cols()
)

nextrain_freqs <- nextrain_freqs_raw %>%
  mutate(
    clade = standardise_clades(clade),
    name = tolower(name)
  )

# We have a lot of viruses that are not in that table
setdiff(cdc_viruses_obj1$virus_full, nextrain_freqs$name)

clade_frequencies <- nextrain_freqs %>%
  filter(clade != "unassigned", year >= 2016, year < 2020) %>%
  group_by(year, clade) %>%
  summarise(freq = sum(freq), .groups = "drop") %>%
  mutate(
    # Got this by looking through viruses and checking our clade and
    # Nextstrain's
    clade = recode(clade, "a1" = "3c2a1", "a2/re" = "3c2a2", "3c" = "3c1") %>%
      str_replace("^a1b/", "3c2a1b+")
  )

# Clades 1 and 2 correspond to 'unnassigned'
compare_vectors(cdc_viruses_obj1$clade, clade_frequencies$clade)

save_data(clade_frequencies, "cdc-clade-frequencies")

# Participants for objective 1 -------------------

cdc_obj1_participants_raw <- read_raw_csv(
  "cdc-obj1/vax_history",
  col_types = cols()
)

# Check that there is one row per participant
cdc_obj1_participants_raw %>%
  group_by(study_id) %>%
  filter(n() > 1)

cdc_obj1_participants <- cdc_obj1_participants_raw %>%
  mutate(
    dob_og = dob,
    dob = if_else(
      is.na(dob),
      lubridate::ymd(paste0(YoB, "-07-01")),
      dob %>% parse_access_date("19")
    ),
  ) %>%
  select(pid = study_id, site = Site, gender = Sex, dob)

# Missing values
cdc_obj1_participants %>% filter(!complete.cases(.))

# Vaccination history -----------------------------------------

cdc_obj1_vax_hist <- cdc_obj1_participants_raw %>%
  select(pid = study_id, contains("Vac")) %>%
  pivot_longer(-pid, names_to = "year", values_to = "status") %>%
  mutate(
    year = str_replace(year, "Vac", "") %>% as.numeric()
  )

# HI titres --------------------------------------

cdc_obj1_hi_raw <- read_raw_csv(
  "cdc-obj1/HI",
  col_types = cols(), guess_max = 1e4
)

cdc_obj1_hi <- cdc_obj1_hi_raw %>%
  mutate(
    sample_id = toupper(Sample_ID) %>% str_replace("NO SAMPLE - ", "")
  ) %>%
  select(sample_id, virus_n = VirusN, titre = Titer) %>%
  filter(!is.na(titre), !str_detect(sample_id, "KK38 2020 D"))

cdc_obj1_hi %>% filter(!complete.cases(.))

# Dates -----------------------------------------------

cdc_obj1_dates_raw1 <- read_raw_csv("cdc-obj1/Obj1Dates", col_types = cols())
cdc_obj1_dates_raw2 <- read_raw_csv("cdc-obj1/Obj1_Dates", col_types = cols())

cdc_obj1_dates1 <- cdc_obj1_dates_raw1 %>%
  mutate(
    timepoint = recode_3_timepoints(Blood_DrawN),
    bleed_date = Blood_Draw_date %>% parse_access_date("20")
  ) %>%
  select(
    pid = PID, sample_id = Specimen_ID, site = Site, timepoint,
    bleed_date, study_year = Specimen_Year
  )

cdc_obj1_dates2 <- cdc_obj1_dates_raw2 %>%
  mutate(
    timepoint = recode_3_timepoints(Blood_Draw),
    bleed_date = Blood_Draw_date %>% parse_access_date("20")
  ) %>%
  select(
    pid = study_id, sample_id = Specimen_ID, site = Site, timepoint,
    bleed_date, study_year = Specimen_Year
  )

dates_comp_result <- compare_vectors(
  cdc_obj1_dates1$sample_id, cdc_obj1_dates2$sample_id
)

cdc_obj1_dates <- cdc_obj1_dates1 %>%
  # Add to 1 what's in 2 but not in 1
  bind_rows(cdc_obj1_dates2 %>% filter(sample_id %in% dates_comp_result$in2)) %>%
  filter(complete.cases(.))

# Now add the required columns to the various tables -------------

# HI needs bleed dates

# This shows the sample ids for which we have hi data but not dates
cdc_obj1_hi_no_dates <- compare_vectors(
  cdc_obj1_dates$sample_id, cdc_obj1_hi$sample_id, "dates", "hi"
) %>%
  filter(!is.na(hi)) %>%
  select(hi)

compare_vectors(cdc_obj1_hi$virus_n, cdc_viruses$virus_n)

cdc_obj1_hi_extra <- cdc_obj1_hi %>%
  inner_join(cdc_obj1_dates, "sample_id") %>%
  inner_join(cdc_viruses, "virus_n") %>%
  select(pid, study_year, timepoint, bleed_date, virus_full, titre)

# One year's worth of data for each participant
cdc_obj1_hi_extra %>%
  count(pid, timepoint, virus_full) %>%
  filter(n > 1)

# First bleeds
cdc_obj1_first_bleeds <- cdc_obj1_dates %>%
  group_by(pid) %>%
  summarise(
    first_bleed = min(bleed_date),
    recruitment_year = min(study_year),
    .groups = "drop"
  )

cdc_obj1_first_bleeds %>% filter(!complete.cases(.))

# Prior vaccinations
compare_vectors(cdc_obj1_vax_hist$pid, cdc_obj1_participants$pid)
cdc_obj1_prior_vacs <- cdc_obj1_vax_hist %>%
  inner_join(cdc_obj1_first_bleeds, "pid") %>%
  # Look only at 5 years prior to recruitment
  # E.g. recruited in 2016 (study_year 1) -> 2011-2015 inclusive
  filter(
    year < recruitment_year + 2015, year >= recruitment_year + 2015 - 5,
    !is.na(status)
  ) %>%
  group_by(pid) %>%
  summarise(prior_vacs = sum(status), .groups = "drop")

# Time difference b/w prevax postvax postseason
cdc_obj1_pre_postvax_time <- cdc_obj1_dates %>%
  select(-sample_id) %>%
  # There are duplicate rows
  distinct() %>%
  pivot_wider(names_from = "timepoint", values_from = "bleed_date") %>%
  mutate(
    pre_post_vax_days = (`Post-vax` - `Pre-vax`) / lubridate::ddays(1),
    post_vax_post_season_days =
      (`Post-season` - `Post-vax`) / lubridate::ddays(1),
  ) %>%
  select(pid, pre_post_vax_days, post_vax_post_season_days)

cdc_obj1_pre_postvax_time %>%
  count(pid) %>%
  filter(n > 1)

cdc_obj1_participants_extra <- cdc_obj1_participants %>%
  inner_join(cdc_obj1_prior_vacs, "pid") %>%
  inner_join(cdc_obj1_first_bleeds, "pid") %>%
  inner_join(cdc_obj1_pre_postvax_time, "pid") %>%
  mutate(age_first_bleed = (first_bleed - dob) / lubridate::dyears(1)) %>%
  select(-first_bleed)

# There are some participants without corresponding HI data
compare_vectors(cdc_obj1_participants_extra$pid, cdc_obj1_hi_extra$pid)

cdc_obj1_participants_extra %>% filter(!complete.cases(.))
cdc_obj1_hi_extra %>% filter(!complete.cases(.))

cdc_obj1_participants_extra %>%
  count(pid) %>%
  filter(n > 1)

save_data(cdc_obj1_participants_extra, "cdc-obj1-participant")
save_data(cdc_obj1_hi_extra, "cdc-obj1-hi")
save_data(cdc_obj1_vax_hist, "cdc-obj1-vax-hist")

# Participants for objective 2 -----------------------------------------

cdc_obj2_participants_raw <- read_raw_csv(
  "cdc-obj2/prior_vacc",
  col_types = cols()
)

cdc_obj2_participants <- cdc_obj2_participants_raw %>%
  mutate(
    dob = if_else(
      is.na(dob),
      lubridate::ymd(paste0(yob, "-07-01")),
      dob %>% parse_access_date("19")
    )
  ) %>%
  select(pid = study_id, gender = Sex, dob = dob, site = Site)

cdc_obj2_participants %>%
  count(pid) %>%
  filter(n > 1)
cdc_obj2_participants %>% filter(!complete.cases(.))

# Prior vaccinations for objective 2 ------------

cdc_obj2_vax_hist <- cdc_obj2_participants_raw %>%
  select(pid = study_id, matches("Vacc_\\d{4}c$")) %>%
  pivot_longer(contains("Vacc"), names_to = "year", values_to = "status") %>%
  mutate(
    year = str_replace(year, "Vacc_(\\d{4})c", "\\1") %>% as.numeric()
  )

cdc_obj2_vax_hist %>% filter(!complete.cases(.))

compare_vectors(cdc_obj2_vax_hist$pid, cdc_obj2_participants$pid)

# HI for objective 2 ------------------------

cdc_obj2_hi_raw <- read_raw_csv("cdc-obj2/HI", col_types = cols())

cdc_obj2_hi <- cdc_obj2_hi_raw %>%
  select(sample_id = Sample_ID, titre = Titer, virus_n = VirusN) %>%
  filter(!is.na(titre))

cdc_obj2_hi %>% filter(!complete.cases(.))

# Dates of bleeds ------------------------------

cdc_obj2_dates_raw <- read_raw_csv("cdc-obj2/Dates", col_types = cols())

cdc_obj2_dates <- cdc_obj2_dates_raw %>%
  mutate(
    study_year = str_replace(Specimen_Year, "Year ", "") %>% as.numeric(),
    timepoint = recode_3_timepoints(Blood_Draw_c),
    bleed_date = parse_access_date(Blood_Draw_date, "20")
  ) %>%
  select(
    pid = study_id, sample_id = Specimen_ID, study_year,
    timepoint, bleed_date
  ) %>%
  # There are rows where everything is empty
  filter(
    !(is.na(pid) & is.na(sample_id) &
      is.na(study_year) & is.na(timepoint) & is.na(bleed_date))
  )

cdc_obj2_dates %>%
  count(sample_id) %>%
  filter(n > 1)

# Everyone should have been first bled in study year 1 (2016)
cdc_obj2_dates %>%
  group_by(pid) %>%
  summarise(date_first_bleed = min(bleed_date), .groups = "drop") %>%
  filter(lubridate::year(date_first_bleed) != 2016)

compare_vectors(cdc_obj2_participants$pid, cdc_obj2_dates$pid)

# @FOLLOWUP
# Specimen id in hi data does not match to dates
cdc_obj2_hi_no_dates <- compare_vectors(
  cdc_obj2_hi$sample_id, cdc_obj2_dates$sample_id, "hi", "dates"
) %>%
  select(hi) %>%
  filter(!is.na(hi))

# Add required columns ----------------------------

# Age at first bleed for participants

cdc_obj2_first_bleeds <- cdc_obj2_dates %>%
  group_by(pid) %>%
  summarise(
    date_first_bleed = min(bleed_date),
    recruitment_year = min(study_year),
    .groups = "drop"
  )

cdc_obj2_first_bleeds %>%
  count(pid) %>%
  filter(n > 1)
compare_vectors(cdc_obj2_first_bleeds$pid, cdc_obj2_participants$pid)

unique(cdc_obj2_vax_hist$year)
unique(cdc_obj2_first_bleeds$recruitment_year)

# Vaccinations in 5 years prior to recrutment
cdc_obj2_prior_vacs <- cdc_obj2_vax_hist %>%
  # We know everyone was recruited in 2016
  # And we know all of this vaccination history is relevant (2011-2015)
  group_by(pid) %>%
  summarise(prior_vacs = sum(status), .groups = "drop")

cdc_obj2_participants_extra <- cdc_obj2_participants %>%
  inner_join(cdc_obj2_first_bleeds, "pid") %>%
  inner_join(cdc_obj2_prior_vacs, "pid") %>%
  mutate(
    age_first_bleed = (date_first_bleed - dob) / lubridate::dyears(1)
  ) %>%
  select(-date_first_bleed)

# HI needs pids, years, timepoints, dates and virus names

cdc_obj2_hi_extra <- cdc_obj2_hi %>%
  inner_join(cdc_viruses, "virus_n") %>%
  # Expect to lose some here while there are sample ids that don't match to
  # the ones in dates
  inner_join(cdc_obj2_dates, "sample_id") %>%
  select(pid, study_year, timepoint, bleed_date, virus_full, titre)

cdc_obj2_hi_extra %>%
  count(pid, timepoint, study_year, virus_full) %>%
  filter(n > 1)

save_data(cdc_obj2_participants_extra, "cdc-obj2-participant")
save_data(cdc_obj2_hi_extra, "cdc-obj2-hi")
save_data(cdc_obj2_vax_hist, "cdc-obj2-vax-hist")

# Objective 3 ----------------------------

# Participants

cdc_obj3_participants_raw <- read_raw_csv(
  name = "cdc-obj3/Obj3_Vacc_History",
  col_types = cols(),
  na = c("NA", "", ".")
)

cdc_obj3_participants <- cdc_obj3_participants_raw %>%
  mutate(
    dob_og = dob,
    dob = if_else(
      is.na(dob),
      lubridate::ymd(paste0(dob_range1, "-07-01")),
      parse_access_date(dob, "19")
    )
  ) %>%
  select(pid = study_id, dob, gender = Sex, site = Site)

# Infections info

cdc_obj3_infection_raw <- read_raw_csv(
  "cdc-obj3/HCP_CDC_Obj3_All",
  col_types = cols()
)

compare_vectors(cdc_obj3_infection_raw$`Study ID`, cdc_obj3_participants$pid)

cdc_obj3_infections <- cdc_obj3_infection_raw %>%
  mutate(
    infected = `Case/Control` == "Case",
    infection_year = `Infection Year`
  ) %>%
  filter(infected) %>%
  select(pid = `Study ID`, infection_year) %>%
  distinct()

save_data(cdc_obj3_infections, "cdc-obj3-infections")

# Vaccination history

cdc_obj3_vax_hist <- cdc_obj3_participants_raw %>%
  select(pid = study_id, contains("V_")) %>%
  pivot_longer(contains("V_"), names_to = "year", values_to = "status") %>%
  mutate(year = str_replace(year, "V_", "") %>% as.integer()) %>%
  filter(!is.na(status))

# Dates

cdc_obj3_dates_raw <- read_raw_csv(
  name = "cdc-obj3/Obj3_Vacc_Dates",
  col_types = cols()
) %>%
  # There is some phantom data reading as NA
  filter(!is.na(study_id))

cdc_obj3_dates_raw2 <- read_raw_csv(
  name = "cdc-obj3/HCP_CDC_Obj3_All",
  col_types = cols()
)

cdc_obj3_dates1 <- cdc_obj3_dates_raw %>%
  mutate(
    timepoint = recode_3_timepoints(Blood_Draw_n),
    date = parse_access_date(Blood_Draw_date, "20")
  ) %>%
  select(
    pid = study_id, study_year = Specimen_Year, timepoint,
    date, sample_id = Specimen_ID
  )

cdc_obj3_dates2 <- cdc_obj3_dates_raw2 %>%
  mutate(timepoint = recode_3_timepoints(`Blood Draw c`)) %>%
  select(
    pid = `Study ID`, sample_id = `Specimen ID`,
    timepoint, study_year = Specimen_Year
  )

cdc_obj3_dates <- cdc_obj3_dates2 %>%
  left_join(
    cdc_obj3_dates1, c("pid", "sample_id", "study_year", "timepoint")
  ) %>%
  group_by(pid, timepoint, study_year) %>%
  mutate(date = na.omit(date)) %>%
  ungroup()

# HI

cdc_obj3_hi_raw <- read_raw_csv("cdc-obj3/HI", col_types = cols())

cdc_obj3_hi <- cdc_obj3_hi_raw %>%
  select(sample_id = Sample_ID, virus_n = VirusN, titre = Titer)

# Join the necessary

# First bleeds
cdc_obj3_first_bleeds <- cdc_obj3_dates %>%
  group_by(pid) %>%
  summarise(
    date_first_bleed = min(date),
    recruitment_year = min(study_year),
    .groups = "drop"
  )

# Vaccinations in the past 5 years
cdc_obj3_prior_vacs <- cdc_obj3_vax_hist %>%
  inner_join(cdc_obj3_first_bleeds, "pid") %>%
  filter(
    year < recruitment_year + 2015, year >= recruitment_year + 2015 - 5
  ) %>%
  group_by(pid) %>%
  summarise(prior_vacs = sum(status), .groups = "drop")

cdc_obj3_participants_extra <- cdc_obj3_participants %>%
  inner_join(cdc_obj3_prior_vacs, "pid") %>%
  inner_join(cdc_obj3_first_bleeds, "pid") %>%
  mutate(
    age_first_bleed = (date_first_bleed - dob) / lubridate::dyears(1),
  ) %>%
  select(-date_first_bleed)

# HI

cdc_obj3_hi_no_dates <- compare_vectors(
  cdc_obj3_hi$sample_id, cdc_obj3_dates$sample_id, "hi", "dates"
) %>%
  select(hi) %>%
  filter(!is.na(hi))

cdc_obj3_hi %>%
  count(sample_id, virus_n) %>%
  filter(n > 1)
cdc_obj3_dates %>%
  count(sample_id) %>%
  filter(n > 1)

cdc_obj3_hi_extra <- cdc_obj3_hi %>%
  inner_join(cdc_viruses, "virus_n") %>%
  # A lot of unexpected loss
  inner_join(cdc_obj3_dates, "sample_id") %>%
  select(pid, study_year, virus_full, timepoint, titre)

cdc_obj3_hi_extra %>%
  count(pid, study_year, virus_full, timepoint) %>%
  filter(n > 1)

cdc_obj3_participants_extra %>%
  count(pid) %>%
  filter(n > 1)

cdc_obj3_vax_hist %>%
  count(pid, year) %>%
  filter(n > 1)

save_data(cdc_obj3_participants_extra, "cdc-obj3-participant")
save_data(cdc_obj3_vax_hist, "cdc-obj3-vax-hist")
save_data(cdc_obj3_hi_extra, "cdc-obj3-hi")

# Objective 4 -----------------------

# Participants
cdc_obj4_participants_raw <- read_raw_csv(
  name = "cdc-obj4/Vax History", col_types = cols(), na = c("NA", "", ".")
)

cdc_obj4_participants <- cdc_obj4_participants_raw %>%
  mutate(
    dob = if_else(
      is.na(dob),
      lubridate::ymd(paste0(dob_range1, "-07-01")),
      parse_access_date(dob, "19")
    )
  ) %>%
  select(pid = study_id, site = Site, gender = Sex, dob)

# Vaccination history
cdc_obj4_vax_hist <- cdc_obj4_participants_raw %>%
  select(pid = study_id, matches("^Vac")) %>%
  pivot_longer(-pid, names_to = "year", values_to = "status") %>%
  mutate(year = str_replace(year, "Vac", "") %>% as.integer()) %>%
  filter(!is.na(status))

# Dates
cdc_obj4_dates_raw <- read_raw_csv("cdc-obj4/Dates", col_types = cols())

cdc_obj4_dates <- cdc_obj4_dates_raw %>%
  mutate(
    timepoint = recode_3_timepoints(Blood_Draw),
    bleed_date = parse_access_date(Blood_Draw_date, "20")
  ) %>%
  select(
    pid = study_id, study_year = Sp_Year, sample_id = Specimen_ID,
    timepoint, bleed_date
  )

# There is one year's worth of data for each participant
cdc_obj4_dates %>%
  count(pid, timepoint) %>%
  filter(n > 1)

# HI
cdc_obj4_hi_raw <- read_raw_csv("cdc-obj4/HI", col_types = cols())

cdc_obj4_hi <- cdc_obj4_hi_raw %>%
  select(sample_id = Sample_ID, virus_n = VirusN, titre = Titer)

# Join up

# First bleeds
cdc_obj4_first_bleeds <- cdc_obj4_dates %>%
  group_by(pid) %>%
  summarise(
    date_first_bleed = min(bleed_date),
    recruitment_year = min(study_year),
    .groups = "drop"
  )

# Prior vaccinations
cdc_obj4_prior_vacs <- cdc_obj4_vax_hist %>%
  inner_join(cdc_obj4_first_bleeds, "pid") %>%
  filter(
    year < recruitment_year + 2015, year >= recruitment_year + 2015 - 5
  ) %>%
  group_by(pid) %>%
  summarise(prior_vacs = sum(status), .groups = "drop")

cdc_obj4_participants_extra <- cdc_obj4_participants %>%
  inner_join(cdc_obj4_prior_vacs, "pid") %>%
  inner_join(cdc_obj4_first_bleeds, "pid") %>%
  mutate(age_first_bleed = (date_first_bleed - dob) / lubridate::dyears(1)) %>%
  select(-date_first_bleed)

cdc_obj4_hi_no_dates <- compare_vectors(
  cdc_obj4_hi$sample_id, cdc_obj4_dates$sample_id, "hi", "dates"
) %>%
  select(hi) %>%
  filter(!is.na(hi))

cdc_obj4_hi_extra <- cdc_obj4_hi %>%
  inner_join(cdc_viruses, "virus_n") %>%
  inner_join(cdc_obj4_dates, "sample_id") %>%
  select(pid, study_year, timepoint, virus_full, titre)

cdc_obj4_hi_extra %>%
  count(pid, timepoint, virus_full) %>%
  filter(n > 1)

save_data(cdc_obj4_participants_extra, "cdc-obj4-participant")
save_data(cdc_obj4_vax_hist, "cdc-obj4-vax-hist")
save_data(cdc_obj4_hi_extra, "cdc-obj4-hi")

# Not matching sample ids -------------------

# Seem to have solved this problem
bind_rows(
  cdc_obj1_hi_no_dates %>% mutate(objective = 1),
  cdc_obj2_hi_no_dates %>% mutate(objective = 2),
  cdc_obj3_hi_no_dates %>% mutate(objective = 3),
  cdc_obj4_hi_no_dates %>% mutate(objective = 4),
) %>%
  rename(sample_id_no_date_match = hi) %>%
  save_csv("cdc-hi-no-date")

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
