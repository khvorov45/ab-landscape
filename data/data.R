cat("extract data for analysis")

library(tidyverse)

data_dir <- "data"
data_raw_dir <- "data-raw"

# Functions ===================================================================

read_raw <- function(name, ...) {
  readxl::read_excel(file.path(data_raw_dir, paste0(name, ".xlsx")), ...)
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

# Script ======================================================================

# The raw data

viruses_raw <- read_raw("Viruses")
dilutions_raw <- read_raw("VirusDiln_BT")
sera_raw <- read_raw("Sera List")
samples_raw <- read_raw("Obj1_Sample_info")
hi_raw <- read_raw("HI")

# Select the required columns

sera <- select(sera_raw, sample = Sample_ID, pid = PID, timepoint = Timepoint)
hi <- select(
  hi_raw,
  sample = Sample_ID, virus = Virus, virus_n = VirusN, titre = Titer
)
viruses <- select(
  viruses_raw,
  virus = Virus_Name, virus_n = VirusN, virus_year = Year
) %>%
  mutate(virus_year = as.integer(virus_year))
participants <- select(
  samples_raw,
  pid = PID, group = `Case/Control`, sex = Sex, age = Age
) %>%
  distinct(pid, .keep_all = TRUE)

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

# Trust viruses for names
hi_no_virname <- select(hi, -virus)

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
