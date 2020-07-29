cat("extract data for analysis")

library(tidyverse)

data_dir <- "data"
data_raw_dir <- "data-raw"

# Functions ===================================================================

read_raw <- function(name, ...) {
  readxl::read_excel(file.path(data_raw_dir, paste0(name, ".xlsx")), ...)
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

# HI results

hi <- hi_raw %>%
  select(id = Sample_ID, virus = Virus, titre = Titer)

save_csv(hi, "hi")
