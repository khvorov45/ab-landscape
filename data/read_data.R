read_data <- function(name) {
  readr::read_csv(
    file.path("data", paste0(name, ".csv")),
    col_types = cols(
      titre = col_integer(),
      timepoint = col_integer(),
      virus_year = col_integer()
    )
  )
}
