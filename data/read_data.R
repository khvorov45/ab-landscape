read_data <- function(name) {
  if (name == "hi-hanam") {
    columns <- cols(
      titre = col_integer(),
      timepoint = col_factor(c("BL", "d4", "d7", "d14", "d21", "d280")),
      timepoint_num = col_integer(),
      virus_year = col_integer(),
      clade = col_character()
    )
  } else {
    columns <- cols(
      titre = col_integer(),
      timepoint = col_factor(c("Pre-vax", "Post-vax", "Post-season")),
      timepoint_num = col_integer(),
      virus_year = col_integer(),
      clade = col_character()
    )
  }
  readr::read_csv(
    file.path("data", paste0(name, ".csv")),
    col_types = columns
  )
}
