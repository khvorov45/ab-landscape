read_data <- function(name) {
  common <- list(
    titre = col_integer(),
    timepoint_num = col_integer(),
    virus_year = col_integer(),
    clade = col_character()
  )
  specific <- list(
    "hi-hanam" = list(
      timepoint = col_factor(c("BL", "d4", "d7", "d14", "d21", "d280"))
    ),
    "hi-obj2" = list(
      timepoint = col_factor(c("Pre-vax", "Post-vax", "Post-season"))
    ),
    "hi" = list(
      timepoint = col_factor(c("Pre-vax", "Post-vax", "Post-season"))
    ),
    "hi-rmh-hcw" = list(
      timepoint = col_factor(c("Pre-vax", "Post-vax", "Post-season")),
      group = col_factor(c("Frequent", "Moderate", "Infrequent"))
    )
  )
  readr::read_csv(
    file.path("data", paste0(name, ".csv")),
    col_types = do.call(cols, c(common, specific[[name]]))
  )
}
