read_data <- function(name) {
  common <- list()
  specific <- list(
    "hi-hanam" = list(
      timepoint = col_factor(c("BL", "d4", "d7", "d14", "d21", "d280"))
    ),
    "cdc-hi-obj1" = list(
      timepoint = col_factor(c("prevax", "postvax", "postseas"))
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
