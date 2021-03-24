read_data <- function(name) {
  common <- list()
  cdc_timepoint <- col_factor(c("prevax", "postvax", "postseas"))
  specific <- list(
    "hi-hanam" = list(
      timepoint = col_factor(c("BL", "d4", "d7", "d14", "d21", "d280"))
    ),
    "cdc-obj1-hi" = list(
      timepoint = cdc_timepoint
    ),
    "cdc-hi-obj2" = list(
      timepoint = cdc_timepoint
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
