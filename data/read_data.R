read_data <- function(name) {
  common <- list()
  cdc_timepoint <- col_factor(c("prevax", "postvax", "postseas"))
  cdc_timepoint2 <- col_factor(c("Pre-vax", "Post-vax", "Post-season"))

  specific <- list(
    "hi-hanam" = list(
      timepoint = col_factor(c("BL", "d4", "d7", "d14", "d21", "d280"))
    ),
    "cdc-obj1-hi" = list(
      timepoint = cdc_timepoint2
    ),
    "cdc-obj2-hi" = list(
      timepoint = cdc_timepoint2
    ),
    "cdc-obj3-hi" = list(
      timepoint = cdc_timepoint2
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
