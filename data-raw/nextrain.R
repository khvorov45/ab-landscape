library(tidyverse)

nextstrain_tree <- httr::GET(
  "https://nextstrain.org/charon/getDataset?prefix=/flu/seasonal/h3n2/ha/12y"
) %>%
  httr::content()

process_child <- function(child) {
  children <- tibble()
  if (!is.null(child$children)) {
    children <- map_dfr(child$children, process_child)
  }
  bind_rows(
    tibble(name = child$name, clade = child$node_attrs$clade_membership$value),
    children
  )
}

nextstrain_viruses <- map_dfr(nextstrain_tree$tree$children, process_child) %>%
  filter(!str_starts(name, "NODE"))

nextstain_freqs <- httr::GET(
  "https://nextstrain.org/charon/getDataset?prefix=/flu/seasonal/h3n2/ha/12y&type=tip-frequencies"
) %>%
  httr::content()

process_freq <- function(freq, name, pivots) {
  if (name == "generated_by" | name == "pivots") {
    return(tibble())
  }
  imap_dfr(
    freq$frequencies,
    ~ tibble(name = name, n = .y, freq = .x, year = pivots[[.y]])
  )
}

nextstrain_freq_table <- imap_dfr(
  nextstain_freqs, process_freq, nextstain_freqs$pivots
)

setdiff(nextstrain_freq_table$name, nextstrain_viruses$name)
setdiff(nextstrain_viruses$name, nextstrain_freq_table$name)

freq_table_extra <- nextstrain_freq_table %>%
  inner_join(nextstrain_viruses, "name")

clade_frequencies <- freq_table_extra %>%
  group_by(year, clade) %>%
  summarise(freq = sum(freq), .groups = "drop")

write_csv(clade_frequencies, "data-raw/nexstrain-clade-frequencies.csv")
