library(tidyverse)

# Functions ===================================================================

source("data/read_data.R")

fit_mixak <- function(data,
                      outcome, fixed_formula, n_clusters,
                      burn, keep, thin) {
  y <- data %>%
    select(!!!rlang::syms(outcome)) %>%
    as.data.frame()
  # Subset the complete model matrix because we need a separate one for each
  # outcome
  x_all <- model.matrix(fixed_formula, data)[, -1]
  outcome_names <- colnames(y)
  names(outcome_names) <- colnames(y)
  x <- map(
    outcome_names,
    function(outcome_name) {
      outcome_formula <- fixed_formula %>%
        paste(collapse = "") %>%
        paste0("`", outcome_name, "`", .) %>%
        as.formula()
      relevant_cols <- model.matrix(outcome_formula, data)[, -1] %>% colnames()
      x_all[, relevant_cols]
    }
  )

  # This produces warnings because someone doesn't know how to work the
  # %in% operator
  fit <- mixAK::GLMM_MCMC(
    # Outcome
    y = y,
    dist = rep("gaussian", ncol(y)),

    # Random effect clusters
    id = data %>% pull(pid),

    # Fixed effects
    x = x,

    # Random effects
    z = rep(list("empty"), ncol(y)),

    # Random intercept
    random.intercept = rep(TRUE, ncol(y)),
    prior.b = list(
      # Clusters
      Kmax = n_clusters
    ),

    # MCMC settings
    nMCMC = c(burn = burn, keep = keep, thin = thin, info = 100),
    parallel = TRUE
  ) %>%
    mixAK::NMixRelabel(type = "stephens", keep.comp.prob = TRUE)
  fit$outcome_names <- colnames(y)
  fit$fixed_names <- map(x, colnames)
  fit$ids <- tibble(pid = unique(data$pid)) %>% mutate(index = row_number())
  fit$n_clusters <- n_clusters
  fit$mcmc_settings <- tibble(burn = burn, keep = keep, thin = thin)
  fit
}

# Get the random effects from the horrible fit object
pull_random <- function(fit) {
  map_dfr(
    1:2,
    ~ fit[[.x]]$mixture_b %>%
      as_tibble() %>%
      select(-contains("Corr")) %>%
      rename_with(
        function(x) {
          paste(
            fit$outcome_names, "intercept",
            tolower(x) %>%
              str_replace("^b\\.", "") %>%
              str_replace("\\.\\d$", ""),
            sep = "_"
          )
        }
      ) %>%
      mutate(chain = .x, iteration = row_number())
  )
}

# Get the fixed effects from the horrible fit object
pull_fixed <- function(fit) {
  map_dfr(
    1:2,
    ~ fit[[.x]]$alpha %>%
      as_tibble() %>%
      rename_with(
        function(x) {
          imap(fit$fixed_names, ~ paste(.y, .x, sep = "_")) %>%
            flatten() %>%
            as.character()
        }
      ) %>%
      mutate(chain = .x, iteration = row_number())
  )
}

lengthen_estimates <- function(data) {
  data %>%
    pivot_longer(
      c(-chain, -iteration),
      names_to = c("virus", "parameter"),
      values_to = "estimate",
      names_pattern = "([^_]*)_(.*)"
    )
}

tidy_mixak <- function(fit) {
  bind_rows(
    pull_fixed(fit) %>% lengthen_estimates() %>% mutate(effect = "fixed"),
    pull_random(fit) %>% lengthen_estimates() %>% mutate(effect = "random")
  )
}

cluster_assignment_mixak <- function(fit) {
  map_dfr(1:2, function(c) {
    fit[[c]]$comp.prob %>%
      as_tibble() %>%
      mutate(iteration = row_number()) %>%
      pivot_longer(contains("P"), names_to = "of", values_to = "prob") %>%
      mutate(
        index = str_replace(of, "^P\\((\\d+),\\d+\\)$", "\\1") %>%
          as.integer(),
        cluster = str_replace(of, "^P\\(\\d+,(\\d+)\\)$", "\\1") %>%
          as.integer(),
        chain = c,
      ) %>%
      inner_join(fit$ids, by = "index") %>%
      select(-of, -index)
  })
}

pull_diagnostics <- function(fit) {
  tibble(
    outcome = list(fit$outcome_names),
    fixed = list(fit$fixed_names),
    n_clusters = fit$n_clusters,
    ped = fit$PED[["PED"]],
  ) %>%
    bind_cols(fit$mcmc_settings)
}

save_data <- function(data, name) {
  write_csv(
    data %>%
      mutate(across(where(is.list), ~ map_chr(., jsonlite::toJSON))),
    glue::glue("cluster/{name}.csv")
  )
  data
}

# Script ======================================================================

cdc_obj1_hi <- read_data("cdc-obj1-hi")
cdc_obj1_partcipant <- read_data("cdc-obj1-participant")

cdc_obj1_titres_wide <- cdc_obj1_hi %>%
  inner_join(cdc_obj1_partcipant, "pid") %>%
  mutate(
    logtitremid = if_else(titre == 5, log(5), log(titre) + log(2) / 2),
    timepoint_n = as.integer(timepoint) %>% as.factor(),
    prior_vacs_3cat = case_when(
      prior_vacs == 0 ~ "0",
      prior_vacs == 5 ~ "5",
      TRUE ~ "1_4"
    )
  ) %>%
  select(-titre, -timepoint, -prior_vacs) %>%
  pivot_wider(names_from = "virus_full", values_from = "logtitremid")

# Fit with different clusters
future::plan(future::multisession)
diff_cl <- furrr::future_map(
  1:2,
  ~ fit_mixak(
    cdc_obj1_titres_wide,
    c(
      "a/uruguay/716/07e", "a/newcastle/30/16", "a/victoria/511/04",
      "a/thailand/409/05", "a/singapore/16-0019/16e", "a/new caledonia/104/14",
      "a/perth/16/09", "a/texas/50/12e", "a/perth/16/09e",
      "a/switzerland/9715293/13e", "a/hong kong/4801/14e", "a/philippines/2/82",
      "a/netherlands/620/89", "a/netherlands/179/93", "a/netherlands/178/95",
      "a/tasmania/1/97", "a/kansas/14/17", "a/sydney/22/18", "a/kansas/14/17e",
      "a/switzerland/8060/17e", "a/brisbane/117/18", "a/victoria/361/11e",
      "a/townsville/2/99", "a/brisbane/60/18", "a/south australia/34/19",
      "a/bilthoven/16190/68", "a/bilthoven/21793/72", "a/philippines/472/02",
      "a/brisbane/10/07", "a/bilthoven/1761/76", "a/hanoi/eli444/10",
      "a/victoria/361/11", "a/texas/50/12", "a/switzerland/9715293/13",
      "a/bilthoven/2271/76"
    ),
    ~ timepoint_n + prior_vacs_3cat + age_first_bleed,
    n_clusters = .x,
    burn = 100,
    keep = 500,
    thin = 1
  )
)

diff_cl %>%
  map_dfr(pull_diagnostics) %>%
  save_data("diag")

map_with_cl <- function(.x, .f, ...) {
  map_dfr(.x, ~ .f(.x) %>% mutate(n_clusters = .x$n_clusters), ...)
}
diff_cl %>%
  map_with_cl(tidy_mixak) %>%
  save_data("parameters")

diff_cl %>%
  map_with_cl(cluster_assignment_mixak) %>%
  group_by(cluster, chain, pid, n_clusters) %>%
  summarise(
    prob_mn = mean(prob),
    prob_q025 = quantile(prob, 0.025),
    prob_q25 = quantile(prob, 0.25),
    prob_q50 = quantile(prob, 0.5),
    prob_q75 = quantile(prob, 0.75),
    prob_q975 = quantile(prob, 0.975),
    .groups = "drop"
  ) %>%
  save_data("clusters")
