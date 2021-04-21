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

sim_n_sample <- 1e2
sim_n_timepoints <- 3
random_effect_sd <- 0.3
residual_sd <- 0.6
timepoint_2_effect <- 3
timepoint_3_effect <- 2
b0 <- 3
sim_data <- tibble(
  pid = 1:sim_n_sample,
  cluster = if_else(pid %% 2 == 0, "1", "2"),
  random_intercept = if_else(
    cluster == "1",
    rnorm(sim_n_sample, -1.5, random_effect_sd),
    rnorm(sim_n_sample, 1.5, random_effect_sd)
  ),
) %>%
  slice(rep(1:sim_n_sample, each = sim_n_timepoints)) %>%
  mutate(
    timepoint = rep(1:sim_n_timepoints, times = sim_n_sample),
    logtitre_expected = b0 + random_intercept + case_when(
      timepoint == 2 ~ timepoint_2_effect,
      timepoint == 3 ~ timepoint_3_effect,
      TRUE ~ 0
    ),
    timepoint_factor = as.factor(timepoint),
    logtitre = rnorm(n(), logtitre_expected, residual_sd),
    logtitre2 = rnorm(n(), logtitre_expected, residual_sd),
    logtitre3 = rnorm(n(), logtitre_expected, residual_sd),
  )

sim_data %>%
  ggplot(aes(timepoint, exp(logtitre), col = as.factor(pid))) +
  ggdark::dark_theme_bw(verbose = FALSE) +
  theme(legend.position = "none") +
  scale_y_log10(breaks = 5 * 2^(0:20)) +
  geom_line(alpha = 0.5)

future::plan(future::multisession)
mixak_result <- furrr::future_map(
  1:3, ~ fit_mixak(
    data = sim_data,
    outcome = c("logtitre", "logtitre2"),
    fixed_formula = ~timepoint_factor,
    n_clusters = .x,
    burn = 500, keep = 2000, thin = 1
  ),
  .options = furrr::furrr_options(seed = TRUE)
)

map_dfr(mixak_result, pull_diagnostics)

map_with_cl <- function(.x, .f, ...) {
  map_dfr(.x, ~ .f(.x) %>% mutate(n_clusters = .x$n_clusters), ...)
}
mixak_clusters <- mixak_result %>%
  map_with_cl(cluster_assignment_mixak) %>%
  group_by(cluster, chain, pid, n_clusters) %>%
  summarise(
    prob_q50 = quantile(prob, 0.5),
    .groups = "drop"
  ) %>%
  filter(cluster == 1, n_clusters == 2, chain == 1) %>%
  mutate(cluster_assigned = if_else(prob_q50 > 0.5, "1", "2"))

mixak_clusters %>%
  select(pid, cluster_mixak = cluster_assigned) %>%
  inner_join(
    sim_data %>% select(pid, cluster_true = cluster) %>% distinct(),
    "pid"
  ) %>%
  summarise(
    right = sum(cluster_mixak == cluster_true) / n()
  )
