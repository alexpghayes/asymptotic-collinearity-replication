library(targets)
library(tarchetypes)
library(crew)

tar_config_set(
  seconds_reporter = 0.5,
  seconds_meta_append = 15,
  reporter_make = "timestamp_positives",
  script = "_targets.R",
  store = "_targets",
  garbage_collection = TRUE
)

tar_option_set(
  packages = c(
    "broom",
    "car",
    "dplyr",
    "fastRG",
    "forcats",
    "glue",
    "ggraph",
    "ggplot2",
    "ggrepel",
    "ggtext",
    "graphlayouts",
    "here",
    "ids",
    "igraph",
    "ivreg",
    "Matrix",
    "paletteer",
    "purrr",
    "rlang",
    "scales",
    "scoringutils",
    "spatialreg",
    "spdep",
    "stringr",
    "tibble",
    "tidyr",
    "tidygraph"
  ),

  # https://books.ropensci.org/targets/performance.html
  controller = crew_controller_local(workers = 8),
  format = "qs",
  memory = "transient",
  garbage_collection = TRUE,
  storage = "worker",
  retrieval = "worker"
)

tar_source()

models <- c(
  "model_lim_tiid",
  "model_lim_tiid2",
  "model_lim_tiid3",
  "model_lim_rdpg_fdelta",
  "model_lim_rdpg_fdelta_drop1"
  # "model_lim_rdpg_gyonly_nodc"
)

static_branch_targets <- tar_map(
  unlist = FALSE,
  values = tibble::tibble(model = rlang::syms(models)),
  tar_target(
    population,
    purrr::map(1:chunk_size, ~ model(n = params$n, k = params$rank)),
    iteration = "list",
    pattern = cross(chunk_indices, map(params))
  ),
  tar_target(
    tbl_graph,
    purrr::map(population, sample_tidygraph),
    iteration = "list",
    pattern = map(population)
  ),
  tar_target(
    estimates,
    purrr::map2(tbl_graph, population, get_estimates),
    iteration = "list",
    pattern = map(tbl_graph, population)
  ),
  tar_target(
    combined_estimates,
    purrr::map_dfr(
      estimates,
      params_helper,
      params = params,
      chunk = chunk_indices,
      .id = "rep_within_chunk"
    ),
    pattern = map(estimates, cross(chunk_indices, map(params)))
  )
)

list(
  tar_target(chunk_size, 20),
  tar_target(num_chunks, 5),
  tar_target(n, c(100, 163, 264, 430, 698, 1135, 1845)),
  tar_target(rank, 4),
  tar_target(chunk_indices, 1:num_chunks),
  tar_target(
    params,
    tibble::tibble(
      n = n,
      rank = rank
    ),
    pattern = cross(n, rank)
  ),
  static_branch_targets,
  tar_combine(
    all_estimates,
    static_branch_targets[[4]],
    command = bind_rows(!!!.x)
  ),
  tar_target(
    biometrika_plots_pdf,
    plot_biometrika_figures(
      all_estimates
    ),
    format = "file"
  )
)
