plot_biometrika_figures <- function(all_estimates) {
  if (!dir.exists(here("figures/simulations/"))) {
    dir.create(here("figures/simulations/"))
  }

  estimates_labelled <- all_estimates |>
    filter(model != "rdpg_gyonly_nodc") |>
    mutate(
      i = str_extract(term, "[0-9]"),
      label = case_when(
        term == "intercept" ~ "widehat(alpha)",
        term == "Gy" ~ "widehat(beta)",
        str_detect(term, "Gtrt") ~ paste0("widehat(", str_replace(term, "Gtrt", "delta)["), "]"),
        str_detect(term, "^trt") ~ paste0("widehat(", str_replace(term, "trt", "gamma)["), "]")
      ),
      label_nice = case_when(
        term == "intercept" ~ "widehat(alpha)",
        term == "Gy" ~ "widehat(beta)",
        str_detect(term, "Gtrt") & model %in% c("tiid", "tiid2", "tiid3") ~ "widehat(delta)",
        str_detect(term, "^trt") & model %in% c("tiid", "tiid2", "tiid3") ~ "widehat(gamma)",
        str_detect(term, "Gtrt") ~ paste0("widehat(", str_replace(term, "Gtrt", "delta)["), "]"),
        str_detect(term, "^trt") ~ paste0("widehat(", str_replace(term, "trt", "gamma)["), "]")
      ),
      asymptotically_identified = case_when(
        model %in% c("tiid", "tiid2", "tiid3") & str_detect(label, "alpha|beta|delta") ~ FALSE,
        model %in% c("tiid", "tiid2", "tiid3") & str_detect(label, "gamma") ~ TRUE,
        model == "rdpg_fdelta" & str_detect(label, "alpha|beta|delta") ~ FALSE,
        model == "rdpg_fdelta" & str_detect(label, "gamma") ~ TRUE,
        model == "rdpg_gyonly_nodc" ~ FALSE,
        model == "rdpg_fdelta_drop1" ~ TRUE
      ),
      color_hack = if_else(!asymptotically_identified, "firebrick", NA), # label_nice, NA),
      model_name = case_match(
        model,
        "rdpg_fdelta" ~ "Unrestricted",
        "rdpg_fdelta_drop1" ~ "Restricted",
        "rdpg_gyonly_nodc" ~ "SBM: Contagion only",
        "tiid" ~ "Bernoulli",
        "tiid2" ~ "Bernoulli 0.5",
        "tiid3" ~ "Bernoulli 0.3"
      ),
      estimator_name = case_match(
        estimator,
        "ols" ~ "Ordinary Least Squares",
        "tsls" ~ "Two-Stage Least Squares",
        "qmle" ~ "Quasi-Maximum Likelihood"
      ),
      estimator_name = fct_relevel(
        estimator_name,
        "Ordinary Least Squares",
        "Two-Stage Least Squares",
        "Quasi-Maximum Likelihood"
      ),
      model_name = as.factor(model_name),
      model_fct = fct_relevel(
        model_name,
        "Bernoulli",
        "Bernoulli 0.5",
        "Bernoulli 0.3",
        "Unrestricted",
        "Restricted"
      ),
      interval_score = interval_score(
        true_values = true,
        lower = conf.low,
        upper = conf.high,
        interval_range = 80,
        weigh = TRUE,
        separate_results = FALSE
      )
    )

  summarized_df <- estimates_labelled |>
    summarize(
      mean_squared_error = mean(squared_error, na.rm = TRUE),
      sd_squared_error = sd(squared_error, na.rm = TRUE),
      mean_standard_error = mean(std.error, na.rm = TRUE),
      mean_interval_score = mean(interval_score, na.rm = TRUE),
      mean_cond_num = mean(cond_num, na.rm = TRUE),
      coverage = mean(covered, na.rm = TRUE),
      mean_interval_length = mean(abs(conf.high - conf.low), na.rm = TRUE),
      mean_relative_error = mean(abs(bias) / true, na.rm = TRUE),
      mean_vif = mean(vif, na.rm = TRUE),
      sd_vif = sd(vif, na.rm = TRUE),
      perc_na = mean(is.na(estimate)),
      .by = c(label, label_nice, n, model, model_name, model_fct, rank, asymptotically_identified, color_hack, estimator, estimator_name)
    ) |>
    mutate(
      mse_lo = mean_squared_error - 2 * sd_squared_error,
      mse_hi = mean_squared_error + 2 * sd_squared_error,
      vif_lo = mean_vif - 2 * sd_vif,
      vif_hi = mean_vif + 2 * sd_vif
    )

  plot1 <- summarized_df |>
    filter(model_name == "Bernoulli") |>
    ggplot(
      aes(
        x = n,
        # ymin = pmin(mse_lo, 0),
        y = mean_squared_error,
        # ymax = mse_hi,
        color = asymptotically_identified,
        group = label,
        label = label_nice,
        linetype = asymptotically_identified
      )
    ) +
    geom_label_repel(
      parse = TRUE,
      data = filter(summarized_df, n == min(n), model_name == "Bernoulli"),
      segment.linetype = "dotted",
      segment.alpha = 0.8,
      nudge_x = -0.25,
      show.legend = FALSE
    ) +
    # geom_point() +
    # geom_ribbon() +
    # geom_errorbar() +
    geom_line() +
    scale_y_log10(labels = label_log(digits = 2)) +
    scale_x_log10(labels = label_log(digits = 2)) +
    scale_color_manual(
      values = c("#cb4154", "grey50"),
      guide = "none"
    ) +
    # scale_color_viridis_d(labels = label_parse(), guide = "none", na.value = "grey", end = 0.85) +
    facet_wrap(
      vars(estimator_name)
    ) +
    labs(
      x = "Number of nodes (log scale)",
      y = "Mean squared error\n(log scale)"
    ) +
    theme_minimal(
      base_size = 11,
      base_family = "TeXGyre Termes"
    ) +
    theme(
      legend.position = "none"
    )

  plot1

  path1 <- here(
    glue("figures/simulations/biometrika-mse.pdf")
  )

  ggsave(
    path1,
    plot = plot1,
    width = 146,
    height = 64,
    units = "mm",
    device = cairo_pdf,
    # https://developer.r-project.org/Blog/public/2020/04/17/changes-to-symbol-fonts-for-cairo-graphics-devices/
    symbolfamily = cairoSymbolFont("TeXGyre Termes")
  )

  plot2 <- summarized_df |>
    filter(estimator == "ols", !(model %in% c("tiid2", "tiid3"))) |>
    ggplot(
      aes(
        x = n,
        # ymin = vif_lo,
        y = mean_vif,
        # ymax = vif_hi,
        color = asymptotically_identified,
        group = label,
        linetype = asymptotically_identified
        # label = label_nice
      )
    ) +
    # geom_label_repel(
    #   parse = TRUE,
    #   aes(label = label_nice),
    #   data = filter(summarized_df, n == min(n), model == "Experiment", estimator == "ols"),
    #   segment.linetype = "dotted",
    #   segment.alpha = 0.8,
    #   nudge_x = -0.2,
    #   show.legend = FALSE
    # ) +
    # geom_point() +
    geom_line() +
    scale_x_log10(labels = label_log(digits = 2)) +
    scale_y_log10(labels = label_log(digits = 2)) +
    # scale_color_viridis_d(labels = label_parse(), guide = "none", na.value = "grey", end = 0.85) +
    scale_color_manual(
      values = c("#cb4154", "grey50"),
      guide = "none"
    ) +
    facet_grid(
      cols = vars(model_fct)
    ) +
    labs(
      x = "Number of nodes (log scale)",
      y = "Average variance inflation factor\n(log scale)",
      parse = TRUE,
      color = "Asymptotically"
    ) +
    theme_minimal(
      base_size = 11,
      base_family = "TeXGyre Termes"
    ) +
    theme(
      legend.position = "none"
    )

  path2 <- here("figures", "simulations", "biometrika-vif.pdf")

  ggsave(
    path2,
    plot = plot2,
    width = 146,
    height = 64,
    units = "mm",
    device = cairo_pdf,
    symbolfamily = cairoSymbolFont("TeXGyre Termes")
  )

  plot3 <- summarized_df |>
    filter(!(model %in% c("tiid2", "tiid3"))) |>
    mutate(Model = model_fct, Estimator = toupper(estimator)) |>
    ggplot(
      aes(
        x = n,
        # ymin = mse_lo,
        y = mean_squared_error,
        # ymax = mse_hi,
        color = asymptotically_identified,
        group = label,
        label = label,
        linetype = asymptotically_identified
      )
    ) +
    # geom_label_repel(
    #   parse = TRUE,
    #   data = filter(summarized_df, n == min(n)),
    #   segment.linetype = "dotted",
    #   segment.alpha = 0.8,
    #   nudge_x = -0.75
    # ) +
    # geom_point() +
    geom_line() +
    scale_y_log10(labels = label_log(digits = 2)) +
    scale_x_log10(labels = label_log(digits = 2)) +
    # scale_color_viridis_d(labels = label_parse(), guide = "none", na.value = "grey", end = 0.85) +
    scale_color_manual(
      values = c("#cb4154", "grey50"),
      guide = "none"
    ) +
    facet_grid(
      rows = vars(Model),
      cols = vars(estimator_name)
    ) +
    labs(
      x = "Number of nodes (log scale)",
      y = "Mean squared error (log scale)",
      color = "Estimate"
    ) +
    theme_minimal(
      base_size = 11,
      base_family = "TeXGyre Termes"
    ) +
    theme(
      legend.position = "none"
    )

  plot3

  path3 <- here(
    glue("figures/simulations/biometrika-mse-all.pdf")
  )

  ggsave(
    path3,
    plot = plot3,
    width = 146,
    height = 95,
    units = "mm",
    device = cairo_pdf,
    symbolfamily = cairoSymbolFont("TeXGyre Termes")
  )

  c(path1, path2, path3)
}

plot_single_replication <- function() {
  # we simulate from very dense graphs to mimic the asymptotic regime
  # fastRG::plot_sparse_matrix(A)

  library(tidyverse)
  library(tidygraph)
  library(igraph)
  library(paletteer)
  library(ggraph)
  library(glue)
  library(here)

  tbl_graph <- tar_read(tbl_graph_model_lim_tiid_7245db95)

  gs <- tbl_graph[[1]] |>
    igraph::simplify() |>
    as_tbl_graph()

  set.seed(27)

  bb <- layout_as_backbone(gs, keep = 0.2)

  gs <- gs |>
    activate(nodes) |>
    mutate(
      grp = as.factor(group_fluid(4) %% 2 == 0),
      condition = as.factor(ifelse(trt1, "treated", "control"))
    )

  E(gs)$backbone <- FALSE
  E(gs)$backbone[bb$backbone] <- TRUE

  plot1 <- ggraph(gs, layout = "manual", x = bb$xy[, 1], y = bb$xy[, 2]) +
    geom_edge_link0(aes(color = backbone), width = 0.2) +
    geom_node_point(aes(color = condition)) +
    scale_color_paletteer_d("ggsci::category10_d3", guide = "none") +
    labs(
      # title = "Randomized experiment on a stochastic blockmodel",
      caption = glue("Treatments are assigned by coin flip and {100 - round(transitivity(gs), 2) * 100}% of triangles are open")
    ) +
    scale_edge_color_manual(
      values = c(rgb(0, 0, 0, 0.1), rgb(0, 0, 0, 0.45)),
      guide = "none"
    ) +
    theme_void(
      base_size = 10.5,
      base_family = "TeXGyre Termes"
    )

  path1 <- here(
    glue("figures/simulations/biometrika-backbone.pdf")
  )

  width <- 74

  ggsave(
    path1,
    plot = plot1,
    width = 146,
    height = 72,
    units = "mm",
    device = cairo_pdf,
    symbolfamily = cairoSymbolFont("TeXGyre Termes"),
    dpi = 600
  )


  tbl_graph2 <- tar_read(tbl_graph_model_lim_rdpg_fdelta_drop1_90de3dd8)

  gs2 <- tbl_graph2[[1]] |>
    simplify() |>
    as_tbl_graph()

  bb2 <- layout_as_backbone(gs2, keep = 0.2)

  gs2 <- gs2 |>
    activate(nodes) |>
    mutate(
      condition = trt2
    )

  E(gs2)$backbone <- FALSE
  E(gs2)$backbone[bb2$backbone] <- TRUE

  plot2 <- ggraph(gs2, layout = "manual", x = bb$xy[, 1], y = bb$xy[, 2]) +
    geom_edge_link0(aes(color = backbone), width = 0.2) +
    geom_node_point(aes(color = condition)) +
    # scale_color_gradient2(high = "blue", low = "red") +
    scale_color_distiller(palette = "RdYlBu") +
    labs(
      # subtitle = "Treatments dependent on block membership",
      caption = "SBM has four blocks and mild degree correction",
      color = "X"
    ) +
    scale_edge_color_manual(
      values = c(rgb(0, 0, 0, 0.1), rgb(0, 0, 0, 0.45)),
      guide = "none"
    ) +
    theme_void(
      base_size = 10.5,
      base_family = "TeXGyre Termes"
    )

  path2 <- here(
    glue("figures/simulations/biometrika-backbone-dependent.pdf")
  )

  width2 <- 82

  ggsave(
    path2,
    plot = plot2,
    width = width2 * 16 / 9,
    height = width2,
    units = "mm",
    device = cairo_pdf,
    symbolfamily = cairoSymbolFont("TeXGyre Termes"),
    dpi = 600
  )

  c(path1, path2)
}
