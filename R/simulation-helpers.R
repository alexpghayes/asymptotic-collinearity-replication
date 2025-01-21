left_padded_sequence <- function(x) {
  
  original <- withr::with_options(
    c(scipen = 999),
    as.character(x)
  )
  
  max_digits <- max(vapply(original, nchar, integer(1)))
  formatC(x, width = max_digits, format = "d", flag = "0")
}

ASE <- function(ufm) {
  s <- fastRG::svds(ufm)
  k <- ufm$k
  
  if (k > 1) {
    S <- diag(sqrt(s$d * 2))
  } else {
    S <- matrix(sqrt(s$d * 2))
  }
  
  US <- s$u %*% S
  colnames(US) <- paste0("US", 1:ncol(US))
  US
}

US <- function(A, rank, ...) {
  s <- RSpectra::svds(A, k = rank, nu = rank, nv = rank, ...)
  if (rank > 1) {
    us <- s$u %*% diag(sqrt(s$d))
  } else {
    us <- s$u %*% matrix(sqrt(s$d))
  }
  
  colnames(us) <- as.character(1:rank)
  us
}

get_estimates <- function(tbl_graph, population) {
  
  true_df <- tibble(
    term = c(
      "intercept",
      paste0("trt", 1:length(population$gamma)),
      paste0("Gtrt", 1:length(population$delta)),
      "Gy"
    ),
    true = c(population$alpha, population$gamma, population$delta, population$beta)
  )
  
  A <- igraph::as_adj(tbl_graph)
  G <- Matrix::rowScale(A, 1 / rowSums(A))
  
  node_data <- as_tibble(tbl_graph)
  
  ols_fit <- node_data |> 
      dplyr::select(matches("^trt"), intercept, matches("^Gtrt"), Gy, y) |> 
      lm(y ~ . + 0, data = _)
  
  # we do some slighly hackier stuff for the 2SLS because partial
  # identifiability in the restricted model means we can't include all
  # of Gtrt. but including all of G G trt is fine because having too
  # many instruments shouldn't be a problem (even if co-linear??). tbd
  
  node_cols <- colnames(node_data)
  Gtrt_names <- paste0(node_cols[grep("^Gtrt", node_cols)], collapse = " + ")
  
  GGtrt <- G %*% (G %*% population$trt)
  
  colnames(GGtrt) <- paste0("GGtrt", 1:ncol(GGtrt))
  GGtrt_df <- as.data.frame(as.matrix(GGtrt))
  
  trts <- colnames(population$trt)
  
  trt_names <- paste0(trts, collapse = " + ")
  GGtrt_names <- paste0("GG", trts, collapse = " + ")
  
  formula <- as.formula(
    glue(
      "y ~ intercept + Gy + {trt_names} + {Gtrt_names} + 0 | {trt_names} + {Gtrt_names} + {GGtrt_names}"
    )
  )
  
  tsls_fit <- node_data |> 
    bind_cols(
      !!!GGtrt_df
    ) |> 
    ivreg(
      formula,
      data = _
    )
  
  # listw <- mat2listw(A, style = "W", zero.policy = TRUE)
  
  # mle_fit <- node_data |> 
  #   dplyr::select(matches("^trt"), y) |> 
  #   lagsarlm(y ~ ., data = _, listw = listw, Durbin = TRUE)
  # 
  # tidied_mle <- mle_fit |> 
  #   tidy(conf.int = TRUE, conf.level = 0.8) |> 
  #   mutate(
  #     estimator = "mle",
  #     term = str_replace_all(term, "lag\\.", "G"),
  #     term = str_replace_all(term, "rho", "Gy"),
  #     term = str_replace_all(term, "\\(Intercept\\)", "intercept")
  #   )
  
  Z <- select(node_data, intercept, matches("^trt"), matches("^Gtrt")) |> 
    as.matrix()
  
  Y <- node_data$y
  
  qmle_fit <- sarMLclassicZ(
    Y = Y,
    L = as.matrix(G),
    X = as.matrix(A),
    Z = Z
  )
  
  alpha <- 0.2
  crit_val <- qnorm(1 - alpha / 2)
  
  tidied_qmle <- tibble(
    term = c(rownames(qmle_fit$betacovariate), "Gy"),
    estimate = c(qmle_fit$betacovariate, qmle_fit$influence),
    std.error = c(qmle_fit$SEbetacovariate, qmle_fit$SEinfluence)
  ) |> 
    mutate(
      conf.low = estimate - crit_val * std.error,
      conf.high = estimate + crit_val * std.error,
      estimator = "qmle"
    )
  
  tidied_tsls <- tsls_fit |> 
    tidy(conf.int = TRUE, conf.level = 0.8) |> 
    mutate(
      estimator = "tsls"
    )
    
  tidied_ols <- ols_fit |> 
    tidy(conf.int = TRUE, conf.level = 0.8) |> 
    mutate(
      estimator = "ols",
    )
  
  tidied <- bind_rows(
    tidied_qmle,
    tidied_tsls,
    tidied_ols
  ) |> 
    full_join(true_df, by = "term") |> 
    mutate(
      bias = estimate - true,
      squared_error = (estimate - true)^2,
      model = population$model_name,
      covered = conf.low <= true & true <= conf.high,
      cond_num = kappa(ols_fit)
    )
  
  # IGNORE: In vif.default(fit) : No intercept: vifs may not be sensible.
  vifs <- try(suppressWarnings(vif(ols_fit)))
  
  if (inherits(vifs, "try-error")) {
    # catch colinearity error (distinct from warning suppressed above!)
    tidied_plus <- tidied
    tidied_plus$vif <- -1
  } else {
    vif_df <- enframe(vifs, name = "term", value = "vif")
    tidied_plus <- left_join(tidied, vif_df, by = join_by(term))
  }
  
  tidied_plus
}

params_helper <- function(tbl, params, chunk) {
  tbl$n <- params$n
  tbl$rank <- params$rank
  tbl$chunk <- chunk
  tbl
}
