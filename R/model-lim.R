model_lim <- function(A_model, trt, subclass, ..., alpha = 3, beta = 0.2,
                      gamma = rep(4, ncol(trt)), delta = 2, sigma = 0.1,
                      model_name = NULL) {
  stopifnot(inherits(A_model, "undirected_factor_model"))
  stopifnot(is.matrix(trt)) # should have reasonable column names
  stopifnot(length(alpha) == 1)
  stopifnot(length(beta) == 1)
  stopifnot(abs(beta) < 1)
  stopifnot(ncol(trt) == length(gamma))
  # stopifnot(length(delta) == 1)

  # ### for univariate trt
  #
  # Y = alpha + beta GY + gamma trt + delta G trt + varepsilon
  #
  # ### for multivariate trt
  #
  # Y = alpha + beta GY + gamma trt + delta G trt[, 1] + varepsilon
  #
  # ### generative process
  #
  # Y = (I - beta G)^{-1} (alpha + gamma trt + delta G trt[, 1] + varepsilon)
  #
  # where varepsilon ~ N(0, sigma^2)

  check_dots_unnamed()
  
  if (is.null(model_name)) {
    model_name <- subclass
  }
  
  gamma <- seq_along(gamma) + 0.5

  model <- list(
    A_model = A_model,
    trt = trt,
    alpha = alpha,
    beta = beta,
    gamma = gamma,
    delta = delta,
    sigma = sigma,
    model_name = model_name,
    ...
  )
  
  model$id <- ids::uuid()

  class(model) <- c(subclass, "lim")
  model
}

# sample_tidygraph generic from fastRG
sample_tidygraph.lim <- function(model, ...) {
  
  trt_df <- tibble::as_tibble(model$trt)

  graph <- sample_tidygraph(
    model$A_model
  ) |>
    tidygraph::arrange(as.numeric(name)) |>
    tidygraph::mutate(!!!trt_df)
  
  # outcome before contagion
  
  G <- Matrix::rowScale(igraph::as_adj(graph), 1 / igraph::degree(graph))
  Gtrt1 <- G %*% model$trt[, 1, drop = FALSE]
  varepsilon <- stats::rnorm(model$A_model$n, sd = model$sigma)
  
  undiffused <- model$alpha + model$trt %*% model$gamma + model$delta * Gtrt1 + varepsilon
  
  # diffusion to create contagion
  
  I <- Matrix::Diagonal(model$A_model$n, 1)
  y <- qr.solve(I - model$beta * G, undiffused)
  Gy <- G %*% y
  
  graph |>
    tidygraph::mutate(
      intercept = rep(1, model$A_model$n),
      y = drop(y),
      Gy = drop(Gy),
      Gtrt1 = drop(Gtrt1)
    )
}

# sample_tidygraph generic from fastRG
sample_tidygraph.rdpg_fdelta <- function(model, ...) {
  
  trt_df <- tibble::as_tibble(model$trt)
  
  graph <- sample_tidygraph(
    model$A_model
  ) |>
    tidygraph::arrange(as.numeric(name)) |>
    tidygraph::mutate(!!!trt_df)
  
  # outcome before contagion
  
  G <- Matrix::rowScale(igraph::as_adj(graph), 1 / igraph::degree(graph))
  Gtrt <- G %*% model$trt
  colnames(Gtrt) <- paste0("Gtrt", 1:ncol(Gtrt))
  varepsilon <- stats::rnorm(model$A_model$n, sd = model$sigma)
  
  undiffused <- model$alpha + model$trt %*% model$gamma + Gtrt %*% model$delta + varepsilon
  
  # diffusion to create contagion
  
  I <- Matrix::Diagonal(model$A_model$n, 1)
  y <- qr.solve(I - model$beta * G, undiffused)
  Gy <- G %*% y
  
  Gtrt_df <- Gtrt |> 
    as.matrix() |> 
    as_tibble()
  
  graph |>
    mutate(
      intercept = rep(1, model$A_model$n),
      y = drop(y),
      Gy = drop(Gy)
    ) |> 
    mutate(
      !!!Gtrt_df
    )
}

# sample_tidygraph generic from fastRG
sample_tidygraph.rdpg_fdelta_drop1 <- function(model, ...) {
  
  trt_df <- tibble::as_tibble(model$trt)
  
  graph <- sample_tidygraph(
    model$A_model
  ) |>
    tidygraph::arrange(as.numeric(name)) |>
    tidygraph::mutate(!!!trt_df)
  
  # outcome before contagion
  
  G <- Matrix::rowScale(igraph::as_adj(graph), 1 / igraph::degree(graph))
  Gtrt <- G %*% model$trt[, -c(1:2), drop = FALSE]  # for identifiability
  colnames(Gtrt) <- paste0("Gtrt", 1:ncol(Gtrt))
  varepsilon <- stats::rnorm(model$A_model$n, sd = model$sigma)
  
  undiffused <- model$alpha + model$trt %*% model$gamma + Gtrt %*% model$delta + varepsilon
  
  # diffusion to create contagion
  
  I <- Matrix::Diagonal(model$A_model$n, 1)
  y <- qr.solve(I - model$beta * G, undiffused)
  Gy <- G %*% y
  
  Gtrt_df <- Gtrt |> 
    as.matrix() |> 
    as_tibble()
  
  graph |>
    mutate(
      intercept = rep(1, model$A_model$n),
      y = drop(y),
      Gy = drop(Gy)
    ) |> 
    mutate(
      !!!Gtrt_df
    )
}


# sample_tidygraph generic from fastRG
sample_tidygraph.rdpg_fdelta_only <- function(model, ...) {

  graph <- sample_tidygraph(
    model$A_model
  ) |>
    tidygraph::arrange(as.numeric(name)) 
  
  # outcome before contagion
  
  G <- Matrix::rowScale(igraph::as_adj(graph), 1 / igraph::degree(graph))
  Gtrt <- G %*% model$trt
  colnames(Gtrt) <- paste0("Gtrt", 1:ncol(Gtrt))
  varepsilon <- stats::rnorm(model$A_model$n, sd = model$sigma)
  
  y <- model$alpha + Gtrt %*% model$delta + varepsilon
  
  Gtrt_df <- Gtrt |> 
    as.matrix() |> 
    as_tibble()
  
  graph |>
    mutate(
      intercept = rep(1, model$A_model$n),
      y = drop(y)
    ) |> 
    mutate(
      !!!Gtrt_df
    )
}

# sample_tidygraph generic from fastRG
sample_tidygraph.rdpg_gyonly <- function(model, ...) {
  
  trt_df <- tibble::as_tibble(model$trt)
  
  graph <- sample_tidygraph(
    model$A_model
  ) |>
    tidygraph::arrange(as.numeric(name)) |>
    tidygraph::mutate(!!!trt_df)
  
  # outcome before contagion
  
  G <- Matrix::rowScale(igraph::as_adj(graph), 1 / igraph::degree(graph))
  varepsilon <- stats::rnorm(model$A_model$n, sd = model$sigma)
  
  undiffused <- model$alpha + model$trt %*% model$gamma + varepsilon
  
  # diffusion to create contagion
  
  I <- Matrix::Diagonal(model$A_model$n, 1)
  y <- qr.solve(I - model$beta * G, undiffused)
  Gy <- G %*% y
  
  graph |>
    tidygraph::mutate(
      intercept = rep(1, model$A_model$n),
      y = drop(y),
      Gy = drop(Gy)
    )
}
