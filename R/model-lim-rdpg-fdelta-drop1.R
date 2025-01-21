#' @examples
#'
#' set.seed(26)
#'
#' population <- model_lim_rdpg_fdelta_drop1(n = 1000, k = 4)
#'
#' tbl_graph <- sample_tidygraph(population)
#' graph
#' 
#' fit <- graph |> 
#'   as_tibble() |> 
#'   select(-name) |> 
#'   lm(y ~ ., data = _)
#'   
#' alias(fit, partial = TRUE)
#' 
#' decomp <- graph |>
#'   as_tibble() |> 
#'   select(-name) |> 
#'   skimr::skim()
#'   as.matrix() |> 
#'   qr()
#'
#' graph |>
#'   as_tibble() |> 
#'   select(-name) |> 
#'   
#'   cov() |> 
#'   ggcorrplot::ggcorrplot()
#'   
#' library(profvis)
#' 
#' profvis(sample_tidygraph(b))
#'
#' hist(igraph::degree(graph))
#'
#'
model_lim_rdpg_fdelta_drop1 <- function(n, k = 5) {
  
  diagonal <- 0.5
  off_diagonal <- 0.1
  
  B <- matrix(off_diagonal, nrow = k, ncol = k)
  diag(B) <- diagonal
  
  pi <- rep(1, k) / k
  
  theta <- runif(n, 1, 2) # rexp(n)^1.5 + 1
  
  A_model <- dcsbm(
    theta = theta,  # heterogeneous degrees
    B = B,
    pi = pi,
    allow_self_loops = TRUE,
    poisson_edges = TRUE,
    expected_degree = 2 * n^0.7
  )
  
  X <- ASE(A_model)
  
  # note: intercept can be in or close to the space of X, so drop first column
  # of X just in case
  trt <- X
  colnames(trt) <- paste0("trt", left_padded_sequence(1:ncol(trt)))
  
  model_lim(
    A_model = A_model,
    trt = trt,
    subclass = "rdpg_fdelta_drop1",
    delta = rep(2, ncol(trt) - 2)
  )
}
