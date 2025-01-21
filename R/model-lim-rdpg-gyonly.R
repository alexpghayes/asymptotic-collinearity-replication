#' @examples
#'
#' set.seed(26)
#'
#' b <- model_lim_rdpg(n = 1000, k = 5)
#'
#' graph <- sample_tidygraph(b)
#' graph
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
model_lim_rdpg_gyonly <- function(n, k = 5) {
  
  diagonal <- 0.5
  off_diagonal <- 0.05
  
  B <- matrix(off_diagonal, nrow = k, ncol = k)
  diag(B) <- diagonal
  
  pi <- rep(1, k) / k
  
  theta <- runif(n, 1, 5) # rexp(n)^1.5 + 1
  
  A_model <- dcsbm(
    theta = theta,  # heterogeneous degrees
    B = B,
    pi = pi,
    allow_self_loops = FALSE,
    poisson_edges = TRUE,
    expected_degree = 2 * n^0.6
  )
  
  trt <- ASE(A_model)
  
  # note: intercept can be in or close to the space of X, so drop first column
  # of X just in case
  colnames(trt) <- paste0("trt", left_padded_sequence(1:ncol(trt)))
  
  model_lim(
    A_model = A_model,
    trt = trt,
    subclass = "rdpg_gyonly"
  )
}
