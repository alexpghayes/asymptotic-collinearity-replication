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
model_lim_rdpg_gyonly_nodc_rankdef <- function(n, k = NULL) {
  
  diagonal <- 0.5
  off_diagonal <- 0.05
  
  # NOTE: k argument is ignored!!
  k <- 4
  
  # B from Keith Warren paper
  U <- matrix(c(0.7, 0.1, 0.2, 0.5, 0.2, 0.6, 0.2, 0.5), 4, 2)
  B <- U %*% t(U)
  
  
  pi <- rep(1, k) / k

  A_model <- dcsbm(
    theta = rep(1, n),
    B = B,
    pi = pi,
    allow_self_loops = FALSE,
    poisson_edges = TRUE,
    expected_degree = 2 * n^0.6
  )
  
  # fastRG doesn't play nice with rank-deficient models
  A_model$k <- 2
  
  trt <- ASE(A_model)
  
  # note: intercept can be in or close to the space of X, so drop first column
  # of X just in case
  colnames(trt) <- paste0("trt", left_padded_sequence(1:ncol(trt)))
  
  model_lim(
    A_model = A_model,
    trt = trt,
    subclass = "rdpg_gyonly",
    model_name = "rdpg_gyonly_nodc_rankdef"
  )
}
