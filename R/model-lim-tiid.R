#' @examples
#'
#' set.seed(26)
#'
#' b <- model_lim_tiid(n = 1000, k = 5)
#'
#' graph <- sample_tidygraph(b)
#' graph
#' 
#' node_data <- as_tibble(graph)
#'   
#' ols_fit <- lm(y ~ trt1 + Gy + Gtrt1, data = node_data)
#' summary(ols_fit)
#'   
#'
#' # fit the model with two stage least squares
#' # NOTE: requires dev version of sphet
#' library(sphet)  
#' 
#' A <- igraph::as_adj(graph)
#' G <- Matrix::rowScale(A, 1 / rowSums(A))
#' 
#' #' 
#' node_data$trt2 <- node_data$trt1 * 2 + rnorm(1000)
#' 
#' library(spatialreg)
#' 
#' 
#' tsls_fit <- spreg(
#'   y ~ trt1,
#'   data = node_data,
#'   listw = G,
#'   model = "lag",
#'   Durbin = TRUE
#' )
#' 
#' summary(tsls_fit)
#' 
#' 
#' decomp <- graph |>
#'   as_tibble() |> 
#'   select(-name) |> 
#'   as.matrix() |> 
#'   qr()
#'
#' graph |>
#'   as_tibble() |> 
#'   select(-name) |> 
#'   skimr::skim()
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
model_lim_tiid <- function(n, k = 5) {
  
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
  
  trt <- matrix(
    rbinom(n, size = 1, prob = 0.5),
    ncol = 1
  )
  
  colnames(trt) <- "trt1"
  
  model_lim(
    A_model = A_model,
    trt = trt,
    subclass = "tiid"
  )
}

model_lim_tiid2 <- function(n, k = 5) {
  
  diagonal <- 0.5
  off_diagonal <- 0.1
  
  B <- matrix(off_diagonal, nrow = k, ncol = k)
  diag(B) <- diagonal
  
  pi <- rep(1, k) / k
  
  theta <- runif(n, 1, 2) #rexp(n)^1.5 + 1
  
  A_model <- dcsbm(
    theta = theta,  # heterogeneous degrees
    B = B,
    pi = pi,
    allow_self_loops = TRUE,
    poisson_edges = TRUE,
    expected_degree = 2 * n^0.5
  )
  
  trt <- matrix(
    rbinom(n, size = 1, prob = 0.5),
    ncol = 1
  )
  
  colnames(trt) <- "trt1"
  
  model_lim(
    A_model = A_model,
    trt = trt,
    subclass = "tiid2"
  )
}

model_lim_tiid3 <- function(n, k = 5) {
  
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
    expected_degree = 2 * log(n)
  )
  
  trt <- matrix(
    rbinom(n, size = 1, prob = 0.5),
    ncol = 1
  )
  
  colnames(trt) <- "trt1"
  
  model_lim(
    A_model = A_model,
    trt = trt,
    subclass = "tiid3"
  )
}

