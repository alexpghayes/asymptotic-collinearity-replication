# implementation taken from Keith Warren therapeutic community pape

### Utility functions

### Compute Derivatives and Hessians of the objective function for bias corrected and uncorrected model

funcder <- function(rho, YL, Y, lambdas, M) {
  n <- length(Y)
  s2 <- t(Y) %*% M %*% Y - 2 * rho * t(Y) %*% M %*% YL + rho^2 * t(YL) %*% M %*% YL
  derval <- (2 / n) * sum(lambdas / (1 - rho * lambdas)) + 2 * (rho * t(YL) %*% M %*% YL - t(Y) %*% M %*% YL) / s2
  return(derval)
}

funchess <- function(rho, YL, Y, lambdas, M) {
  n <- length(Y)
  s2 <- t(Y) %*% M %*% Y - 2 * rho * t(Y) %*% M %*% YL + rho^2 * t(YL) %*% M %*% YL
  hessval <- (2 / n) * sum((lambdas / (1 - rho * lambdas))^2) + 2 * t(YL) %*% M %*% YL / s2 - 4 * (rho * t(YL) %*% M %*% YL - t(Y) %*% M %*% YL)^2 / s2^2
  return(hessval)
}

# A simple matrix trace function
tr <- function(X) {
  return(sum(diag(X)))
}


### No latent variable, only predictor Z

## Y is univariate response
## X is binary adjacency matrix (undirected)
## L is the Laplacian or the row-normalized adjacency (can be weighted) matrix
## d is the number of dimensions - can either be given or automatically estimated
## Z is the set of additional predictors.


# might not support sparse eigen decomp
sarMLclassicZ <- function(Y, L, X, Z) {
  n <- length(Y)
  d2 <- dim(as.matrix(Z))[2]
  lambdas <- eigen(L)$values
  YL <- L %*% Y
  nhat <- Z
  M <- diag(n) - nhat %*% solve(t(nhat) %*% (nhat)) %*% t(nhat)
  rholast <- (t(Y) %*% M %*% YL) / (t(Y) %*% M %*% Y)
  rhodiff <- 1
  t <- 0
  
  print(summary(lambdas))
  
  ### optimize for rho
  
  # NOTE: this optimization step seems to fail sometimes, and it seems to do so
  # based on bad initializations. instead of totally fixing this rho search
  # process in this code i didn't write, i'm going to cheat and initialize
  # rho at the true value so that this grad descent procedure has a chance,
  # which seems to work well in practice (i.e. lead to estimates of rho similar
  # to those estimated by different methods)
  rholast <- 0.2 # as.numeric(rholast)
  # print(glue("rholast: {rholast},  1 / lambda_max : {1 / max(lambdas)}, 1 / lambda_min: {1 / min(lambdas)}"))
  while (rhodiff > 0.001 & t < 50) {
    rho <- rholast - funcder(rholast, YL, Y, lambdas, M) / funchess(rholast, YL, Y, lambdas, M)
    rho <- as.numeric(rho)
    rhodiff <- abs(rho - rholast)
    rholast <- rho
    t <- t + 1
    # print(glue("rholast: {rholast},  1 / lambda_max : {1 / max(lambdas)}, 1 / lambda_min: {1 / min(lambdas)}"))
  }
  
  # try to handle type hell -- coerce 1x1 Matrix to length-1 numeric vector
  rho <- as.numeric(rho)
  
  ### compute estimates for beta and sigmasq
  hatZ <- (diag(n) - rho * L) %*% Y
  beta <- solve(t(nhat) %*% nhat) %*% t(nhat) %*% hatZ
  sigmasqhat <- (t(Y) %*% M %*% Y - 2 * rho * t(Y) %*% M %*% YL + rho^2 * t(YL) %*% M %*% YL) / n
  
  G <- L %*% solve(diag(n) - rho * L)
  H <- G %*% (nhat %*% beta)
  
  ## compute the Fisher information matrix
  
  Infmatrix <- matrix(NA, d2 + 2, d2 + 2)
  Infmatrix[1:d2, ] <- cbind(t(nhat) %*% nhat, t(nhat) %*% H, rep(0, d2))
  Infmatrix[d2 + 1, ] <- cbind(t(H) %*% nhat, t(H) %*% H + sigmasqhat * tr((G + t(G)) %*% G), tr(G))
  Infmatrix[d2 + 2, ] <- cbind(t(rep(0, d2)), tr(G), n / (2 * sigmasqhat))
  
  Infmatrix <- Infmatrix / c(n * sigmasqhat)
  
  ## compute the standard errors
  
  cov <- solve(Infmatrix)
  serho <- sqrt(diag(cov)[d2 + 1] / n)
  sebeta <- sqrt(diag(cov)[1:d2] / n)
  
  prho <- 2 * (1 - pnorm(abs(rho) / serho))
  pbeta <- 2 * (1 - pnorm(abs(beta) / sebeta))
  
  ### Return as list the estimates, SEs, and p-values of influence (rho) and coefficeint of Z (beta), along with estimate of variance of error term
  
  return(list(influence = rho, SEinfluence = serho, pvaluerho = prho, betacovariate = beta, SEbetacovariate = sebeta, pvaluecovariate = pbeta, sigmasqhat = sigmasqhat))
}
