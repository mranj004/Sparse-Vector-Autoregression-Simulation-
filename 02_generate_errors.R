
# 02_generate_errors.R
# Gaussian VAR(1) innovations following
# Basu & Michailidis (2015), Section 6.2.
#
# Covariance structures:
# 1. Block-I
# 2. Block-II
# 3. Toeplitz
#
# rho = 0.5, 0.7, or 0.9 in the paper.
# epsilon_t ~ N_d(0, Sigma_e), independent over time.


err_cov <- function(
    d,
    type = c("block_I", "block_II", "toeplitz"),
    rho = 0.5
) {
  
  type <- match.arg(type)
  
  stopifnot(
    d >= 2,
    d %% 2 == 0,
    rho >= 0,
    rho < 1
  )
  
  Sigma <- diag(d)
  
  first <- seq_len(d / 2)
  second <- (d / 2 + 1):d
  
  if (type == "block_I") {
    
    # First half of variables are correlated
    Sigma[first, first] <- rho
    diag(Sigma) <- 1
    
  } else if (type == "block_II") {
    
    # Correlation within both blocks
    # No correlation between blocks
    Sigma[first, first] <- rho
    Sigma[second, second] <- rho
    diag(Sigma) <- 1
    
  } else if (type == "toeplitz") {
    
    # Correlation decreases with |i-j|
    ind <- seq_len(d)
    
    Sigma <- outer(
      ind, ind,
      function(i, j) rho^abs(i - j)
    )
  }
  
  return(Sigma)
}


generate_errors <- function(
    n,
    d,
    type = c("block_I", "block_II", "toeplitz"),
    rho = 0.5,
    seed = NULL,
    Z = NULL
) {
  
  type <- match.arg(type)
  
  Sigma <- err_cov(d, type, rho)
  
  if (is.null(Z)) {
    
    if (!is.null(seed)) {
      set.seed(seed)
    }
    
    Z <- matrix(
      rnorm(n * d),
      nrow = n,
      ncol = d
    )
    
  } else {
    
    stopifnot(
      is.matrix(Z),
      is.numeric(Z),
      all(dim(Z) == c(n, d)),
      all(is.finite(Z))
    )
  }
  
  # Transform standard normal observations
  # so each row has covariance Sigma
  errors <- Z %*% chol(Sigma)
  
  return(errors)
}
