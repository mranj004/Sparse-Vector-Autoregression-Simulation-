# 03_generate_VAR.R
# Simulate VAR(1) using a fixed A and generated innovations

generate_VAR <- function(
    A,
    errors,
    burn_in = 200,
    y0 = NULL
) {
  
  # Check dimensions
  d <- nrow(A)
  
  
  # Verify stability
  spec_rad <- max(Mod(eigen(A, only.values = TRUE)$values))

  # Initialize starting observation
  if (is.null(y0)) {
    y0 <- rep(0, d)
  }
  
  # Number of innovations
  n <- nrow(errors)
  
  # Initialize time series
  Y <- matrix(0, nrow = n + 1, ncol = d)
  
  Y[1, ] <- y0
  
  # Generate observations recursively
  for (t in seq_len(n)) {
    
    Y[t + 1, ] <- as.vector(
      A %*% Y[t, ]
    ) + errors[t, ]
    
  }
  
  # Remove burn-in observations
  Y <- Y[(burn_in + 2):(n + 1), , drop = FALSE]
  
  return(Y)
}
