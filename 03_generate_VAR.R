# 03_generate_VAR.R
# Simulate VAR(1) using a fixed A and generated innovations

generate_VAR <- function(
    A,
    errors,
    burn_in = 200,
    y0 = NULL
) {
  
  d <- nrow(A)  # Check dimensions

  spec_rad <- max(Mod(eigen(A, only.values = TRUE)$values))  # Verify stability

  if (is.null(y0)) {  # Initialize starting observation
    y0 <- rep(0, d)
  }
  
  n <- nrow(errors)  # Number of innovations
  Y <- matrix(0, nrow = n + 1, ncol = d)  # Initialize time series
  Y[1, ] <- y0
  
  for (t in seq_len(n)) {   # Generate observations recursively
    Y[t + 1, ] <- as.vector(
      A %*% Y[t, ]
    ) + errors[t, ]
  }
  
  Y <- Y[(burn_in + 2):(n + 1), , drop = FALSE]  # Remove burn-in observations
  return(Y)
}
