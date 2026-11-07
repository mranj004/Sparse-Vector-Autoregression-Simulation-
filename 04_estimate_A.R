
# Estimate A in a sparse VAR(1) using row-wise LASSO

estLASSO <- function(Y, lambda = 0.05) {
  
  T <- nrow(Y)
  d <- ncol(Y)
  
  # Lagged predictors
  X <- Y[1:(T-1), , drop = FALSE]
  
  # Response variables
  repY <- Y[2:T, , drop = FALSE]
  
  A_hat <- matrix(0, nrow = d, ncol = d)
  
  lambda_glm <- lambda * T / (T-1)
  
  for (i in seq_len(d)) {
    
    fit <- glmnet::glmnet(
      x = X,
      y = repY[, i],
      family = "gaussian",
      alpha = 1,
      lambda = lambda_glm * c(2, 1, 0.5),
      intercept = FALSE,
      standardize = FALSE
    )
    
    A_hat[i, ] <- as.numeric(
      coef(fit, s = lambda_glm)
    )[-1]
  }
  
  return(A_hat)
}
