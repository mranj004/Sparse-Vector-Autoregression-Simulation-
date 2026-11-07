# Estimate A in a sparse VAR(1) using row-wise LASSO using 
#GLMNET and SparseVAR


#Using a fixed lambda via glmnet
estLASSO <- function(Y, lambda = 0.05) {
  T <- nrow(Y)
  d <- ncol(Y)
  
  X <- Y[1:(T-1), , drop = FALSE]  # Lagged predictors
  repY <- Y[2:T, , drop = FALSE]  # Response variables
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


# Adaptive LASSO using sparsevar
estSparseVAR <- function(Y) {
  fit <- sparsevar::fit_var(
    data = Y,
    p = 1,
    penalty = "ENET",
    method = "timeSlice",
    alpha = 1,
    center = FALSE, #avoids centering simulated zero-mean VAR
    scale = FALSE,
    threshold = FALSE,
    leaveOut = 10
  )
  A_hat <- as.matrix(fit$A[[1]])
  
  return(list(
    A_hat = A_hat,
    lambda = fit$lambda
  ))
}

