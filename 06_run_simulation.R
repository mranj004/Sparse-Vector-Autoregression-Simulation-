
# 06_run_simulations.R
# Sparse VAR(1) coefficient recovery under
# different Gaussian innovation covariance structures

# Load functions
source("01_generate_A.R")
source("02_generate_errors.R")
source("03_generate_VAR.R")
source("04_estimate_A.R")
source("05_evaluate_A.R")


runSim <- function(
    T = 300,
    d = 10,
    spar = 0.8,
    rho_A = 0.8,
    lambda = 0.05,
    burn_in = 200,
    n_rep = 30,
    rhoVal = c(0.5, 0.7, 0.9),
    seedA = 293,
    seedSim = 123,
    saveResults = TRUE
) {
  
  stopifnot(
    T >= 3,
    d >= 2,
    d %% 2 == 0,
    n_rep >= 1,
    burn_in >= 0,
    lambda > 0,
    length(rhoVal) > 0,
    all(is.finite(rhoVal)),
    all(rhoVal >= 0 & rhoVal < 1)
  )
  
  # Innovation covariance structures
  errorTypes <- c(
    "block_I",
    "block_II",
    "toeplitz"
  )
  
  # Generate true A once
  A <- generate_sparse_A(
    d = d,
    spar = spar,
    rho = rho_A,
    seed = seedA
  )
  
  # Set random seed for simulation
  set.seed(seedSim)
  
  # Store simulation results
  resultsList <- list()
  k <- 1
  
  # Repeat experiments
  for (rep in seq_len(n_rep)) {
    
    # Generate standard Gaussian random variables
    # Reuse Z across covariance structures
    Z <- matrix(
      rnorm((T + burn_in) * d),
      nrow = T + burn_in,
      ncol = d
    )
    
    for (type in errorTypes) {
      
      for (rho in rhoVal) {
        
        # Step 1: Generate innovations
        E <- generate_errors(
          n = T + burn_in,
          d = d,
          type = type,
          rho = rho,
          Z = Z
        )
        
        # Step 2: Generate VAR(1) observations
        Y <- generate_VAR(
          A = A,
          errors = E,
          burn_in = burn_in
        )
        
        # Step 3: Estimate A using LASSO
        A_hat <- estLASSO(
          Y = Y,
          lambda = lambda
        )
        
        # Step 4: Evaluate estimation accuracy
        metrics <- eval_A(A, A_hat)
        
        # Step 5: Store results
        resultsList[[k]] <- data.frame(
          rep = rep,
          errorType = type,
          rho = rho,
          T = T,
          d = d,
          lambda = lambda,
          metrics
        )
        
        k <- k + 1
      }
    }
  }
  
  # Combine all simulation results
  results <- do.call(rbind, resultsList)
  rownames(results) <- NULL
  
  # Mean excluding missing values
  safeMean <- function(x) {
    x <- x[is.finite(x)]
    
    if (length(x) == 0) {
      return(NA_real_)
    }
    
    return(mean(x))
  }
  
  # Standard deviation excluding missing values
  safeSD <- function(x) {
    x <- x[is.finite(x)]
    
    if (length(x) < 2) {
      return(NA_real_)
    }
    
    return(sd(x))
  }
  
  # Group results by error structure and rho
  groups <- split(
    results,
    interaction(
      results$errorType,
      results$rho,
      drop = TRUE
    )
  )
  
  # Calculate summary statistics
  summaryList <- lapply(groups, function(x) {
    
    data.frame(
      errorType = x$errorType[1],
      rho = x$rho[1],
      meanRelError = safeMean(x$relError),
      sdRelError = safeSD(x$relError),
      meanMSE = safeMean(x$mse),
      meanPrecision = safeMean(x$precision),
      meanRecall = safeMean(x$recall),
      meanF1 = safeMean(x$f1)
    )
  })
  
  summary <- do.call(rbind, summaryList)
  
  # Sort summary by error structure and rho
  summary <- summary[
    order(
      match(summary$errorType, errorTypes),
      summary$rho
    ),
    ,
    drop = FALSE
  ]
  
  rownames(summary) <- NULL
  
  # Save results
  if (saveResults) {
    
    dir.create(
      "results",
      showWarnings = FALSE,
      recursive = TRUE
    )
    
    write.csv(
      results,
      "results/simulation_results.csv",
      row.names = FALSE
    )
    
    write.csv(
      summary,
      "results/simulation_summary.csv",
      row.names = FALSE
    )
    
    write.csv(
      A,
      "results/true_A.csv",
      row.names = FALSE
    )
  }
  
  return(list(
    A = A,
    results = results,
    summary = summary
  ))
}
