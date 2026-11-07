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
  
  errorTypes <- c(   # Innovation covariance structures
    "block_I",
    "block_II",
    "toeplitz"
  )
  
  A <- generate_sparse_A(  # Generate true A once
    d = d,
    spar = spar,
    rho_A = rho_A,
    seed = seedA
  )
  
  set.seed(seedSim)  # Set random seed for simulation
  
  resultsList <- list()  # Store simulation results
  k <- 1
  
  for (rep in seq_len(n_rep)) {  # Repeat experiments
    # Generate standard Gaussian random variables
    # Reuse Z across covariance structures
    Z <- matrix(
      rnorm((T + burn_in) * d),
      nrow = T + burn_in,
      ncol = d
    )
    
    for (type in errorTypes) {
      for (rho in rhoVal) {
        
        #Generate innovations
        E <- generate_errors(
          n = T + burn_in,
          d = d,
          type = type,
          rho = rho,
          Z = Z
        )
        
        # Generate VAR(1) observations
        Y <- generate_VAR(
          A = A,
          errors = E,
          burn_in = burn_in
        )
        
        #Estimate A using baseline lambda=0.05
        A_hatFixed <- estLASSO(
          Y = Y,
          lambda = lambda
        )
        
        #Adaptive Lambda sparsevar
        fitSparse <- estSparseVAR(Y)
        A_hatSparse <- fitSparse$A_hat
        lambdaSparse <- fitSparse$lambda
        
        #Evaluate estimation accuracy
        metricsFixed <- eval_A(A, A_hatFixed)
        metricsSparse <- eval_A(A, A_hatSparse)
        
        #Storing results
        resultsList[[k]] <-  rbind(
          data.frame(
            rep = rep,
            errorType = type,
            rho = rho,
            T = T,
            d = d,
            method = "glmnet_fixed",
            lambda = lambda,
            metricsFixed
          ),
          data.frame(
            rep = rep,
            errorType = type,
            rho = rho,
            T = T,
            d = d,
            method = "sparsevar_adaptive",
            lambda = lambdaSparse,
            metricsSparse
          )
        )
        k <- k + 1
      }
    }
  }
  
  results <- do.call(rbind, resultsList)  # Combine all simulation results
  rownames(results) <- NULL
  
  safeMean <- function(x) {   # Mean excluding missing values
    x <- x[is.finite(x)]
    
    if (length(x) == 0) {
      return(NA_real_)
    }
    return(mean(x))
  }
  
  safeSD <- function(x) {   #Standard deviation excluding missing values
    x <- x[is.finite(x)]
    
    if (length(x) < 2) {
      return(NA_real_)
    }
    
    return(sd(x))
  }
  
  groups <- split(
    results,
    interaction(
      results$method,
      results$errorType,
      results$rho,
      drop = TRUE
    )
  )
  
  summaryList <- lapply(groups, function(x) {
    data.frame(
      method = x$method[1],
      errorType = x$errorType[1],
      rho = x$rho[1],
      meanLambda = safeMean(x$lambda),
      meanRelError = safeMean(x$relError),
      sdRelError = safeSD(x$relError),
      meanMSE = safeMean(x$mse),
      meanPrecision = safeMean(x$precision),
      meanRecall = safeMean(x$recall),
      meanF1 = safeMean(x$f1)
    )
  })
  
  summary <- do.call(rbind, summaryList)
  summary <- summary[
    order(
      match(summary$errorType, errorTypes),
      summary$rho,
      summary$method
    ),
    ,
    drop = FALSE
  ]
  rownames(summary) <- NULL
  
  # Construct paired comparisons
  fixed <- subset(
    results,
    method == "glmnet_fixed",
    select = c(rep, errorType, rho, relError, mse)
  )
  
  adaptive <- subset(
    results,
    method == "sparsevar_adaptive",
    select = c(rep, errorType, rho, relError, mse)
  )
  
  paired <- merge(
    fixed,
    adaptive,
    by = c("rep", "errorType", "rho"),
    suffixes = c("Fixed", "Adaptive")
  )
  
  # Negative delta means adaptive sparsevar performed better
  paired$deltaRelError <- (
    paired$relErrorAdaptive - paired$relErrorFixed
  )
  
  paired$deltaMSE <- (
    paired$mseAdaptive - paired$mseFixed
  )
  
  if (saveResults) { #saving results
    
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
      paired,
      "results/paired_comparison.csv",
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
    summary = summary,
    paired = paired
  ))
}