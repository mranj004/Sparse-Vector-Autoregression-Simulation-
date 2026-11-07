# Generate a sparse, stable coefficient matrix for VAR(1)

generate_sparse_A <- function(
  d = 10,
  spar = 0.8, #this controls off-diagonal sparsity
  rho_A = 0.8,
  seed = NULL #completely optional
){
  
  if (!is.null(seed)) {
      set.seed(seed)
    }
  A <- matrix(0, nrow = d, ncol = d)
  
  diag(A) <- runif(d, min= 0.2, max=0.6) #generating diag via Uniform
  
  offDiag <- which(row(A) != col(A), arr.ind = TRUE)
  
  nnz <- round((1-spar) * nrow(offDiag))
  
  if (nnz > 0){ #Giving nonzero coeff
    selected <- sample.int(nrow(offDiag), nnz)
    
    signs <- sample(c(-1,1), nnz, replace = TRUE)
    
    val <- runif(nnz, min=0.1, max=0.5)
    
    A[offDiag[selected, ,drop = FALSE]] <- val * signs
  }
  
  spec_rad <- max(Mod(eigen(A, only.values =TRUE)$values))
  
  A <- A * (rho_A / spec_rad)
  
  return(A)
}

