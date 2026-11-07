#Evaluating how well $\hat{A}$ estimates A
eval_A <- function(
    A,
    A_hat,
    tol =1e-6
){
  relError <- norm(A_hat - A, "F")/ norm(A, "F") #using Frobenious estimation from Section 6.2
  
  mse <- mean((A_hat - A)^2) #MSE
  
  #Finding the non-zero coefficients
  trueNZ <- abs(A) > tol
  estNZ <- abs(A_hat) > tol
  
  # confusion matrix countrs
  TP <- sum(trueNZ & estNZ)
  FP <- sum(!trueNZ & estNZ)
  FN <- sum(trueNZ & !estNZ)
  TN <- sum(!trueNZ & !estNZ)
  
  precision <- if ((TP + FP) > 0) { #calculating precision
    TP / (TP + FP)
  } else{
    NA_real_
  }
  
  recall <- if ((TP + FN) > 0){
    TP / (TP + FN)
  } else{
    NA_real_
  }
  
  f1 <- if (!is.na(precision) && 
            !is.na(recall) &&
            (precision + recall) > 0){
    2 * precision * recall / (precision + recall)
  } else{
    NA_real_
  }
  return(data.frame( # Return evaluation metrics
    relError = relError,
    mse = mse,
    precision = precision,
    recall = recall,
    f1 = f1,
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN
  ))
}