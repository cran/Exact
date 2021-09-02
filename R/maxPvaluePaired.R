maxPvaluePaired <-
function(Mat, N, int, beta, delta){
  
  index <- 1
  prob <- rep(NA, length(int))
  Tbls <- which(Mat==1, arr.ind = TRUE) - 1
  for (probVal in int) {
    prob[index] <- sum(trinom(Tbls[ , 1], Tbls[ , 2], N, p12=probVal, p21=probVal, delta=delta))
    index <- index + 1
  }
  prob <- signif(prob + beta, 12) #Remove rounding errors
  
  # Cutoff probabilities at 1 (only relevant if beta != 0) #
  prob[!is.na(prob) & prob > 1] <- 1
  
  np <- int[which(prob==max(prob, na.rm=TRUE))]
  pvalue <- max(prob, na.rm=TRUE)
  
  return(list(prob=prob, pvalue=pvalue, np=np))
}
