maxPvalue <-
function(Mat, Ns, int, beta, delta, doublePvalue){
  Mat[is.na(Mat)] <- 0
  index <- 1
  prob <- rep(NA, length(int))
  
  # Instead of calculating all binomials, just calculate the binomials that are needed
  Tbls <- which(Mat==1, arr.ind = TRUE) - 1
  maxX1 <- max(Tbls[ , 1])  #Always have upper triangle, so minX1 is always 0
  minX2 <- min(Tbls[ , 2])  #Always have upper triangle, so maxX2 is always Ns[2]
  
  for (probVal in int) {
    #prob[index] <- sum(dbinom(0:Ns[1], Ns[1], probVal)*(Mat %*% dbinom(0:Ns[2], Ns[2], probVal)))
    prob[index] <- suppressWarnings(sum(dbinom(0:maxX1, Ns[1], probVal+delta)*(Mat[1:(maxX1+1), (minX2+1):(Ns[2]+1), drop=FALSE] %*% dbinom(minX2:Ns[2], Ns[2], probVal))))
    index <- index + 1
  }
  prob <- signif((1+doublePvalue)*prob + beta, 12) #Remove rounding errors and double probabilities if applicable
  
  # Cutoff probabilities at 1 (only relevant if beta != 0) #
  prob[!is.na(prob) & prob > 1] <- 1
  
  np <- int[which(prob==max(prob, na.rm=TRUE))]
  pvalue <- max(prob, na.rm=TRUE)
    
  return(list(prob=prob, pvalue=pvalue, np=np))
}
