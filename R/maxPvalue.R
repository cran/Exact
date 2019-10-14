maxPvalue <-
function(Mat, Ns, int, beta, delta){
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
  prob <- signif(prob, 12) #Remove rounding errors
  
  np <- int[which(prob==max(prob, na.rm=TRUE))]
  pvalue <- max(prob, na.rm=TRUE) + beta
  
  return(list(prob=prob, pvalue=pvalue, np=np))
}
