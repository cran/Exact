chisq_TX <-
function(data, Ns, yates) {
  stopifnot(is.logical(yates))
  
  if (!is.null(data)) {
    x <- data[1,1]
    y <- data[1,2]
  } else {
    x <- rep(0:Ns[1], each=(Ns[2]+1))
    y <- rep.int(0:Ns[2], Ns[1]+1)
  }
  N <- sum(Ns)
  # Note: when Yates is false, this statistic is equivalent to the squared z-pooled statistic
  if (!yates) {
    numerator <- N*(x*(Ns[2]-y) - y*(Ns[1]-x))^2
  } else {
    sqTerm <- abs(x*(Ns[2]-y) - y*(Ns[1]-x)) - 1/2*N
    sqTerm[sqTerm < 0] <- 0
    numerator <- N*(sqTerm)^2
  }
  denominator <- Ns[1]*Ns[2]*(x+y)*(N-x-y)
  
  TX <- numerator / denominator
  TX[numerator == 0 & denominator == 0] <- 0
  return(cbind(x, y, TX, deparse.level=0))
}
