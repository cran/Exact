chisq_TX <-
function(data, Ns, yates) {
  stopifnot(is.logical(yates))
  N <- sum(Ns)
  if (!is.null(data)) {
    x <- data[1,1]
    y <- data[1,2]
    if (!yates) {
      TX <- N*(x*(Ns[2]-y) - y*(Ns[1]-x))^2/(Ns[1]*Ns[2]*(x+y)*(N-x-y))
    } else {
      sqTerm <- abs(x*(Ns[2]-y) - y*(Ns[1]-x)) - 1/2*N
      sqTerm[sqTerm < 0] <- 0
      TX <- N*(sqTerm)^2/(Ns[1]*Ns[2]*(x+y)*(N-x-y))
    }
  } else {
    N <- sum(Ns)
    x <- rep(0:Ns[1], each=(Ns[2]+1))
    y <- rep.int(0:Ns[2], Ns[1]+1)
    # Note: when Yates is false, this statistic is equivalent to the squared z-pooled statistic
    if (!yates) {
      TX <- matrix(c(x,y,N*(x*(Ns[2]-y) - y*(Ns[1]-x))^2/(Ns[1]*Ns[2]*(x+y)*(N-x-y))),(Ns[1]+1)*(Ns[2]+1),3)
    } else {
      sqTerm <- abs(x*(Ns[2]-y) - y*(Ns[1]-x)) - 1/2*N
      sqTerm[sqTerm < 0] <- 0
      TX <- matrix(c(x,y,N*(sqTerm)^2/(Ns[1]*Ns[2]*(x+y)*(N-x-y))),(Ns[1]+1)*(Ns[2]+1),3)
    }
  }
  return(TX)
}
