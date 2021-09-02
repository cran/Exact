mcnemar_TX <-
function(data, N, delta, CC) {
  
  stopifnot(is.logical(CC))
  
  if (!is.null(data)) {
    x <- data[1,2]
    y <- data[2,1]
  } else {
    # Find all discordant possibilities #
    x <- unlist(lapply(0:N, function(z) {0:(N-z)}))
    y <- rep(0:N, N-0:N+1)
  }
  
  # See Tango (1998) papers to calculate Z-statistic with non-zero delta.  Reverse sign of delta
  delta <- -delta
  numerator <- x - y + N*delta
  
  if (CC) {
    # If numerator is between -1 and 1, set to 0
    numerator[-1 < numerator & numerator < 1] <- 0
    # If numerator is >= 1, then subtract 1
    numerator[numerator >= 1] <- numerator[numerator >= 1] - 1
    # If numerator is <= -1, then add 1
    numerator[numerator <= -1] <- numerator[numerator <= -1] + 1
  }
  
  if (delta==0) {
    denominator <- sqrt(x + y)
  } else {
    A <- 2*N
    B <- -x - y - (2*N-x+y)*delta
    C <- y*delta*(delta + 1)
    # Due to floating point issues, it's possible B^2 - 4*A*C is very slightly <0
    inSqrtRoot <- B^2 - 4*A*C
    inSqrtRoot[inSqrtRoot < 0] <- 0
    q21D <- (-B + sqrt(inSqrtRoot)) / (2*A)
    denominator <- sqrt( N*(2*q21D - delta*(delta+1)) )
  }
  
  TX <- numerator / denominator
  TX[numerator == 0 & denominator == 0] <- 0
  return(cbind(x, y, TX, deparse.level=0))
}
