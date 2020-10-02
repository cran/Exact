zunpooled_TX <-
function(data, Ns, delta) {
  N <- sum(Ns)
  if (!is.null(data)) {
    x <- data[1,1]
    y <- data[1,2]
  } else {
    x <- rep(0:Ns[1], each=(Ns[2]+1))
    y <- rep.int(0:Ns[2], Ns[1]+1)
  }
  p1 <- x/Ns[1]
  p2 <- y/Ns[2]
  
  numerator <- p1 - p2 - delta
  denominator <- sqrt(p2*(1-p2)/Ns[2]+(p1)*(1-p1)/Ns[1])
  TX <- numerator / denominator
  TX[numerator == 0 & denominator == 0] <- 0
  return(cbind(x, y, TX, deparse.level=0))
}
