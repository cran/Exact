zunpooled_TX <-
function(data, Ns) {
  N <- sum(Ns)
  if (!is.null(data)) {
    p1 <- data[1,1]/Ns[1]
    p2 <- data[1,2]/Ns[2]
    TX <- (p1-p2)/sqrt(p2*(1-p2)/Ns[2]+(p1)*(1-p1)/Ns[1])
  } else {
    x <- rep(0:Ns[1], each=(Ns[2]+1))
    y <- rep.int(0:Ns[2], Ns[1]+1)
    p1 <- x/Ns[1]
    p2 <- y/Ns[2]
    TX <- matrix(c(x, y, (p1-p2)/sqrt(p2*(1-p2)/Ns[2]+(p1)*(1-p1)/Ns[1])), nrow=(Ns[1]+1)*(Ns[2]+1), ncol=3)
  }
  return(TX)
}
