santner_TX <-
function(Ns){
  N <- sum(Ns)
  x <- rep(0:Ns[1], each=(Ns[2]+1))
  y <- rep.int(0:Ns[2], Ns[1]+1)
  p1 <- x/Ns[1]
  p2 <- y/Ns[2]
  TX <- matrix(c(x, y, p1-p2), nrow=(Ns[1]+1)*(Ns[2]+1), ncol=3)
}