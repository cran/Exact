zpooled_TX <-
function(Ns){
  N <- sum(Ns)
  x <- rep(0:Ns[1], each=(Ns[2]+1))
  y <- rep.int(0:Ns[2], Ns[1]+1)
  p1 <- x/Ns[1]
  p2 <- y/Ns[2]
  return(matrix(c(x, y, (p1-p2)/sqrt(((x+y)/N)*(1-((x+y)/N))*sum(1/Ns))), nrow=(Ns[1]+1)*(Ns[2]+1), ncol=3))
}
