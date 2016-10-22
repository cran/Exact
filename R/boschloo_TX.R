boschloo_TX <-
function(Ns, alternative){
  N <- sum(Ns)
  x <- rep(0:Ns[1], each=(Ns[2]+1))
  y <- rep.int(0:Ns[2], Ns[1]+1)
  p1 <- x/Ns[1]
  p2 <- y/Ns[2]
  pval <- apply(matrix(c(x, Ns[1]-x, y, Ns[2]-y), (Ns[1]+1)*(Ns[2]+1), 4), 1,
                FUN=function(tbls){fisher.2x2(matrix(tbls,2,2), alternative=alternative)})
  return(matrix(c(x, y, pval), nrow=(Ns[1]+1)*(Ns[2]+1), ncol=3))
}
