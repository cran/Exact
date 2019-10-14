dbinomCalc <-
function(Ns, int, delta) {
  lookupArray <- list(matrix(nrow=Ns[1] + 1, ncol=length(int)),
                      matrix(nrow=Ns[2] + 1, ncol=length(int)))
  cellOne <- rep.int(1:(Ns[1]+1), length(int))+(Ns[1]+1)*rep(seq_along(int)-1, each=Ns[1]+1)
  lookupArray[[1]][cellOne] <- suppressWarnings(dbinom(rep.int(0:Ns[1], length(int)), Ns[1], rep(int+delta, each=Ns[1]+1)))
  
  cellTwo <- rep.int(1:(Ns[2]+1), length(int))+(Ns[2]+1)*rep(seq_along(int)-1, each=Ns[2]+1)
  lookupArray[[2]][cellTwo] <- suppressWarnings(dbinom(rep.int(0:Ns[2], length(int)), Ns[2], rep(int, each=Ns[2]+1)))
  
  return(lookupArray)
}
