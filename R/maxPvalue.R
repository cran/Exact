maxPvalue <-
function(tbls, Ns, npNumbers, int, beta){
  #Instead of calculating the binomial probabilities several times for a more extreme cell,
  #calculate the probability once and then combine with other more extreme cells:
  xTbls <- tbls[,1]
  yTbls <- tbls[,2]
  nTbls <- length(xTbls) # Number of 'as or more extreme' tables
  
  x.unique <- unique(xTbls)
  y.unique <- unique(yTbls)
  lxu <- length(x.unique)
  lyu <- length(y.unique)
  
  xnr <- max(xTbls)+1
  A <- matrix(nrow=xnr, ncol=npNumbers)
  cellnr <- rep.int(x.unique+1, npNumbers)+xnr*rep(seq(npNumbers)-1, each=lxu)
  A[cellnr] <- dbinom(rep.int(x.unique, npNumbers), Ns[1], rep(int, each=lxu))
  
  ynr <- max(yTbls)+1
  B <- matrix(nrow=ynr, ncol=npNumbers)
  cellnr <- rep.int(y.unique+1, npNumbers)+ynr*rep(seq(npNumbers)-1, each=lyu)
  B[cellnr] <- dbinom(rep.int(y.unique, npNumbers), Ns[2], rep(int, each=lyu))
  
  prob <- matrix(A[xTbls+1,]*B[yTbls+1,], ncol=length(int))
  prob <- colSums(prob, nTbls)
  
  np <- int[which(prob==max(prob))]
  pvalue <- max(prob) + beta
  
  return(list(prob=prob, pvalue=pvalue, np=np))
}
