maxPvalueLookup <-
function(tbls, int, lookupArray){
  #Instead of calculating the binomial probabilities several times for a more extreme cell,
  #calculate the probability once and then combine with other more extreme cells:
  xTbls <- tbls[,1]
  yTbls <- tbls[,2]
  nTbls <- length(xTbls) # Number of 'as or more extreme' tables
  
  prob <- matrix(lookupArray[[1]][xTbls+1,]*lookupArray[[2]][yTbls+1,], ncol=length(int))
  prob <- .colSums(prob, nTbls, length(int))
  prob <- signif(prob, 12) #Remove rounding errors
  
  np <- int[which(prob==max(prob, na.rm=TRUE))]
  pvalue <- max(prob, na.rm=TRUE)
  
  return(list(prob=prob, pvalue=pvalue, np=np))
}
