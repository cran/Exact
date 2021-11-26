maxPvaluePairedLookup <-
function(tbls, int, lookupArray, doublePvalue){
  #Instead of calculating the trinomial probabilities several times for a more extreme cell,
  #calculate the probability once and then combine with other more extreme cells:
  prob <- matrix(0, ncol=length(int))
  for (row in 1:nrow(tbls)) {
    prob <- prob + lookupArray[[tbls[row,1]+1]][tbls[row,2]+1, ]
  }
  
  prob <- signif((1+doublePvalue)*prob, 12) #Remove rounding errors and double probabilities if applicable
  
  np <- int[which(prob==max(prob, na.rm=TRUE))]
  pvalue <- max(prob, na.rm=TRUE)
  
  return(list(prob=prob, pvalue=pvalue, np=np))
}
