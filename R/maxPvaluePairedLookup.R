maxPvaluePairedLookup <-
function(tbls, int, lookupArray){
  #Instead of calculating the trinomial probabilities several times for a more extreme cell,
  #calculate the probability once and then combine with other more extreme cells:
  prob <- matrix(0, ncol=length(int))
  for (row in 1:nrow(tbls)) {
    prob <- prob + lookupArray[[tbls[row,1]+1]][tbls[row,2]+1, ]
  }
  
  prob <- signif(prob, 12) #Remove rounding errors
  
  np <- int[which(prob==max(prob, na.rm=TRUE))]
  pvalue <- max(prob, na.rm=TRUE)
  
  return(list(prob=prob, pvalue=pvalue, np=np))
}
