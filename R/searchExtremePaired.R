searchExtremePaired <-
function(TX, N, alternative, method, int, delta, alpha, lookupArray) {
  TXunique <- unique(TX[is.na(TX[ , 4]), 3])
  if (length(TXunique) == 0) {return(TX[as.logical(TX[ , 4]), 1:2, drop=FALSE])}
  m <- floor(length(TXunique)/2) + 1  #Very slightly faster
  s <- TXunique[m]
  
  if (method %in% c("am", "amcc")) {
    # Many ways to calculate Asymptotic McNemar test, but use z-distribution for one-sided tests
    pvalue <- ifelse(alternative=="two.sided",2,1)*pnorm(TX[TX[,3] == s, 3][[1]])
    
    #Tbls <- TX[TX[,3] == s, , drop=FALSE][1, 1:2]
    #Tbls <- matrix(c(N-Tbls[1]-Tbls[2],Tbls[1],Tbls[2],0), byrow=TRUE, ncol=2)
    #mcnemar.test(Tbls, correct=FALSE)$p.value
  } else {
    Tbls <- TX[TX[,3] <= s, 1:2, drop=FALSE]
    pvalue <- maxPvaluePairedLookup(Tbls, int=int, lookupArray=lookupArray, doublePvalue=FALSE)$pvalue
  }  
  if (pvalue <= alpha){ TX[TX[,3] <= s, 4] <- TRUE
  } else { TX[TX[,3] >= s, 4] <- FALSE }
  
  return(searchExtremePaired(TX = TX, N = N, alternative = alternative, method = method, int = int, delta = delta, alpha = alpha,
                             lookupArray = lookupArray))
}
