searchExtreme <-
function(TX, n1, n2, alternative, method, int, delta, alpha, lookupArray) {
  TXunique <- unique(TX[is.na(TX[ , 4]), 3])
  if (length(TXunique) == 0) {return(TX[as.logical(TX[ , 4]), 1:2, drop=FALSE])}
  m <- floor(length(TXunique)/2) + 1  #Very slightly faster
  s <- TXunique[m]
  
  if (method %in% c("chisq", "yates chisq", "fisher")) {
    Tbls <- TX[TX[,3] == s, , drop=FALSE][1, 1:2]
    Tbls <- matrix(c(Tbls[1],n1-Tbls[1],Tbls[2],n2-Tbls[2]), byrow=TRUE, ncol=2)
    if (method == "fisher") { pvalue <- fisher.2x2(Tbls, NULL, alternative=alternative)
    } else {
      pvalue <- suppressWarnings(prop.test(Tbls, alternative=alternative, correct=(method=="yates chisq"))$p.value)
    }
  } else {
    Tbls <- TX[TX[,3] <= s, 1:2, drop=FALSE]
    pvalue <- maxPvalueLookup(Tbls, int=int, lookupArray=lookupArray)$pvalue
  }
  if (pvalue <= alpha){ TX[TX[,3] <= s, 4] <- TRUE
  } else { TX[TX[,3] >= s, 4] <- FALSE }
  
  return(searchExtreme(TX = TX, n1 = n1, n2 = n2, alternative = alternative, method = method, int = int, delta = delta, alpha = alpha, lookupArray = lookupArray))
}
