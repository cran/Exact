csmApprox_TX <-
function(data, Ns, alternative, int, lookupArray){
  
  if (!is.null(data)) {
    x <- data[1,1]
    y <- data[1,2]
  } else {
    x <- rep(0:Ns[1], each=(Ns[2]+1))
    y <- rep.int(0:Ns[2], Ns[1]+1)
  }

  Tbls <- cbind(x, y, deparse.level=0)
  TX <- apply(Tbls, 1, function(x) {
    if (alternative == "two.sided" && (x[1] != Ns[1]-x[1] || x[2] != Ns[2]-x[2])) {
      x <- rbind(x, c(Ns[1]-x[1], Ns[2]-x[2]))
    } else { x <- rbind(x) }
    maxPvalueLookup(x, int=int, lookupArray=lookupArray, doublePvalue=FALSE)$pvalue
  })
  
  return(cbind(x, y, TX, deparse.level=0))
}
