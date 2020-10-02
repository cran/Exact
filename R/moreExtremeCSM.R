moreExtremeCSM <-
function(data, Ns, alternative, int, doublePvalue, delta, reject.alpha){
  
  # There are cases where maxPvalue is faster than using the lookup tables in maxPvalueLookup
  lookupArray <- dbinomCalc(Ns, int, delta)
  
  # Initialize first extreme table and considered tables:
  moreExtremeMat <- prevMoreExtremeMat <- matrix(NA, Ns[1]+1, Ns[2]+1, dimnames=list(0:Ns[1], 0:Ns[2]))
  
  moreExtremeMat[1, Ns[2]+1] <- 1
  # Faster to search for 0's than NA's
  moreExtremeMat[1, Ns[2]] <- 0
  moreExtremeMat[2, Ns[2]+1] <- 0
  
  # If two.sided and delta = 0, then opposite side is also most extreme
  # If two.sided and delta != 0, then opposite side must be a candidate to be picked
  if (alternative == "two.sided") { moreExtremeMat[Ns[1]+1, 1] <- (delta == 0) }
  
  # If data is the most extreme table then simply return most extreme table
  if (!is.null(data)) {
    addRow <- which(moreExtremeMat == 1, arr.ind = TRUE)
    for (j in 1:nrow(addRow)) {
      if (all(addRow[j, ]-1 == data[1,])) {
        moreExtremeMat[is.na(moreExtremeMat)] <- 0
        return(list(TXO=NA, moreExtremeMat=moreExtremeMat))
      }
    }
  }
  
  # The two.sided and delta != 0 has to be programmed differently
  if (alternative == "two.sided" && delta != 0) {
    moreExtremeMat <- csmTemp2sidedDelta(data, moreExtremeMat, Ns, int, alternative, lookupArray, doublePvalue, delta, reject.alpha, checkPrev = TRUE, prevMoreExtremeMat)
  } else {
    moreExtremeMat <- csmTemp(data, moreExtremeMat, Ns, int, alternative, lookupArray, doublePvalue, reject.alpha, checkPrev = TRUE, prevMoreExtremeMat)
  }
  # If moreExtremeMat is just FALSE, then return FALSE
  if (isFALSE(moreExtremeMat)) { return(moreExtremeMat) }
  # If reject.alpha is NULL, replace NA to 0
  moreExtremeMat[is.na(moreExtremeMat)] <- 0
  return(list(TXO=NA, moreExtremeMat=moreExtremeMat))
}
