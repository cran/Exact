moreExtremePairedCSM <-
function(data, N, alternative, int, delta, reject.alpha){
  
  # There are cases where maxPvalue is faster than using the lookup tables in maxPvalueLookup
  lookupArray <- trinomCalc(N, int, delta)
  
  # Initialize first extreme table and considered tables:
  moreExtremeMat <- matrix(NA, N+1, N+1, dimnames=list(0:N, 0:N))
  moreExtremeMat <- lowerMatVal(moreExtremeMat)
  prevMoreExtremeMat <- moreExtremeMat
  
  moreExtremeMat[1, N+1] <- 1
  # Faster to search for 0's than NA's
  moreExtremeMat[1, N] <- 0
  
  # If two.sided and delta = 0, then opposite side is also most extreme
  # If two.sided and delta != 0, then opposite side must be a candidate to be picked
  if (alternative == "two.sided") { moreExtremeMat[N+1, 1] <- (delta == 0) }
  
  # If data is the most extreme table then simply return most extreme table
  if (!is.null(data)) {
    addRow <- which(moreExtremeMat == 1, arr.ind = TRUE)
    for (j in 1:nrow(addRow)) {
      if (all(addRow[j, ]-1 == c(data[1,2], data[2,1]))) {
        moreExtremeMat[is.na(moreExtremeMat)] <- 0
        return(list(TXO=NA, moreExtremeMat=moreExtremeMat))
      }
    }
  }
  
  # The two.sided and delta != 0 has to be programmed differently
  if (alternative == "two.sided" && delta != 0) {
    moreExtremeMat <- csmPairedTemp2sidedDelta(data, moreExtremeMat, N, int, alternative, lookupArray, delta, reject.alpha, checkPrev = TRUE, prevMoreExtremeMat)
  } else {
    moreExtremeMat <- csmPairedTemp(data, moreExtremeMat, N, int, alternative, lookupArray, delta, reject.alpha, checkPrev = TRUE, prevMoreExtremeMat)
  }
  # If moreExtremeMat is just FALSE, then return FALSE
  if (isFALSE(moreExtremeMat)) { return(moreExtremeMat) }
  # If reject.alpha is NULL, replace NA to 0
  moreExtremeMat[is.na(moreExtremeMat)] <- 0
  moreExtremeMat <- lowerMatVal(moreExtremeMat, NA)
  return(list(TXO=NA, moreExtremeMat=moreExtremeMat))
}
