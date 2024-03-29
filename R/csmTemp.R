csmTemp <-
function(data, moreExtremeMat, Ns, int, alternative, lookupArray, doublePvalue, delta, reject.alpha, checkPrev, prevMoreExtremeMat){
  #If observed AC is in Tbls, then stop
  #if (!is.null(data) && sum(apply(Tbls, 1, function(x) all(x == data[1,])))) {return(Tbls)}
  
  # Only use z-pooled for ties
  TX <- zpooled_TX(NULL, Ns, delta)
  TX[, 3] <- signif(TX[ , 3], 12)  #Remove rounding errors
  TX <- TX[order(TX[,1], TX[,2]), ]
  if (alternative == "two.sided") { TX[,3] <- -abs(TX[,3]) }
  nIter <- 1

  # for (i in 1:278) {
  # Use a while loop (instead of recursive loop) to prevent error: "node stack overflow; no more error handlers available" #
  while (TRUE) {
  
    #AC <- which(is.na(moreExtremeMat), arr.ind = TRUE) - 1
    AC <- which(moreExtremeMat==0, arr.ind = TRUE) - 1
    AC <- AC[order(AC[,1],-AC[,2]), , drop=FALSE]
    AC <- AC[!duplicated(AC[,1]), , drop=FALSE]
    AC <- AC[order(AC[,2],AC[,1]), , drop=FALSE]
    AC <- AC[!duplicated(AC[,2]), , drop=FALSE]
    
    #Calculate the possible more extreme test statistic:
    Tbls <- which(moreExtremeMat==1, arr.ind = TRUE) - 1
    
    CcondAC <- rep(0, nrow(AC))
    for (j in 1:nrow(AC)) {
      if (alternative == 'two.sided') {
        if (all(AC[j,] == c(Ns[1]-AC[j,1], Ns[2]-AC[j,2]))) {
          CcondAC[j] <- maxPvalueLookup(rbind(Tbls, AC[j,]),
                                        int=int, lookupArray=lookupArray, doublePvalue=doublePvalue)$pvalue
        } else {
          CcondAC[j] <- maxPvalueLookup(rbind(Tbls, AC[j,], c(Ns[1]-AC[j,1], Ns[2]-AC[j,2])),
                                        int=int, lookupArray=lookupArray, doublePvalue=doublePvalue)$pvalue
        }
      } else {
        CcondAC[j] <- maxPvalueLookup(rbind(Tbls, AC[j,]),
                                      int=int, lookupArray=lookupArray, doublePvalue=doublePvalue)$pvalue
      }
    }
    
    smallestPvalue <- min(round(CcondAC, digits=12))
    
    if (!is.null(reject.alpha) && smallestPvalue > reject.alpha) {
      
      # If looking at a specific dataset, then just return FALSE; otherwise, trying to form rejection region
      if (!is.null(data)) { return(FALSE) }
      
      # There are 2 cases where moreExtremeMat may be incorrect and needs to be updated:
      # (1) if no tables have been added and even most extreme table is not significant (unlikely)
      # (2) if previously added two tables where individually the p-values are < alpha, but together are larger than alpha (possible)
      if (checkPrev && maxPvalueLookup(Tbls, int=int, lookupArray=lookupArray, doublePvalue=doublePvalue)$pvalue > reject.alpha) {
        moreExtremeMat <- prevMoreExtremeMat
      }
      
      return(moreExtremeMat)
    }
    
    # Update moreExtremeMat
    addRow <- AC[which(round(CcondAC, digits=12) == smallestPvalue), , drop=FALSE] + 1
    # If there are ties, use Z-test to break ties
    if (nrow(addRow) > 1) {
      TXties <- cbind(addRow, apply(addRow, 1, function(x) { TX[TX[ , 1] == (x[1]-1) & TX[ , 2] == (x[2]-1), 3] }))
      TXties <- TXties[order(TXties[,3]), ]
      addRow <- TXties[TXties[ , 3] <= TXties[1,3], 1:2, drop=FALSE]
    }
    
    checkPrev <- (nrow(addRow) > 1)
    prevMoreExtremeMat <- moreExtremeMat
    
    for (j in 1:nrow(addRow)) {
      moreExtremeMat[addRow[j,1], addRow[j,2]] <- 1
      if (alternative == 'two.sided') { moreExtremeMat[Ns[1] + 2 - addRow[j,1], Ns[2] + 2 - addRow[j,2]] <- 1 }
      if (length(moreExtremeMat[addRow[j,1], addRow[j,2]-1]) > 0 &&
          is.na(moreExtremeMat[addRow[j,1], addRow[j,2]-1])) { moreExtremeMat[addRow[j,1], addRow[j,2]-1] <- 0 }
      if (addRow[j,1] <= Ns[1] && length(moreExtremeMat[addRow[j,1]+1, addRow[j,2]]) > 0 &&
          is.na(moreExtremeMat[addRow[j,1]+1, addRow[j,2]])) { moreExtremeMat[addRow[j,1]+1, addRow[j,2]] <- 0 }
    }
    
    # Check if added row includes data
    if (!is.null(data)) {
      for (j in 1:nrow(addRow)) {
        if (all(addRow[j, ]-1 == data[1,]) || (alternative == "two.sided" && all(c(Ns[1] + 2 - addRow[j,1], Ns[2] + 2 - addRow[j,2])-1 == data[1,]))) {
          return(moreExtremeMat)
        }
      }
    }
    
    nIter <- nIter + 1
    if (nIter %% 5000 == 0) {
      print(paste0("CSM added ", nIter, " more extreme tables so far; may be too computationally intensive and suggest aborting"))
    }
  }
  
  #Perform recursive loop
  #csmTemp(data, moreExtremeMat, Ns, int, alternative, lookupArray, doublePvalue, delta, reject.alpha, checkPrev, prevMoreExtremeMat)
}
