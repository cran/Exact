csmPairedTemp2sidedDelta <-
function(data, moreExtremeMat, N, int, alternative, lookupArray, doublePvalue, delta, reject.alpha, checkPrev, prevMoreExtremeMat){
  
  # Only use Asymptotic McNemar for ties
  TX <- mcnemar_TX(NULL, N, delta=delta, CC=FALSE)
  TX[, 3] <- signif(TX[ , 3], 12)  #Remove rounding errors
  TX <- TX[order(TX[,1], TX[,2]), ]
  TX[,3] <- -abs(TX[,3])
  nIter <- 1
  
  # Use a while loop (instead of recursive loop) to prevent error: "node stack overflow; no more error handlers available" #
  while (TRUE) {
    
    AC <- which(moreExtremeMat==0, arr.ind = TRUE) - 1
    
    AC <- cbind(AC, (AC[,1]/N - AC[,2]/N) < delta)
    AC_LT <- AC[AC[,3] == 0, 1:2, drop=FALSE]
    AC_LT <- AC_LT[order(AC_LT[,1],AC_LT[,2]), , drop=FALSE]
    AC_LT <- AC_LT[!duplicated(AC_LT[,1]), , drop=FALSE]
    AC_LT <- AC_LT[order(AC_LT[,2],-AC_LT[,1]), , drop=FALSE]
    AC_LT <- AC_LT[!duplicated(AC_LT[,2]), , drop=FALSE]
    
    AC_UT <- AC[AC[,3] == 1, 1:2, drop=FALSE]
    AC_UT <- AC_UT[order(AC_UT[,1],-AC_UT[,2]), , drop=FALSE]
    AC_UT <- AC_UT[!duplicated(AC_UT[,1]), , drop=FALSE]
    AC_UT <- AC_UT[order(AC_UT[,2],AC_UT[,1]), , drop=FALSE]
    AC_UT <- AC_UT[!duplicated(AC_UT[,2]), , drop=FALSE]
    
    AC <- rbind(AC_LT, AC_UT)
    
    #Calculate the possible more extreme test statistic:
    Tbls <- which(moreExtremeMat==1, arr.ind = TRUE) - 1
    
    CcondAC <- rep(0, nrow(AC))
    for (j in 1:nrow(AC)) {
      CcondAC[j] <- maxPvaluePairedLookup(rbind(Tbls, AC[j,]),
                                    int=int, lookupArray=lookupArray, doublePvalue=doublePvalue)$pvalue
    }
    
    smallestPvalue <- min(round(CcondAC, digits=12))
    
    if (!is.null(reject.alpha) && smallestPvalue > reject.alpha) {
      
      # If looking at a specific dataset, then just return FALSE; otherwise, trying to form rejection region
      if (!is.null(data)) { return(FALSE) }
      
      # There are 2 cases where moreExtremeMat may be incorrect and needs to be updated:
      # (1) if no tables have been added and even most extreme table is not significant (unlikely)
      # (2) if previously added two tables where individually the p-values are < alpha, but together are larger than alpha (possible)
      if (checkPrev && maxPvaluePairedLookup(Tbls, int=int, lookupArray=lookupArray, doublePvalue=doublePvalue)$pvalue > reject.alpha) {
        moreExtremeMat <- prevMoreExtremeMat
      }
      return(moreExtremeMat)
    }
    
    # Update moreExtremeMat
    addRow <- AC[which(round(CcondAC, digits=12) == smallestPvalue), , drop=FALSE] + 1
    # If there are ties, use McNemar's test to break ties
    if (nrow(addRow) > 1) {
      TXties <- cbind(addRow, apply(addRow, 1, function(x) { TX[TX[ , 1] == (x[1]-1) & TX[ , 2] == (x[2]-1), 3] }))
      TXties <- TXties[order(TXties[,3]), ]
      addRow <- TXties[TXties[ , 3] <= TXties[1,3], 1:2, drop=FALSE]
    }
    
    checkPrev <- (nrow(addRow) > 1)
    prevMoreExtremeMat <- moreExtremeMat
    
    moreExtremeMat <- updateMat(moreExtremeMat, addRow)
    
    if (!is.null(data)) {
      for (j in 1:nrow(addRow)) {
        if (all(addRow[j, ]-1 == c(data[1,2], data[2,1]))) {
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
  #csmPairedTemp2sidedDelta(data, moreExtremeMat, N, int, alternative, lookupArray, doublePvalue, delta, reject.alpha, checkPrev, prevMoreExtremeMat)
}
