csmTemp2sidedDelta <-
function(data, moreExtremeMat, Ns, int, alternative, lookupArray, delta, reject.alpha){
  
  AC <- which(moreExtremeMat==0, arr.ind = TRUE) - 1
  
  AC <- cbind(AC, (AC[,1]/Ns[1] - AC[,2]/Ns[2]) < delta)
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
    CcondAC[j] <- maxPvalueLookup(rbind(Tbls, AC[j,]),
                                  int=int, lookupArray=lookupArray)$pvalue
  }
  
  smallestPvalue <- min(round(CcondAC, digits=12))
  addRow <- AC[which(round(CcondAC, digits=12) == smallestPvalue), , drop=FALSE] + 1
  
  if (!is.null(reject.alpha) && smallestPvalue > reject.alpha) {
    
    # If looking at a specific dataset, then just return FALSE; otherwise, trying to form rejection region
    if (!is.null(data)) { return(FALSE) }
    
    # Return moreExtremeMat.  If you haven't added any tables, you have to check the most extreme table is still rejected
    if (sum(moreExtremeMat, na.rm=TRUE) == 1 && 
        maxPvalueLookup(Tbls, int=int, lookupArray=lookupArray)$pvalue > reject.alpha) {
      moreExtremeMat[which(moreExtremeMat==1, arr.ind = TRUE)] <- 0
    }
    return(moreExtremeMat)
  }
  
  # Update moreExtremeMat
  for (j in 1:nrow(addRow)) {
    moreExtremeMat[addRow[j,1], addRow[j,2]] <- 1
    #if (alternative == 'two.sided') { moreExtremeMat[Ns[1] + 2 - addRow[j,1], Ns[2] + 2 - addRow[j,2]] <- 1 }
    if (addRow[j,2] >= 2 && is.na(moreExtremeMat[addRow[j,1], addRow[j,2]-1])) { moreExtremeMat[addRow[j,1], addRow[j,2]-1] <- 0 }
    if (addRow[j,1] <= Ns[1] && is.na(moreExtremeMat[addRow[j,1]+1, addRow[j,2]])) { moreExtremeMat[addRow[j,1]+1, addRow[j,2]] <- 0 }
    if (addRow[j,2] <= Ns[2] && is.na(moreExtremeMat[addRow[j,1], addRow[j,2]+1])) { moreExtremeMat[addRow[j,1], addRow[j,2]+1] <- 0 }
    if (addRow[j,1] >= 2 && is.na(moreExtremeMat[addRow[j,1]-1, addRow[j,2]])) { moreExtremeMat[addRow[j,1]-1, addRow[j,2]] <- 0 }
  }
  
  if (!is.null(data)) {
    for (j in 1:nrow(addRow)) {
      if (all(addRow[j, ]-1 == data[1,])) { return(moreExtremeMat) }
    }
  }
  
  #Perform recursive loop
  csmTemp2sidedDelta(data, moreExtremeMat, Ns, int, alternative, lookupArray, delta, reject.alpha)
}
