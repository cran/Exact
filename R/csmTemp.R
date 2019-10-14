csmTemp <-
function(data, moreExtremeMat, Ns, int, alternative, lookupArray, reject.alpha){
  #If observed AC is in Tbls, then stop
  #if (!is.null(data) && sum(apply(Tbls, 1, function(x) all(x == data[1,])))) {return(Tbls)}
  
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
      CcondAC[j] <- maxPvalueLookup(rbind(Tbls, AC[j,], c(Ns[1]-AC[j,1], Ns[2]-AC[j,2])),
                                    int=int, lookupArray=lookupArray)$pvalue
    } else {
      CcondAC[j] <- maxPvalueLookup(rbind(Tbls, AC[j,]),
                                    int=int, lookupArray=lookupArray)$pvalue
    }
  }
  
  smallestPvalue <- min(round(CcondAC, digits=12))
  addRow <- AC[which(round(CcondAC, digits=12) == smallestPvalue), , drop=FALSE] + 1
  
  if (!is.null(reject.alpha) && smallestPvalue > reject.alpha) {
    
    # If looking at a specific dataset, then just return FALSE; otherwise, trying to form rejection region
    if (!is.null(data)) { return(FALSE) }
    
    # Return moreExtremeMat.  If you haven't added any tables, you have to check the most extreme table is still rejected
    if (alternative=="two.sided" && sum(moreExtremeMat, na.rm=TRUE) == 2 && 
        maxPvalueLookup(Tbls, int=int, lookupArray=lookupArray)$pvalue > reject.alpha) {
      moreExtremeMat[which(moreExtremeMat==1, arr.ind = TRUE)] <- 0
    } else if (sum(moreExtremeMat, na.rm=TRUE) == 1 && 
               maxPvalueLookup(Tbls, int=int, lookupArray=lookupArray)$pvalue > reject.alpha) {
      moreExtremeMat[which(moreExtremeMat==1, arr.ind = TRUE)] <- 0
    }
    return(moreExtremeMat)
  }
  
  # Update moreExtremeMat
  for (j in 1:nrow(addRow)) {
    moreExtremeMat[addRow[j,1], addRow[j,2]] <- 1
    if (alternative == 'two.sided') { moreExtremeMat[Ns[1] + 2 - addRow[j,1], Ns[2] + 2 - addRow[j,2]] <- 1 }
    if (length(moreExtremeMat[addRow[j,1], addRow[j,2]-1]) > 0 &&
        is.na(moreExtremeMat[addRow[j,1], addRow[j,2]-1])) { moreExtremeMat[addRow[j,1], addRow[j,2]-1] <- 0 }
    if (addRow[j,1] <= Ns[1] && length(moreExtremeMat[addRow[j,1]+1, addRow[j,2]]) > 0 &&
        is.na(moreExtremeMat[addRow[j,1]+1, addRow[j,2]])) { moreExtremeMat[addRow[j,1]+1, addRow[j,2]] <- 0 }
  }
  
  if (!is.null(data)) {
    for (j in 1:nrow(addRow)) {
      if (all(addRow[j, ]-1 == data[1,])) { return(moreExtremeMat) }
    }
  }
  
  #Perform recursive loop
  csmTemp(data, moreExtremeMat, Ns, int, alternative, lookupArray, reject.alpha)
}
