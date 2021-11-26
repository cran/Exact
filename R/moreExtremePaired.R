moreExtremePaired <-
function(method, data, N, alternative, int, delta){
  
  # If alternative="less", then only need to calculate boundary in upper triangle (faster than calculating all test statistics).
  # If alternative="two.sided" and delta=0, then can also calculate boundary in upper triangle since the lower triangle is always symmetric
  # If alternative="two.sided" and delta != 0, then unfortunately, need to consider every table.
  # Special case for UCM: Conditional McNemar's p-value doesn't include a delta term, so use same ordering ignoring delta
  # (note: this can't work if delta != 0 and two-sided)
  
  if (alternative != "two.sided" || delta == 0) {
    
    # Important note: ignore delta when ordering "ucm"
    TXO <- switch(method, 
                  "uam" = mcnemar_TX(data, N, delta = delta, CC=FALSE),
                  "ucm" = mcnemar.2x2(data, N, alternative=alternative, pval1ties = TRUE),
                  "uamcc" = mcnemar_TX(data, N, delta = delta, CC=TRUE))[3]
    
    # Doesn't appear this check is needed, but useful to confirm if debugging
    if (is.na(TXO)) { stop("Test statistic is NA; please check code") }
    TXO <- signif(TXO, 12)  #Remove rounding errors
    
    moreExtremeMat <- matrix(NA, N+1, N+1, dimnames=list(0:N, 0:N))
    for (j in (data[2,1]+1):(N+1)) { moreExtremeMat[0:(data[1,2]+1), j] <- 1 }
    
    # Change values that are impossible to 999
    moreExtremeMat <- lowerMatVal(moreExtremeMat, 999)
    
    for (i in 0:N) {  #Go through each row
      # Find first column with NA
      startJ <- which(is.na(moreExtremeMat[i+1, ]))[1] - 1
      if (!is.na(startJ)) {
        for (j in startJ:N) {
          newDat <- matrix(c(N-i-j, i, j, 0), 2, 2, byrow=TRUE)
          newTX <- switch(method,
                          "uam" = mcnemar_TX(newDat, N, delta=delta, CC=FALSE),
                          "ucm" = mcnemar.2x2(newDat, N, alternative=alternative, pval1ties = TRUE),
                          "uamcc" = mcnemar_TX(newDat, N, delta=delta, CC=TRUE))[3]
          
          # Doesn't appear this check is needed, but useful to confirm if debugging
          if (is.na(newTX)) { stop("Test statistic is NA; please check code") }
          
          newTX <- signif(newTX, 12)  #Remove rounding errors
          
          if (method %in% c("uam", "uamcc")) {
            rejectFlg <- switch(alternative,
                                "less" = (newTX <= TXO),
                                "two.sided" = (abs(newTX) >= abs(TXO)))
          } else if (method == "ucm") {
            rejectFlg <- (newTX <= TXO)
          }
          
          if (rejectFlg) {
            # If more extreme test statistic, then know the remaining columns in the row is more extreme
            moreExtremeMat[i+1, (j+1):(N+1)] <- 1
            break
          } else {
            # If less extreme test statistic, then know the remaining rows in the column is less extreme
            moreExtremeMat[(i+1):(N+1), j+1] <- 0
          }
        }
      }
    }
    
    # If delta=0, then two-sided test will be symmetric
    if (alternative == "two.sided") {
      rejectTemp <- which(moreExtremeMat==1, arr.ind = TRUE)
      moreExtremeMat[rejectTemp[ , 2:1, drop=FALSE]] <- 1
    }
    
    # Change values that are impossible to NA at the end
    moreExtremeMat <- lowerMatVal(moreExtremeMat, NA)
    
  } else {  #The only case where we can't just calculate the boundary is if we have a two-sided test with delta != 0
    
    TX <- switch(method, 
                 "uam" = mcnemar_TX(NULL, N, delta=delta, CC=FALSE),
                 "ucm" = mcnemar.2x2(NULL, N, alternative=alternative, pval1ties = TRUE),
                 "uamcc" = mcnemar_TX(NULL, N, delta=delta, CC=TRUE))
    
    # Doesn't appear this check is needed, but useful to confirm if debugging
    if (any(is.na(TX[ , 3]))) { stop("Test statistic is NA; please check code") }
    TX[, 3] <- signif(TX[ , 3], 12)  #Remove rounding errors
    
    TXO <- TX[TX[,1]==data[1,2] & TX[,2]==data[2,1], 3]
    
    # Note: alternative must be "two.sided"
    if (method %in% c("uam", "uamcc")) { rejectFlg <- (abs(TX[,3]) >= abs(TXO))
    } else { rejectFlg <- (TX[,3] <= TXO) }
    
    moreExtremeMat <- matrix(0, N+1, N+1, dimnames=list(0:N, 0:N))
    moreExtremeMat <- lowerMatVal(moreExtremeMat, NA)
    moreExtremeMat[!is.na(moreExtremeMat)] <- rejectFlg
  }
  return(list(TXO=TXO, moreExtremeMat=moreExtremeMat))
}
