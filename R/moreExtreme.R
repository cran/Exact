moreExtreme <-
function(method, data, Ns, alternative, int, delta){
  
  # If alternative="less", then only need to calculate boundary in upper triangle (faster than calculating all test statistics).
  # If alternative="two.sided" and delta=0, then can also calculate boundary in upper triangle since the lower triangle is always symmetric
  # If alternative="two.sided" and delta != 0, then unfortunately, need to consider every table.
  # Special case for Boschloo and Z-unpooled: test statistic doesn't include a delta term, so use same ordering ignoring delta
  # (note: this can't work if delta != 0 and two-sided)
  
  if (alternative != "two.sided" || delta == 0) {
    
    # Important note: ignore delta when ordering "boschloo" or "z-unpooled"
    TXO <- switch(method, 
                  "z-pooled" = zpooled_TX(data, Ns, delta=delta),
                  "z-unpooled" = zunpooled_TX(data, Ns),
                  "boschloo" = fisher.2x2(data, alternative=alternative),
                  "santner and snell" = santner_TX(data, Ns, delta=delta))[3]
    
    # Doesn't appear this check is needed, but useful to confirm if debugging
    if (is.na(TXO)) { stop("Test statistic is NA; please check code") }
    TXO <- signif(TXO, 12)  #Remove rounding errors
    
    moreExtremeMat <- matrix(NA, Ns[1]+1, Ns[2]+1, dimnames=list(0:Ns[1], 0:Ns[2]))
    for (j in (data[1,2]+1):(Ns[2]+1)) { moreExtremeMat[0:(data[1,1]+1), j] <- 1 }
    
    for (i in 0:Ns[1]) {  #Go through each row
      # Find first column with NA
      startJ <- which(is.na(moreExtremeMat[i+1, ]))[1] - 1
      if (!is.na(startJ)) {
        for (j in startJ:Ns[2]) {
          newDat <- matrix(c(i, Ns[1]-i, j, Ns[2]-j), 2, 2)
          newTX <- switch(method,
                          "z-pooled" = zpooled_TX(newDat, Ns, delta=delta),
                          "z-unpooled" = zunpooled_TX(newDat, Ns),
                          "boschloo" = fisher.2x2(newDat, alternative=alternative),
                          "santner and snell" = santner_TX(newDat, Ns, delta=delta))[3]
          
          # Doesn't appear this check is needed, but useful to confirm if debugging
          if (is.na(newTX)) { stop("Test statistic is NA; please check code") }
          
          newTX <- signif(newTX, 12)  #Remove rounding errors
          
          if (method %in% c("z-pooled", "z-unpooled", "santner and snell")) {
            rejectFlg <- switch(alternative,
                                "less" = (newTX <= TXO),
                                "two.sided" = (abs(newTX) >= abs(TXO)))
          } else if (method == "boschloo") {
            rejectFlg <- (newTX <= TXO)
          }
          
          if (rejectFlg) {
            
            # If more extreme test statistic, then often know the remaining columns in the row is more extreme
            if (i != Ns[1] && j != Ns[2] && j != (Ns[2] - 1) &&
                abs(i/(Ns[1]-i) - j/(Ns[2]-j)) < abs(i/(Ns[1]-i) - (j+1)/(Ns[2]-j-1))) {
              moreExtremeMat[i+1, (j+1):(Ns[2]+1)] <- 1
              break
            } else {
              moreExtremeMat[i+1, j+1] <- 1
            }
            
          } else {
            
            # If less extreme test statistic, then often know the remaining rows in the column is less extreme
            if (i != Ns[1] && i != (Ns[1] - 1) && j != Ns[2] &&
                abs(i/(Ns[1]-i) - j/(Ns[2]-j)) > abs((i+1)/(Ns[1]-i-1) - j/(Ns[2]-j))) {
              moreExtremeMat[(i+1):(Ns[1]+1), j+1] <- 0
            } else {
              moreExtremeMat[i+1, j+1] <- 0
            }
          }
        }
      }
    }
    
    # If delta=0, then two-sided test will be symmetric
    if (alternative == "two.sided") {
      rejectTemp <- which(moreExtremeMat==1, arr.ind = TRUE)
      moreExtremeMat[rep(c(nrow(moreExtremeMat) + 1, ncol(moreExtremeMat) + 1), each=nrow(rejectTemp)) - rejectTemp] <- 1
    }
    
  } else {  #The only case where we can't just calculate the boundary is if we have a two-sided test with delta != 0
    
    TX <- switch(method, 
                 "z-pooled" = zpooled_TX(NULL, Ns, delta=delta),
                 "z-unpooled" = zunpooled_TX(NULL, Ns),
                 "boschloo" = fisher.2x2(NULL, Ns, alternative=alternative),
                 "santner and snell" = santner_TX(NULL, Ns, delta=delta))
    
    # Doesn't appear this check is needed, but useful to confirm if debugging
    if (any(is.na(TX[ , 3]))) { stop("Test statistic is NA; please check code") }
    TX[, 3] <- signif(TX[ , 3], 12)  #Remove rounding errors
    
    TXO <- TX[TX[,1]==data[1,1] & TX[,2]==data[1,2], 3]
    
    # Note: alternative must be "two.sided"
    if (method %in% c("z-pooled", "z-unpooled", "santner and snell")) { rejectFlg <- (abs(TX[,3]) >= abs(TXO))
    } else { rejectFlg <- (TX[,3] <= TXO) }
    
    moreExtremeMat <- matrix(rejectFlg, Ns[1]+1, Ns[2]+1, byrow=TRUE, dimnames=list(0:Ns[1], 0:Ns[2]))*1
  }
  return(list(TXO=TXO, moreExtremeMat=moreExtremeMat))
}
