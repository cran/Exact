moreExtreme <-
function(method, data, Ns, alternative, int, delta){
  
  # If alternative="less", then only need to calculate boundary in upper triangle (faster than calculating all test statistics).
  # If alternative="two.sided" and delta=0, then can also calculate boundary in upper triangle since the lower triangle
  # will always be symmetric.  The only case where we can't just calculate the boundary but instead must calculate the test statistic
  # for all tables is when you have a two-sided test with delta != 0
  if (alternative != "two.sided" || delta == 0) {
    
    moreExtremeMat <- matrix(NA, Ns[1]+1, Ns[2]+1, dimnames=list(0:Ns[1], 0:Ns[2]))
    TXO <- switch(method, 
                  "z-pooled" = zpooled_TX(data, Ns, delta=delta),
                  "z-unpooled" = zunpooled_TX(data, Ns),
                  "boschloo" = fisher.2x2(data, NULL, alternative=alternative),
                  "santner and snell" = santner_TX(data, Ns),
                  "csm approximate" = csmApprox_TX(data, Ns, alternative, int))
    
    if (is.na(TXO)) { TXO <- 0 }
    TXO <- signif(TXO, 12)  #Remove rounding errors
    
    #for (j in 0:(data[1,2]+1)) { moreExtremeMat[(data[1,1]+1):(Ns[1]+1), j] <- 0 }
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
                          "boschloo" = fisher.2x2(newDat, Ns = NULL, alternative=alternative),
                          "santner and snell" = santner_TX(newDat, Ns),
                          "csm approximate" = csmApprox_TX(newDat, Ns, alternative, int))
          if (is.na(newTX)) { newTX <- 0 }
          newTX <- signif(newTX, 12)  #Remove rounding errors
          
          if (method %in% c("z-pooled", "z-unpooled", "santner and snell")) {
            rejectFlg <- switch(alternative,
                                "less" = (newTX <= TXO),
                                "two.sided" = (abs(newTX) >= abs(TXO)))
          } else if (method %in% c("boschloo", "csm approximate")) {
            rejectFlg <- (newTX <= TXO)
          }
          
          if (rejectFlg) {
            # If more extreme test statistic, then know the remaining columns in the row is more extreme
            moreExtremeMat[i+1, (j+1):(Ns[2]+1)] <- 1
            break
          } else {
            # If less extreme test statistic, then know the remaining rows in the column is less extreme
            moreExtremeMat[(i+1):(Ns[1]+1), j+1] <- 0
          }
        }
      }
    }
    
    # If delta=0, then two-sided test will be symmetric
    if (alternative=="two.sided") {
      rejectTemp <- which(moreExtremeMat==1, arr.ind = TRUE)
      moreExtremeMat[rep(c(nrow(moreExtremeMat) + 1, ncol(moreExtremeMat) + 1), each=nrow(rejectTemp)) - rejectTemp] <- 1
    }
    
  } else {
    
    if (method != "z-pooled") { stop("Delta != 0 and not z-pooled?") }
    TX <- zpooled_TX(NULL, Ns, delta)
    TX[is.na(TX[, 3]), 3] <- 0
    TX[, 3] <- signif(TX[ , 3], 12)  #Remove rounding errors
    TXO <- TX[TX[,1]==data[1,1] & TX[,2]==data[1,2], 3]
    moreExtremeMat <- matrix(abs(TX[,3]) >= abs(TXO), Ns[1]+1, Ns[2]+1, byrow=TRUE, dimnames=list(0:Ns[1], 0:Ns[2]))*1
  }
  return(list(TXO=TXO, moreExtremeMat=moreExtremeMat))
}
