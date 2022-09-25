multinomCode <-
function(data, alternative, npNumbers, np.interval, beta, method){
  
  N <- sum(data)
  
  # It would be slightly faster to use the convexity property instead of calculating test statistic
  # for all possible tables
  TX <- array(NA, dim=c(N+1, N+1, N+1), dimnames=list(x11=0:N, x12=0:N, x21=0:N))
  for (N1 in 0:N) {
    TXtemp <- switch(method, 
                     "z-pooled" = zpooled_TX(NULL, Ns=c(N1,N-N1), delta=0),
                     "z-unpooled" = zunpooled_TX(NULL, Ns=c(N1,N-N1)),
                     "boschloo" = fisher.2x2(NULL, Ns=c(N1,N-N1), alternative=alternative))
    if (N1 == 0 || N1 == N) {
      if (method=="boschloo") { TXtemp[,3] <- 1
      } else { TXtemp[,3] <- 0 }
    }
    
    for (x11 in 0:N1) {
      x21 <- N1 - x11
      TXtempSubset <- TXtemp[TXtemp[,1] == x11, , drop=FALSE]
      TX[x11+1, TXtempSubset[,2]+1, x21+1] <- TXtempSubset[,3]
    }
  }
  
  TX <- signif(TX, 12)  #Remove rounding errors
  
  #Observed test statistic:
  TXO <- TX[data[1,1]+1, data[1,2]+1, data[2,1]+1]
  
  if (alternative=="less" || method=="boschloo") { moreExtremeArray <- (TX <= TXO)
  } else if (alternative=="greater") { moreExtremeArray <- (TX >= TXO)
  } else if (alternative=="two.sided") { moreExtremeArray <- (abs(TX) >= abs(TXO)) }
  
  #Specify nuisance parameter range
  if (np.interval) {
    tempInt1 <- binom.CI(sum(data[1,]), N, conf.level=1-beta)
    int1 <- seq(max(c(0.00001, tempInt1[1])),min(c(0.99999, tempInt1[2])), length=npNumbers)
    tempInt2 <- binom.CI(sum(data[,1]), N, conf.level=1-beta)
    int2 <- seq(max(c(0.00001, tempInt2[1])),min(c(0.99999, tempInt2[2])), length=npNumbers)
  } else {
    int1 <- seq(0.00001, .99999, length=npNumbers)
    int2 <- seq(0.00001, .99999, length=npNumbers)
    beta <- 0
  }
  
  #Search for the maximum p-value:
  maxP <- maxPvalueMultinom(moreExtremeArray, N, int1, int2, beta)
  
  return(list(method=method, p.value=maxP$pvalue, test.statistic=TXO, np1=maxP$np1, np2=maxP$np2,
              np1.range=c(min(int1),max(int1)), np2.range=c(min(int2), max(int2))))
}
