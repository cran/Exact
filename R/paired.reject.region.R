paired.reject.region <-
function(N, alternative=c("two.sided", "less", "greater"), alpha=0.05, 
                                 npNumbers=100, np.interval=FALSE, beta=0.001,
                                 method=c("uam", "ucm", "uamcc", "csm", "cm", "am", "amcc"),
                                 tsmethod=c("square", "central"), delta=0, convexity=TRUE, useStoredCSM=TRUE){
  
  alternative <- match.arg(tolower(alternative), c("two.sided", "less", "greater"))
  method <- convertMethod(method)
  tsmethod <- match.arg(tolower(tsmethod), c("square", "central"))
  
  #Perform several checks
  checkPairedParam(N=N, alternative=alternative, alpha=alpha, npNumbers=npNumbers, np.interval=np.interval, beta=beta,
                   method=method, tsmethod=tsmethod, delta=delta, convexity=convexity, useStoredCSM=useStoredCSM)
  
  # if two-sided and central, easiest to simply find rejection regions when greater and less and dividing alpha by 2
  if (alternative == "two.sided" && tsmethod == "central") {
    rejectRegionUpper <- paired.reject.region(N=N, alternative="less", alpha=alpha/2, 
                                             npNumbers=npNumbers, np.interval=np.interval, beta=beta,
                                             method=method, tsmethod=tsmethod, delta=delta, convexity=convexity,
                                             useStoredCSM=useStoredCSM)
    rejectRegionLower <- paired.reject.region(N=N, alternative="greater", alpha=alpha/2, 
                                             npNumbers=npNumbers, np.interval=np.interval, beta=beta,
                                             method=method, tsmethod=tsmethod, delta=delta, convexity=convexity,
                                             useStoredCSM=useStoredCSM)
    rejectRegion <- rejectRegionUpper + rejectRegionLower
    if (!all(rejectRegion %in% c(0,1,NA))) { stop("Check rejection region logic") }
    return(rejectRegion)
  }
  
  # Instead of trying to solve both "less" and "greater", just have code only work for "less"
  # Will transpose at the end
  alternativeTemp <- alternative
  
  swapFlg <- (alternative == "greater" || (alternative=="two.sided" && delta < 0))
  if (swapFlg) {
    delta <- -delta
    if (alternative == "greater") { alternative <- "less" }
  } 
  
  # Each test has a different way of minimizing the computation time to calculate the rejection region.
  # Most tests have a test statistic which shows the order of the tables considered.  This can be utilized to greatly speed up time.
  # The CSM test is essentially just calculating the p-value except stopping at alpha
  # The conditional mcnemar test doesn't have an ordered test statistic, so calculate sequentially with convexity property
  # The np.interval approach is a special case since the NP interval depends on the number of successes, which makes it not completely
  # convex.  The precise way is to calculate a p-value for all tables, but in almost all cases the convexity property holds
  
  if (!np.interval && !(method %in% c("csm", "cm"))) {
    
    if (delta == 0) { int <- seq(0.00001,.49999,length=npNumbers)
    } else if (delta > 0) { int <- seq(0.00001, (1-delta)/2 - 0.00001, length=npNumbers)
    } else if (delta < 0) { int <- seq(-delta + 0.00001, (1 - delta)/2 - 0.00001, length=npNumbers)}
    
    lookupArray <- trinomCalc(N, int, delta)
    
    # For the tests with a test statistic:
    TX <- switch(method,
                 "uam" = mcnemar_TX(NULL, N, delta=delta, CC=FALSE),
                 "ucm" = mcnemar.2x2(NULL, N, alternative=alternative, pval1ties = TRUE),
                 "uamcc" = mcnemar_TX(NULL, N, delta=delta, CC=TRUE),
                 "am" = mcnemar_TX(NULL, N, delta=delta, CC=FALSE),
                 "amcc" = mcnemar_TX(NULL, N, delta=delta, CC=TRUE))
    
    # Instead of starting with the most extreme table and adding sequentially, start in the middle and move up or down
    # depending on whether we have a more extreme table.  This can greatly reduce the number of p-values calculated
    
    # If two-sided test, then first table will be the largest absolute value of the test statistic.
    # Take the negative value to match the "less" alternative
    if (alternative == "two.sided" && method != "ucm") {TX[, 3] <- -abs(TX[,3])}
    
    TX <- signif(TX[order(TX[,3]), ], 12)
    TX <- cbind(TX, NA)
    
    # We can just cross out all test statistics >=0 since this p-value would be >= 0.5.
    if (method %in% c("uam", "uamcc", "am", "amcc") && delta==0) { TX[TX[,3] >= 0, 4] <- FALSE }
    
    # For ucm, we know significant if Conditional McNemar's p-value < alpha
    if (method == "ucm" && delta==0) { TX[TX[,3] <= alpha, 4] <- TRUE }
    
    moreExtremeTbls <- searchExtremePaired(TX = TX, N = N, alternative = alternative, method = method, int = int, delta = delta,
                                           alpha = alpha, lookupArray = lookupArray)
    
    rejectRegion <- matrix(0, N+1, N+1)
    rejectRegion <- lowerMatVal(rejectRegion, NA)
    if (nrow(moreExtremeTbls) > 0) { #If at least one more extreme table
      for (i in 1:nrow(moreExtremeTbls)) {
        rejectRegion[moreExtremeTbls[i, 1] + 1, moreExtremeTbls[i, 2] + 1]  <- 1
      }
    }
    rownames(rejectRegion) <- 0:N
    colnames(rejectRegion) <- 0:N
    if (swapFlg){ rejectRegion <- t(rejectRegion) }
    
    return(rejectRegion)
    
  } else if (method == "csm") {
    
    # Use stored ordering matrix if available
    if ( N <= 200 && delta == 0 && useStoredCSM) {
      
      if (!requireNamespace("ExactData", quietly = TRUE)) {
        stop(paste("ExactData R package must be installed when useStoredCSM=TRUE. To install ExactData R package, run:",
                   "`install.packages('ExactData', repos='https://pcalhoun1.github.io/drat/', type='source')`"))
      }
      
      # Note: delta=0
      int <- seq(0.00001, .49999, length=npNumbers) 
      lookupArray <- trinomCalc(N, int, delta)
      
      if (alternative == "two.sided") { orderMat <- ExactData::orderCSMPairedMatTwoSided[[N]]
      } else { orderMat <- ExactData::orderCSMPairedMatOneSided[[N]] }
      
      orderVec <- as.vector(orderMat)[!is.na(orderMat)]
      TX <- matrix(cbind(which(!is.na(orderMat), arr.ind = TRUE) - 1, orderVec), ncol=3, dimnames = NULL)
      TX <- signif(TX[order(TX[,3]), ], 12)
      TX <- cbind(TX, NA)

      moreExtremeTbls <- searchExtremePaired(TX = TX, N = N, alternative = alternative, method = method, int = int, delta = delta,
                                             alpha = alpha, lookupArray = lookupArray)
      
      rejectRegion <- matrix(0, N+1, N+1)
      rejectRegion <- lowerMatVal(rejectRegion, NA)
      if (nrow(moreExtremeTbls) > 0) { #If at least one more extreme table
        for (i in 1:nrow(moreExtremeTbls)) {
          rejectRegion[moreExtremeTbls[i, 1] + 1, moreExtremeTbls[i, 2] + 1]  <- 1
        }
      }
      rownames(rejectRegion) <- 0:N
      colnames(rejectRegion) <- 0:N
      
    } else {
    
      if (delta == 0) { int <- seq(0.00001,.49999,length=npNumbers)
      } else if (delta > 0) { int <- seq(0.00001, (1-delta)/2 - 0.00001, length=npNumbers)
      } else if (delta < 0) { int <- seq(-delta + 0.00001, (1 - delta)/2 - 0.00001, length=npNumbers)}
      
      rejectRegion <- moreExtremePairedCSM(data = NULL, N = N, alternative = alternative,
                                           int = int, doublePvalue = FALSE, delta = delta, reject.alpha = alpha)$moreExtremeMat
      
    }
    
    if (swapFlg) { rejectRegion <- t(rejectRegion) }
    return(rejectRegion)
    
  } else if (method == "cm" || np.interval) {
    
    rejectRegion <- matrix(NA, N+1, N+1)
    # Lower diagonal always starts at 0
    rejectRegion[round((row(rejectRegion)-1)/(nrow(rejectRegion)-1), digits=10) >= round((col(rejectRegion)-1)/(ncol(rejectRegion)-1), digits=10)] <- 0
    rejectRegion <- lowerMatVal(rejectRegion)
    
    # Conditional McNemar's test has the convexity property, so can determine the entire rejection region through the boundary.
    # Essentially, start at first row and move through each column until table rejects.  Then move to next row
    
    # For a given alpha and sample size, the np.interval approach almost always has the convexity property.  However, there are
    # very rare instances when the convexity property doesn't hold and one has to consider all possible tables.  Add an input parameter
    # for whether user is wanting to assume convexity (very fast), or do it correctly by going through all possible tables (very slow)
    
    for (i in 0:N) {  #Go through each row
      # Find first column with NA
      startJ <- which(is.na(rejectRegion[i+1, ]))[1] - 1
      if (!is.na(startJ)) {
        for (j in startJ:(N-i)) {
          tables <- matrix(c(N-i-j, i, j, 0), 2, 2, byrow=TRUE)
          if (method=="cm") {
            rejectRegionTemp <- (mcnemar.2x2(tables, N, alternative=alternative, pval1ties = FALSE)[,3] <= alpha)
          } else {
            rejectRegionTemp <- pairedCode(tables, alternative=alternative, npNumbers=npNumbers,
                                           np.interval=np.interval, beta=beta, method=method, tsmethod=tsmethod, to.plot=FALSE,
                                           ref.pvalue=FALSE, delta=delta, reject.alpha=alpha, useStoredCSM=useStoredCSM)
          }
          
          if (rejectRegionTemp) {
            if (method=="cm" || convexity) {
              # If p-value <= alpha, then know the remaining columns in the row is more extreme
              rejectRegion[i+1, (j+1):(N+1)] <- 1
              if (alternative=="two.sided") { rejectRegion[(j+1):(N+1), i+1] <- 1 }
              break
            } else {
              rejectRegion[i+1, j+1] <- 1
              if (alternative=="two.sided") { rejectRegion[j+1, i+1] <- 1 }
            }
          } else {
            if (method=="cm" || convexity) {
              # If p-value > alpha, then know the remaining rows in the column is less extreme
              rejectRegion[is.na(rejectRegion[ , j+1]), j+1] <- 0
            } else { rejectRegion[i+1, j+1] <- 0 }
          }
        }
      }
    }
    rownames(rejectRegion) <- 0:N
    colnames(rejectRegion) <- 0:N
    rejectRegion <- lowerMatVal(rejectRegion, NA)
    if (swapFlg){ rejectRegion <- t(rejectRegion) }
    return(rejectRegion)
  }
}
