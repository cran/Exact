exact.reject.region <-
function(n1, n2, alternative=c("two.sided", "less", "greater"), alpha=0.05, 
                                npNumbers=100, np.interval=FALSE, beta=0.001,
                                method=c("z-pooled", "z-unpooled", "boschloo", "santner and snell", "csm", "fisher", "pearson chisq", "yates chisq"),
                                tsmethod=c("square", "central"), delta=0, convexity=TRUE, useStoredCSM=TRUE){

  alternative <- match.arg(tolower(alternative), c("two.sided", "less", "greater"))
  # The Z-pooled statistic is actually the Score statistic, which are equivalent when delta = 0
  # The classic Z-pooled statistic is not performed as the performance is inferior when delta != 0
  if (length(method)==1 && tolower(method)=="score") { method <- "z-pooled" }
  method <- match.arg(tolower(method), c("z-pooled", "z-unpooled", "boschloo", "santner and snell", "csm",
                                         "fisher", "pearson chisq", "yates chisq"))
  tsmethod <- match.arg(tolower(tsmethod), c("square", "central"))

  #Perform several checks
  checkParam(n1=n1, n2=n2, alternative=alternative, alpha=alpha, npNumbers=npNumbers, np.interval=np.interval, beta=beta,
             method=method, tsmethod=tsmethod, delta=delta, convexity=convexity, useStoredCSM=useStoredCSM)
    
  # if two-sided and central, easiest to simply find rejection regions when greater and less and dividing alpha by 2
  if (alternative == "two.sided" && tsmethod == "central") {
    rejectRegionUpper <- exact.reject.region(n1=n1, n2=n2, alternative="less", alpha=alpha/2, 
                                             npNumbers=npNumbers, np.interval=np.interval, beta=beta,
                                             method=method, tsmethod=tsmethod, delta=delta, convexity=convexity,
                                             useStoredCSM=useStoredCSM)
    rejectRegionLower <- exact.reject.region(n1=n1, n2=n2, alternative="greater", alpha=alpha/2, 
                                             npNumbers=npNumbers, np.interval=np.interval, beta=beta,
                                             method=method, tsmethod=tsmethod, delta=delta, convexity=convexity,
                                             useStoredCSM=useStoredCSM)
    rejectRegion <- rejectRegionUpper + rejectRegionLower
    if (!all(rejectRegion %in% c(0,1))) { stop("Check rejection region logic") }
    return(rejectRegion)
  }
    
  # Instead of trying to solve both "less" and "greater", just switch n1 and n2 if "greater" and have code only work for "less"
  # Will transpose at the end
  alternativeTemp <- alternative
  
  swapFlg <- (alternative == "greater" || (alternative=="two.sided" && delta < 0))
  if (swapFlg) {
    n2temp <- n2
    n2 <- n1; n1 <- n2temp;
    delta <- -delta
    if (alternative == "greater") { alternative <- "less" }
  } 
  
  # Each test has a different way of minimizing the computation time to calculate the rejection region.
  # Most tests have a test statistic which shows the order of the tables considered.  This can be utilized to greatly speed up time.
  # The CSM test is essentially just calculating the p-value except stopping at alpha
  # The Fisher test is using the fisher.2x2 test, but it's faster to use the fisher.2x2 test sequentially (instead of all 2x2 tables)
  # The np.interval approach is a special case since the NP interval depends on the number of successes, which makes it not completely
  # convex.  The precise way is to calculate a p-value for all tables, but in almost all cases the convexity property holds
  
  if (!np.interval && !(method %in% c("csm", "fisher"))) {
    
    if (delta == 0) { int <- seq(0.00001,.99999,length=npNumbers)
    } else if (delta > 0) { int <- seq(0.00001, 1 - delta - 0.00001, length=npNumbers)
    } else if (delta < 0) { int <- seq(abs(delta) + 0.00001, .99999, length=npNumbers)}
    
    lookupArray <- dbinomCalc(Ns = c(n1, n2), int, delta)
    
    # For the tests with a test statistic:
    TX <- switch(method,
                 "z-pooled" = zpooled_TX(NULL, c(n1, n2), delta=delta),
                 "z-unpooled" = zunpooled_TX(NULL, c(n1, n2), delta=delta),
                 "boschloo" = fisher.2x2(NULL, c(n1, n2), alternative=alternative),
                 "santner and snell" = santner_TX(NULL, c(n1, n2), delta=delta),
                 "pearson chisq" = chisq_TX(NULL, c(n1, n2), yates=FALSE),
                 "yates chisq" = chisq_TX(NULL, c(n1, n2), yates=TRUE))
    
    # Instead of starting with the most extreme table and adding sequentially, perform binary search by starting in the middle and either
    # rejecting or failing to reject more extreme tables based on p-value.  This can greatly reduce the number of p-values calculated
    TX[is.na(TX[,3]), 3] <- 0
    
    # If two-sided test, then first table will be the largest absolute value of the test statistic.
    # Take the negative value to match the "less" alternative
    if ((alternative == "two.sided" && method != "boschloo") || (method %in% c("pearson chisq", "yates chisq"))) {TX[, 3] <- -abs(TX[,3])}
    if (method %in% c("pearson chisq", "yates chisq") && alternative == "less") {TX[TX[,1]/n1 >= TX[,2]/n2, 3] <- -TX[TX[,1]/n1 >= TX[,2]/n2, 3]}
    
    TX <- signif(TX[order(TX[,3]), ], 12)
    TX <- cbind(TX, NA)
    
    # We can just cross out all test statistics >=0 since this p-value would be >= 0.5.
    if (method %in% c("z-pooled", "z-unpooled", "santner and snell", "pearson chisq", "yates chisq") && delta==0) { TX[TX[,3] >= 0, 4] <- FALSE }
    
    # For boschloo, we know significant if Fisher's p-value < alpha
    if (method == "boschloo" && delta==0) { TX[TX[,3] <= alpha, 4] <- TRUE }
    
    moreExtremeTbls <- searchExtreme(TX = TX, n1 = n1, n2 = n2, alternative = alternative, method = method, int = int, delta = delta,
                                     alpha = alpha, lookupArray = lookupArray)
    
    rejectRegion <- matrix(0, n1+1, n2+1)
    if (nrow(moreExtremeTbls) > 0) { #If at least one more extreme table
      for (i in 1:nrow(moreExtremeTbls)) {
        rejectRegion[moreExtremeTbls[i, 1] + 1, moreExtremeTbls[i, 2] + 1]  <- 1
      }
    }
    rownames(rejectRegion) <- 0:n1
    colnames(rejectRegion) <- 0:n2
    if (swapFlg){ rejectRegion <- t(rejectRegion) }
    
    return(rejectRegion)
    
  } else if (method == "csm") {
    
    # Use stored ordering matrix if available
    if (max(c(n1, n2)) <= 100 && delta == 0 && useStoredCSM) {
      
      if (!requireNamespace("ExactData", quietly = TRUE)) {
        stop(paste("ExactData R package must be installed when useStoredCSM=TRUE. To install ExactData R package, run:",
                   "`install.packages('ExactData', repos='https://pcalhoun1.github.io/drat/', type='source')`"))
      }
      
      # Note: delta=0
      int <- seq(0.00001, .99999, length=npNumbers) 
      lookupArray <- dbinomCalc(Ns = c(n1, n2), int, delta)
      
      if (alternative == "two.sided") {
        if (n1 < n2) {  # Must convert known matrix of (n2,n1) to (n1,n2)
          orderMat <- t(ExactData::orderCSMMatTwoSided[[paste0("(",n2,",",n1,")")]])
          orderMat <- orderMat[(n1+1):1, (n2+1):1]
          rownames(orderMat) <- 0:n1
          colnames(orderMat) <- 0:n2
        } else { orderMat <- ExactData::orderCSMMatTwoSided[[paste0("(",n1,",",n2,")")]] }
        
      } else {
        if (n1 < n2) {  # Must convert known matrix of (n2,n1) to (n1,n2)
          orderMat <- t(ExactData::orderCSMMatOneSided[[paste0("(",n2,",",n1,")")]])
          orderMat <- orderMat[(n1+1):1, (n2+1):1]
          rownames(orderMat) <- 0:n1
          colnames(orderMat) <- 0:n2
        } else { orderMat <- ExactData::orderCSMMatOneSided[[paste0("(",n1,",",n2,")")]] }
      }
      
      TX <- matrix(cbind(as.matrix(expand.grid(0:n1, 0:n2)), as.vector(orderMat)), ncol=3, dimnames = NULL)
      TX <- signif(TX[order(TX[,3]), ], 12)
      TX <- cbind(TX, NA)

      moreExtremeTbls <- searchExtreme(TX = TX, n1 = n1, n2 = n2, alternative = alternative, method = method, int = int, delta = delta,
                                       alpha = alpha, lookupArray = lookupArray)
      
      rejectRegion <- matrix(0, n1+1, n2+1)
      if (nrow(moreExtremeTbls) > 0) { #If at least one more extreme table
        for (i in 1:nrow(moreExtremeTbls)) {
          rejectRegion[moreExtremeTbls[i, 1] + 1, moreExtremeTbls[i, 2] + 1]  <- 1
        }
      }
      rownames(rejectRegion) <- 0:n1
      colnames(rejectRegion) <- 0:n2
      
    } else {

      if (delta == 0) { int <- seq(0.00001,.99999,length=npNumbers)
      } else if (delta > 0) { int <- seq(0.00001, 1 - delta - 0.00001, length=npNumbers)
      } else if (delta < 0) { int <- seq(abs(delta) + 0.00001, .99999, length=npNumbers)}
      
      rejectRegion <- moreExtremeCSM(data = NULL, Ns = c(n1, n2), alternative = alternative,
                                     int = int, doublePvalue = FALSE, delta = delta, reject.alpha = alpha)$moreExtremeMat
      
    }
    
    if (swapFlg) { rejectRegion <- t(rejectRegion) }
    return(rejectRegion)
    
  } else if (method == "fisher" || np.interval) {
    
    rejectRegion <- matrix(NA, n1+1, n2+1)
    # Lower diagonal always starts at 0
    rejectRegion[round((row(rejectRegion)-1)/(nrow(rejectRegion)-1), digits=10) >= round((col(rejectRegion)-1)/(ncol(rejectRegion)-1), digits=10)] <- 0
    
    # Fisher's test has the convexity property, so can determine the entire rejection region through the boundary.
    # Essentially, start at first row and move through each column until table rejects.  Then move to next row
    
    # For a given alpha and sample size, the np.interval approach almost always has the convexity property.  However, there are
    # very rare instances when the convexity property doesn't hold and one has to consider all possible tables.  Add an input parameter
    # for whether user is wanting to assume convexity (very fast), or do it correctly by going through all possible tables (very slow)
    
    for (i in 0:n1) {  #Go through each row
      # Find first column with NA
      startJ <- which(is.na(rejectRegion[i+1, ]))[1] - 1
      if (!is.na(startJ)) {
        for (j in startJ:n2) {
          tables <- matrix(c(i, n1-i, j, n2-j), ncol=2, nrow=2)
          if (method=="fisher") {
            rejectRegionTemp <- (fisher.2x2(tables, alternative=alternative)[3] <= alpha)
          } else {
            rejectRegionTemp <- (binomialCode(tables, alternative=alternative, npNumbers=npNumbers, np.interval=np.interval,
                                              beta=beta, method=method, tsmethod=tsmethod, to.plot=FALSE,
                                              ref.pvalue=FALSE, delta=delta, reject.alpha=alpha, useStoredCSM=useStoredCSM))
          }
          
          if (rejectRegionTemp) {
            if (method=="fisher" || convexity) {
              # If p-value <= alpha, then know the remaining columns in the row is more extreme
              rejectRegion[i+1, (j+1):(n2+1)] <- 1
              if (alternative=="two.sided") { rejectRegion[n1-i+1, 1:(n2-j+1)] <- 1 }
              break
            } else {
              rejectRegion[i+1, j+1] <- 1
              if (alternative=="two.sided") { rejectRegion[n1-i+1, n2-j+1] <- 1 }
            }
          } else {
            if (method=="fisher" || convexity) {
              # If p-value > alpha, then know the remaining rows in the column is less extreme
              rejectRegion[is.na(rejectRegion[ , j+1]), j+1] <- 0
            } else { rejectRegion[i+1, j+1] <- 0 }
          }
        }
      }
    }
    rownames(rejectRegion) <- 0:n1
    colnames(rejectRegion) <- 0:n2
    if (swapFlg){ rejectRegion <- t(rejectRegion) }
    return(rejectRegion)
  }
}
