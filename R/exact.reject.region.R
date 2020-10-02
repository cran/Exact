exact.reject.region <-
function(n1, n2, alternative=c("two.sided", "less", "greater"), alpha=0.05, 
                                npNumbers=100, np.interval=FALSE, beta=0.001,
                                method=c("z-pooled", "z-unpooled", "boschloo", "santner and snell", "csm", "csm approximate", "fisher", "chisq", "yates chisq"),
                                tsmethod=c("square", "central"), ref.pvalue=TRUE, delta=0, convexity=TRUE){
  
  alternative <- match.arg(tolower(alternative), c("two.sided", "less", "greater"))
  tsmethod <- match.arg(tolower(tsmethod), c("square", "central"))
  
  # if two-sided and central, easiest to simply find rejection regions when greater and less and dividing alpha by 2
  if (alternative == "two.sided" && tsmethod == "central") {
    rejectRegionUpper <- exact.reject.region(n1=n1, n2=n2, alternative="less", alpha=alpha/2, 
                                             npNumbers=npNumbers, np.interval=np.interval, beta=beta,
                                             method=method, tsmethod=tsmethod, ref.pvalue=ref.pvalue, delta=delta, convexity=convexity)
    rejectRegionLower <- exact.reject.region(n1=n1, n2=n2, alternative="greater", alpha=alpha/2, 
                                             npNumbers=npNumbers, np.interval=np.interval, beta=beta,
                                             method=method, tsmethod=tsmethod, ref.pvalue=ref.pvalue, delta=delta, convexity=convexity)
    rejectRegion <- rejectRegionUpper + rejectRegionLower
    if (!all(rejectRegion %in% c(0,1))) { stop("Check rejection region logic") }
    return(rejectRegion)
  }
  
  #Perform several checks
  stopifnot(is.logical(np.interval) && is.logical(ref.pvalue))
  if (alpha < 0 || alpha >= 0.5) { stop("To improve code efficiency, alpha must be between 0 and 0.5") }
  if (np.interval && (beta < 0 || beta > 1)) { stop("Beta must be between 0 and 1") }
  if (npNumbers < 1) { stop("Total number of nuisance parameters considered must be at least 1") }
  
  # The Z-pooled statistic that calculates the variance using MLE, which is the pooled variance if delta=0.
  # The Z-pooled statistic is also (perhaps better) known as the Score statistic
  # The classic z-pooled statistic is not performed as the performance is inferior when delta != 0
  if (length(method)==1 && tolower(method)=="score") { method <- "z-pooled" }
  #if (length(method)==1 && tolower(method)=="wald") { method <- "z-unpooled" }
  method <- match.arg(tolower(method), c("z-pooled", "z-unpooled", "boschloo", "santner and snell",
                                         "csm", "csm approximate", "fisher", "chisq", "yates chisq"))
  
  if (n1 <= 0 || n2 <= 0) { stop("Fixed sample sizes must be greater than 0") }
  
  if (method %in% c("csm", "csm approximate", "fisher", "chisq", "yates chisq") && np.interval) {
    warning("Interval of nuisance parameter cannot be used with CSM, fisher, or chi-square test; np.interval changed to FALSE")
    np.interval <- FALSE
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
  
  if (!np.interval && !(method %in% c("csm", "csm approximate", "fisher"))) {
    int <- seq(0.00001, .99999, length=npNumbers) 
    lookupArray <- dbinomCalc(Ns = c(n1, n2), int, delta)
    
    # For the tests with a test statistic:
    TX <- switch(method,
                 "z-pooled" = zpooled_TX(NULL, c(n1, n2), delta=delta),
                 "z-unpooled" = zunpooled_TX(NULL, c(n1, n2), delta=delta),
                 "boschloo" = fisher.2x2(NULL, c(n1, n2), alternative=alternative),
                 "santner and snell" = santner_TX(NULL, c(n1, n2), delta=delta),
                 "csm approximate" = csmApprox_TX(NULL, c(n1, n2), alternative, int, lookupArray),
                 "chisq" = chisq_TX(NULL, c(n1, n2), yates=FALSE),
                 "yates chisq" = chisq_TX(NULL, c(n1, n2), yates=TRUE))
    
    # Instead of starting with the most extreme table and adding subquentially, start in the middle and move up or down
    # depending on whether we have a more extreme table.  This can greatly reduce the number of p-values calculated
    TX[is.na(TX[,3]), 3] <- 0
    
    # If two-sided test, then first table will be the largest absolute value of the test statistic.
    # Take the negative value to match the "less" alternative
    if (alternative == "two.sided" && !(method %in% c("boschloo", "csm approximate")) || (method %in% c("chisq", "yates chisq"))) {TX[, 3] <- -abs(TX[,3])}
    if (method %in% c("chisq", "yates chisq") && alternative == "less") {TX[TX[,1]/n1 >= TX[,2]/n2, 3] <- -TX[TX[,1]/n1 >= TX[,2]/n2, 3]}
    TX <- signif(TX[order(TX[,3]), ], 12)
    TX <- cbind(TX, NA)
    # We can just cross out all test statistics >=0 since this p-value would be >= 0.5.
    if (method %in% c("z-pooled", "z-unpooled", "santner and snell", "chisq", "yates chisq") & delta==0) { TX[TX[,3] >= 0, 4] <- FALSE }
    
    moreExtremeTbls <- searchExtreme(TX = TX, n1 = n1, n2 = n2, alternative = alternative, method = method, int = int, delta = delta,
                                     alpha = alpha, lookupArray = lookupArray)
  } else if (method == "csm") {
    int <- seq(0.00001, .99999, length=npNumbers)
    rejectRegion <- moreExtremeCSM(data=NULL, Ns = c(n1, n2), alternative = alternative,
                                   int = int, doublePvalue = FALSE, delta = delta, reject.alpha = alpha)$moreExtremeMat
    if (swapFlg){ rejectRegion <- t(rejectRegion) }
    return(rejectRegion)
  } else if (method %in% c("csm approximate", "fisher") || np.interval) {
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
          tables <- matrix(c(i, n1-i, j, n2-j), 2, 2, byrow=TRUE)
          # Since we changed p1 and p2 if "greater", always set the test to "less"
          if (method=="fisher") {
            rejectRegionTemp <- (fisher.2x2(tables, alternative=alternative)[3] <= alpha)
          } else {
            rejectRegionTemp <- (exact.test(tables, npNumbers=npNumbers, alternative=alternative, np.interval=np.interval,
                                            beta=beta, method=method, tsmethod=tsmethod, to.plot=FALSE, ref.pvalue=ref.pvalue, delta=delta, reject.alpha=alpha))
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
}
