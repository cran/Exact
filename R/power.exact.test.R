power.exact.test <-
function(p1, p2, n1, n2, alternative=c("two.sided", "less", "greater"), alpha=0.05,
                             npNumbers=100, np.interval=FALSE, beta=0.001,
                             method=c("z-pooled", "z-unpooled", "boschloo", "santner and snell", "csm", "fisher", "pearson chisq", "yates chisq"),
                             tsmethod=c("square", "central"), simulation=FALSE, nsim = 100, delta=0, convexity=TRUE, useStoredCSM=TRUE){
  
  alternative <- match.arg(tolower(alternative), c("two.sided", "less", "greater"))
  # The Z-pooled statistic is actually the Score statistic, which are equivalent when delta = 0
  # The classic Z-pooled statistic is not performed as the performance is inferior when delta != 0
  if (length(method)==1 && tolower(method)=="score") { method <- "z-pooled" }
  method <- match.arg(tolower(method), c("z-pooled", "z-unpooled", "boschloo", "santner and snell", "csm",
                                         "fisher", "pearson chisq", "yates chisq"))
  tsmethod <- match.arg(tolower(tsmethod), c("square", "central"))
  
  #Perform several checks
  checkParam(p1=p1, p2=p2, n1=n1, n2=n2, alternative=alternative, alpha=alpha, npNumbers=npNumbers, np.interval=np.interval, beta=beta,
             method=method, tsmethod=tsmethod, simulation=simulation, nsim = nsim, delta=delta, convexity=convexity, useStoredCSM=useStoredCSM)
  
  # Consider all tables:
  if (!simulation) {
    rejectRegion <- exact.reject.region(n1 = n1, n2 = n2, alternative=alternative, alpha=alpha,
                                        npNumbers=npNumbers, np.interval=np.interval, beta=beta,
                                        method=method, tsmethod=tsmethod, delta=delta, convexity=convexity,
                                        useStoredCSM=useStoredCSM)
    # Note: don't use delta for probabilities
    prob <- dbinom(0:n1, n1, p1) %*% t(dbinom(0:n2, n2, p2))
    power <- sum(prob[as.logical(rejectRegion)])
    
  } else if (simulation) {
    #Randomly generate a table based on known proportions
    randA <- rbinom(nsim, size = n1, prob = p1)
    randC <- rbinom(nsim, size = n2, prob = p2)
    randTables <- cbind(randA, n1 - randA, randC, n2 - randC)
    moreExtremeSim <- rep(NA, nsim)
    for (i in 1:nsim) {
      if (delta != 0 || (alternative=="greater" && (randTables[i,1]/n1 - randTables[i,3]/n2) > 0) ||
          (alternative=="less" && (randTables[i,1]/n1 - randTables[i,3]/n2) < 0) || 
          (alternative=="two.sided" && randTables[i,1]/n1 != randTables[i,3]/n2)) {
        if (method=="fisher") {
          moreExtremeSim[i] <- (fisher.2x2(matrix(randTables[i,], ncol=2, nrow=2), alternative=alternative)[3] <= alpha)
        } else if (method %in% c("pearson chisq", "yates chisq")) {
          moreExtremeSim[i] <- suppressWarnings(prop.test(matrix(randTables[i,], nrow=2, ncol=2), alternative=alternative,
                                                       correct=(method == "yates chisq"))$p.value <= alpha)
        } else {
          moreExtremeSim[i] <- (binomialCode(matrix(randTables[i,], ncol=2, nrow=2), alternative=alternative, npNumbers=npNumbers,
                                          np.interval=np.interval, beta=beta, method=method, tsmethod=tsmethod, to.plot=FALSE,
                                          ref.pvalue=FALSE, delta=delta, reject.alpha=alpha, useStoredCSM=useStoredCSM))
        }
      } else { moreExtremeSim[i] <- FALSE }
    }
    power <- mean(moreExtremeSim)
  }
  
  #Convert data to power.htest structure
  methodDescribed <- methodText(method, np.interval)
  
  return(structure(list("n1, n2" = c(n1, n2), "p1, p2" = c(p1, p2), alpha = alpha, 
                        power = power, alternative = alternative, delta = delta, method=methodDescribed), class = "power.htest"))
}
