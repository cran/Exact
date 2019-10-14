power.exact.test <-
function(p1, p2, n1, n2, alternative=c("two.sided", "less", "greater"), alpha=0.05,
                             npNumbers=100, np.interval=FALSE, beta=0.001,
                             method=c("z-pooled", "z-unpooled", "boschloo", "santner and snell", "csm", "csm approximate", "fisher", "chisq", "yates chisq"),
                             ref.pvalue=TRUE, simulation=FALSE, nsim = 100, delta=0, convexity=TRUE){
  
  stopifnot(is.logical(np.interval) && is.logical(ref.pvalue) && is.logical(simulation) && is.logical(convexity))
  if (alpha < 0 || alpha >= 0.5) { stop("To improve code efficiency, alpha must be between 0 and 0.5") }
  if (np.interval && (beta < 0 || beta > 1)) { stop("Beta must be between 0 and 1") }
  if (npNumbers < 1) { stop("Total number of nuisance parameters considered must be at least 1") }
  if (nsim < 1) { stop("Need at least one simulation") }
  
  alternative <- match.arg(tolower(alternative), c("two.sided", "less", "greater"))
  
  #Sometimes Z-pooled is called score and Z-unpooled is called wald statistic
  if (length(method)==1 && tolower(method)=="score") { method <- "z-pooled" }
  if (length(method)==1 && tolower(method)=="wald") { method <- "z-unpooled" }
  method <- match.arg(tolower(method), c("z-pooled", "z-unpooled", "boschloo", "santner and snell",
                                         "csm", "csm approximate", "fisher", "chisq", "yates chisq"))
  
  if (p1 < 0 || p1 > 1 || p2 < 0 || p2 > 1) { stop("Probabilities must be between 0 and 1") }
  if (n1 <= 0 || n2 <= 0) { stop("Fixed sample sizes must be greater than 0") }
  
  if (delta != 0 && !(method %in% c("z-pooled", "csm"))) {
    stop("Delta != 0 only works for Z-pooled and CSM tests")
  }
  
  if (method %in% c("csm", "csm approximate", "fisher", "chisq", "yates chisq") && np.interval) {
    warning("Interval of nuisance parameter cannot be used with CSM, fisher, or chi-square test; np.interval changed to FALSE")
    np.interval <- FALSE
  }
  
  # Consider all tables:
  if (!simulation) {
    rejectRegion <- exact.reject.region(n1 = n1, n2 = n2, alternative=alternative, alpha=alpha, npNumbers=npNumbers, np.interval=np.interval, beta=beta,
                                        method=method, ref.pvalue=ref.pvalue, delta=delta, convexity=convexity)
    prob <- dbinom(0:n1, n1, p1) %*% t(dbinom(0:n2, n2, p2))
    power <- sum(prob[as.logical(rejectRegion)])
    
  } else if (simulation) {
    #Randomly generate a table based on known proportions
    randA <- rbinom(nsim, size = n1, prob = p1)
    randC <- rbinom(nsim, size = n2, prob = p2)
    randTables <- cbind(randA, n1 - randA, randC, n2 - randC)
    moreExtreme <- rep(NA, nsim)
    for (i in 1:nsim) {
      if ((alternative=="greater" && (randTables[i,1]/n1 - randTables[i,3]/n2) > delta) ||
          (alternative=="less" && (randTables[i,1]/n1 - randTables[i,3]/n2) < delta) || 
          (alternative=="two.sided" && delta==0 && randTables[i,1]/n1 != randTables[i,3]/n2)) {
        if (method=="fisher") {
          moreExtreme[i] <- (fisher.2x2(matrix(randTables[i,], 2, 2, byrow=TRUE), alternative=alternative) <= alpha)
        } else if (method %in% c("chisq", "yates chisq")) {
          moreExtreme[i] <- suppressWarnings(prop.test(matrix(randTables[i,], 2, 2, byrow=TRUE),
                                                       alternative=alternative, correct=(method == "yates chisq"))$p.value <= alpha)
        } else {
          moreExtreme[i] <- (exact.test(matrix(randTables[i,], 2, 2, byrow=TRUE), npNumbers=npNumbers, alternative=alternative,
                                        np.interval=np.interval, beta=beta, method=method, to.plot=FALSE, ref.pvalue=ref.pvalue, delta=delta, reject.alpha=alpha))
        }
      } else { moreExtreme[i] <- FALSE }
    }
    power <- mean(moreExtreme)
  }
  
  #Convert data to power.htest structure
  methodDescribed <- methodText(method, np.interval)
  
  return(structure(list("n1, n2" = c(n1, n2), "p1, p2" = c(p1, p2), alpha = alpha, 
                        power = power, alternative = alternative, delta = delta, method=methodDescribed), class = "power.htest"))
}
