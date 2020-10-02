confInt <-
function(conf.level, ESTIMATE, data, alternative, npNumbers, np.interval, beta, method, tsmethod, cond.row, ref.pvalue, delta) {

  alpha <- 1 - conf.level
  
  # Create function that takes the difference between the alpha and p-value and determine when it crosses 0 (can cross multiple times)
  rootFunc <- function(x) {
    vapply(x, function(x){alpha - binomialCode(data, alternative=alternative, npNumbers=npNumbers, np.interval=np.interval, beta=beta,
                                               method=method, tsmethod=tsmethod, cond.row=cond.row, to.plot=FALSE, ref.pvalue=ref.pvalue, delta=x, reject.alpha=NULL)$p.value}, numeric(1))
  }
  
  # uniroot.all function will capture the LB and UB most of the time, but need to check another root isn't missing.  Use while loop until no update
  # pvalueRoots <- uniroot.all(rootFunc, lower=-0.9999, upper=0.9999, n=100)
  
  if (alternative == "greater" || ESTIMATE[[1]] >= 0.9999) { UB <- 1
  } else {
    UB <- 1; prevUB <- 999; lowerVal <- max(c(ESTIMATE, -0.9999))
    while (UB != prevUB) {
      prevUB <- UB  # update prevUB before updating UB
      UBroots <- uniroot.all(rootFunc, lower=lowerVal, upper=0.9999, n=100)
      if (length(UBroots) != 0) { UB <- max(UBroots) }
      lowerVal <- UB + 0.001  #After finding largest root, add a small value and see if there's any other larger root that was missed
    }
  }
  
  if (alternative == "less" || ESTIMATE[[1]] <= -0.9999) { LB <- -1
  } else {
    LB <- -1; prevLB <- 999; lowerVal <- min(c(ESTIMATE, 0.9999))
    while (LB != prevLB) {
      prevLB <- LB  # update prevLB before updating LB
      LBroots <- uniroot.all(rootFunc, lower=-0.9999, upper=lowerVal, n=100)
      if (length(LBroots) != 0) { LB <- min(LBroots) }
      lowerVal <- LB - 0.001  #After finding smallest root, subtract a small value and see if there's any other smaller root that was missed
    }
  }
  
  return(c(LB, UB))
}
