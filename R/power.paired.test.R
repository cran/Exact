power.paired.test <-
function(p12, p21, N, alternative=c("two.sided", "less", "greater"), alpha=0.05,
                              npNumbers=100, np.interval=FALSE, beta=0.001,
                              method=c("uam", "ucm", "uamcc", "csm", "cm", "am", "amcc"),
                              tsmethod=c("square", "central"), simulation=FALSE, nsim = 100, delta=0, convexity=TRUE, useStoredCSM=TRUE){
  
  alternative <- match.arg(tolower(alternative), c("two.sided", "less", "greater"))
  method <- convertMethod(method)
  tsmethod <- match.arg(tolower(tsmethod), c("square", "central"))

  #Perform several checks
  checkPairedParam(p12=p12, p21=p21, N=N, alternative=alternative, alpha=alpha, npNumbers=npNumbers, np.interval=np.interval, beta=beta,
                   method=method, tsmethod=tsmethod, simulation=simulation, nsim = nsim, delta=delta, convexity=convexity, useStoredCSM=useStoredCSM)
  
  # Consider all tables:
  if (!simulation) {
    rejectRegion <- paired.reject.region(N = N, alternative=alternative, alpha=alpha,
                                         npNumbers=npNumbers, np.interval=np.interval, beta=beta,
                                         method=method, tsmethod=tsmethod, delta=delta, convexity=convexity,
                                         useStoredCSM=useStoredCSM)
    prob <- rejectRegion #temporarily set it to be the same matrix style as reject region
    index <- 1
    for (x12 in 0:N) {
      # Note: don't use delta for probabilities
      prob[index, 1:(N+1-x12)] <- trinom(x12, 0:(N-x12), N, p12, p21, delta=0)
      index = index + 1
    }
    power <- sum(prob[!is.na(rejectRegion) & as.logical(rejectRegion)])
    
  } else if (simulation) {
    #Randomly generate a table based on known proportions
    randDiscord <- rbinom(nsim, size = N, prob = p12 + p21)
    randX12 <- rbinom(nsim, size = randDiscord, prob = p12/(p12 + p21))
    randTables <- cbind(N-randDiscord, randX12, randDiscord-randX12, 0, deparse.level = 0)
    moreExtremeSim <- rep(NA, nsim)
    for (i in 1:nsim) {
      if (delta != 0 || (alternative=="greater" && (randTables[i,2]/N - randTables[i,3]/N) > 0) ||
          (alternative=="less" && (randTables[i,2]/N - randTables[i,3]/N) < 0) || 
          (alternative=="two.sided" && randTables[i,2] != randTables[i,3])) {
        
        dat <- matrix(randTables[i,], ncol=2, nrow=2, byrow=TRUE)
        if (method %in% c("am", "amcc")) {
          # Note: mcnemar.test() only works for two.sided test
          TXO <- switch(method,
                        "am" = mcnemar_TX(dat, N, delta=delta, CC=FALSE),
                        "amcc" = mcnemar_TX(dat, N, delta=delta, CC=TRUE))[3]
          moreExtremeSim[i] <- (ifelse(alternative=="two.sided",2,1)*pnorm(TXO) <= alpha)
          #mcnemar.test(dat, correct=FALSE)$p.value
        } else if (method == "cm") {
          moreExtremeSim[i] <- (mcnemar.2x2(dat, N, alternative=alternative, pval1ties = FALSE)[,3] <= alpha)
        } else {
          moreExtremeSim[i] <- (pairedCode(dat, npNumbers=npNumbers, alternative=alternative,
                                        np.interval=np.interval, beta=beta, method=method, tsmethod=tsmethod, to.plot=FALSE,
                                        ref.pvalue=FALSE, delta=delta, reject.alpha=alpha, useStoredCSM=useStoredCSM))
        }
      } else { moreExtremeSim[i] <- FALSE }
    }
    power <- mean(moreExtremeSim)
  }
  
  #Convert data to power.htest structure
  methodDescribed <- pairedMethodText(method, np.interval)
  
  return(structure(list(N = N, "p12, p21" = c(p12, p21), alpha = alpha, 
                        power = power, alternative = alternative, delta = delta, method=methodDescribed), class = "power.htest"))
}
