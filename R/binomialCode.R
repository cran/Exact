binomialCode <-
function(data, alternative, npNumbers, np.interval, beta, method, cond.row, to.plot, ref.pvalue, delta, reject.alpha){
  
  #If conditioning on row, then transpose 2x2 table
  if (cond.row) { data <- t(data) }
  
  Ns <- .colSums(data, 2, 2)
  N <- sum(Ns)
  
  if (any(Ns <= 0)) { stop("Can't have a sample size of 0 for one of the groups") }
  
  # If alternative is "greater" or two.sided test has first group with larger proportion, then swap groups.  This is a major
  # simplification, since code can now only consider cases where the alternative is less than or two.sided with less than always being more
  # extreme.  We do have to change the test statistic though at the end
  swapFlg <- ((alternative == "greater") || (alternative=="two.sided" && (data[1,1]/Ns[1] - data[1,2]/Ns[2]) > delta))
  if (swapFlg) {
    data <- data[ , 2:1]
    Ns <- Ns[2:1]
    delta <- -delta
    if (alternative == "greater") { alternative <- "less" }
  }
  
  #Specify nuisance parameter range
  if (np.interval) {
    tempInt <- binom.CI(sum(data[1,]), N, conf.level=1-beta)
    int <- seq(max(c(0.00001, tempInt[1])),min(c(0.99999, tempInt[2])), length=npNumbers)
  } else {
    #Note: can't use non-zero delta with np.interval approach
    if (delta == 0) { int <- seq(0.00001,.99999,length=npNumbers)
    } else if (delta > 0) { int <- seq(0.00001, 1 - delta - 0.00001, length=npNumbers)
    } else if (delta < 0) { int <- seq(abs(delta) + 0.00001, .99999, length=npNumbers)}
    beta <- 0
  }
  
  #Find tables that have a test statistic as or more extreme than the observed statistic:
  if (method == "csm") {
    findMoreExtreme <- moreExtremeCSM(data=data, Ns=Ns, alternative=alternative, int=int, delta=delta, reject.alpha=reject.alpha)
    # For non-null reject.alpha, then findMoreExtreme can return FALSE.  
    if (is.logical(findMoreExtreme)) { return(FALSE) }
  } else {
    findMoreExtreme <- moreExtreme(method=method, data=data, Ns=Ns, alternative=alternative, int=int, delta=delta)
  }
  
  #Search for the maximum p-value:
  #lookupArray <- dbinomCalc(Ns, int, delta)
  #Tbls <- which(findMoreExtreme$moreExtremeMat==1, arr.ind = TRUE) - 1
  #maxPvalueLookup(Tbls, int=int, lookupArray)
  maxP <- maxPvalue(findMoreExtreme$moreExtremeMat, Ns, int, beta, delta)
  
  prob <- maxP$prob
  pvalue <- maxP$pvalue
  np <- maxP$np
  if (!is.null(reject.alpha) && pvalue > reject.alpha) { return(FALSE) }
  
  #Refine the p-value using the optimise function
  if (ref.pvalue) {
    refPvalue <- rep(0,length(np))
    refNp <- rep(0,length(np))
    for (i in 1:length(np)) {
      ref <- suppressWarnings(
        optimise(f=function(p){sum(matrix(dbinom(0:Ns[1], Ns[1], p+delta),ncol=1)*(findMoreExtreme$moreExtremeMat %*% dbinom(0:Ns[2], Ns[2], p)))},
                 interval=c(max(int[1],np[i]-1/npNumbers),min(int[npNumbers],np[i]+1/npNumbers)), maximum=TRUE))
      refPvalue[i] <- ref$objective+beta
      refNp[i] <- ref$maximum
    }
    refPvalue <- signif(refPvalue, 12) #Remove rounding errors
    if (!all(is.na(refPvalue)) && max(refPvalue, na.rm=TRUE) > pvalue) {
      np <- refNp[refPvalue == max(refPvalue)]
      pvalue <- max(refPvalue)
    }
  }
  
  #Note: if beta is large, can have a p-value greater than 1.  Cap at 1
  pvalue[pvalue > 1] <- 1
  if (!is.null(reject.alpha)) {return(pvalue <= reject.alpha)}
  
  #Plot p-value vs np
  if (to.plot) {
    plot(int, prob+beta, xlim=c(floor(min(int)*10)/10, ceiling(max(int)*10)/10),
         ylim=c(0,max(pvalue)), xlab="np", ylab="P-value", main="P-value as a function of the nuisance parameter")
    points(np, rep(pvalue, length(np)), col="red", pch=21, bg="red")
  }
  
  TXO <- ifelse(swapFlg && method %in% c("z-pooled", "z-unpooled", "santner and snell"), -findMoreExtreme$TXO, findMoreExtreme$TXO)
  return(list(method=method, p.value=pvalue, test.statistic=TXO, np=np, np.range=c(min(int),max(int))))
}
