pairedCode <-
function(data, alternative, npNumbers, np.interval, beta, method, tsmethod, to.plot, ref.pvalue, delta, reject.alpha, useStoredCSM) {

  N <- sum(data)
  
  # If alternative is "greater" or two.sided test has p21 - p12 > delta, then swap groups.  This is a major
  # simplification, since code can now only consider cases where the alternative is less than or two.sided with less than always being more
  # extreme.  Will change the test statistic at the end
  swapFlg <- ((alternative == "greater") || (alternative=="two.sided" && (data[1,2]/N - data[2,1]/N) > delta))
  if (swapFlg) {
    dataTemp <- data
    data[1, 2] <- dataTemp[2, 1]
    data[2, 1] <- dataTemp[1, 2]
    delta <- -delta
    if (alternative == "greater") { alternative <- "less" }
  }

  # If two.sided test with central method, then simply calculate one-sided test and double the p-value
  if (alternative == "two.sided" && tsmethod == "central") {
    alternative <- "less"
    doublePvalue <- TRUE
  } else { doublePvalue <- FALSE }
  
  #Specify nuisance parameter range
  if (np.interval && beta != 0) {
    tempInt <- discordant.CI(data[1,2]+data[2,1], N, conf.level=1-beta)
    if (delta == 0) { int <- seq(max(c(0.00001, tempInt[1])), min(c(0.49999, tempInt[2])), length=npNumbers)
    } else if (delta > 0) {
      # The LB may be impossible if LB + LB + delta > 1.  In that case, change int to be lowest possible value
      if (((1-delta)/2 - 0.00001) < tempInt[1]) { int <- (1-delta)/2 - 0.00001
      } else {
        int <- seq(max(c(0.00001, tempInt[1])), min(c((1-delta)/2 - 0.00001, tempInt[2])), length=npNumbers)
      }
    } else if (delta < 0) {
      # The UB may be impossible if UB + UB + delta < -1.  In that case, change int to be highest possible value
      if ((-delta + 0.00001) > tempInt[2]) { int <- -delta + 0.00001
      } else {
        int <- seq(max(c(-delta + 0.00001, tempInt[1])), min(c((1 - delta)/2 - 0.00001, tempInt[2])), length=npNumbers)
      }
    }
  } else {
    if (delta == 0) { int <- seq(0.00001,.49999,length=npNumbers)
    } else if (delta > 0) { int <- seq(0.00001, (1-delta)/2 - 0.00001, length=npNumbers)
    } else if (delta < 0) { int <- seq(-delta + 0.00001, (1 - delta)/2 - 0.00001, length=npNumbers)}
    beta <- 0
  }
  
  #Find tables that have a test statistic as or more extreme than the observed statistic:
  if (method == "csm") {
    
    # Use stored ordering matrix if available
    if (N <= 200 && delta == 0 && useStoredCSM) {

      if (!requireNamespace("ExactData", quietly = TRUE)) {
        stop(paste("ExactData R package must be installed when useStoredCSM=TRUE. To install ExactData R package, run:",
                   "`install.packages('ExactData', repos='https://pcalhoun1.github.io/drat/', type='source')`"))
      }
      
      if (alternative == "two.sided") { orderMat <- ExactData::orderCSMPairedMatTwoSided[[N]]
      } else { orderMat <- ExactData::orderCSMPairedMatOneSided[[N]] }
      
      orderVec <- as.vector(orderMat)[!is.na(orderMat)]
      TX <- matrix(cbind(which(!is.na(orderMat), arr.ind = TRUE) - 1, orderVec), ncol=3, dimnames = NULL)
      TX[ , 3] <- signif(TX[ , 3], 12)  #Remove rounding errors
      
      TXO <- TX[TX[,1]==data[1,2] & TX[,2]==data[2,1], 3]
      
      rejectFlg <- (TX[,3] <= TXO)
      
      moreExtremeMat <- matrix(0, N+1, N+1, dimnames=list(0:N, 0:N))
      moreExtremeMat <- lowerMatVal(moreExtremeMat, NA)
      moreExtremeMat[!is.na(moreExtremeMat)] <- rejectFlg
      
      # Set TXO to be NA (not ordering step number)
      findMoreExtreme <- list(TXO=NA, moreExtremeMat=moreExtremeMat)
    } else {
      findMoreExtreme <- moreExtremePairedCSM(data = data, N = N, alternative = alternative,
                                              int = int, doublePvalue = doublePvalue, delta = delta, reject.alpha = reject.alpha)
    }
    # if data is empty or findMoreExtreme is just FALSE, then return findMoreExtreme
    if (is.null(data) || isFALSE(findMoreExtreme)) { return(findMoreExtreme) }
  } else {
    findMoreExtreme <- moreExtremePaired(method=method, data=data, N=N, alternative=alternative, int=int, delta=delta)
  }
  
  #Search for the maximum p-value:
  maxP <- maxPvaluePaired(findMoreExtreme$moreExtremeMat, N, int, beta, delta, doublePvalue)
  
  prob <- maxP$prob
  pvalue <- maxP$pvalue
  np <- maxP$np
  if (!is.null(reject.alpha) && pvalue > reject.alpha) { return(FALSE) }
  
  #Refine the p-value using the optimise function
  if (ref.pvalue && length(int) > 1) {
    Tbls <- which(findMoreExtreme$moreExtremeMat==1, arr.ind = TRUE) - 1
    refPvalue <- rep(0,length(np))
    refNp <- rep(0,length(np))
    for (i in 1:length(np)) {
      ref <- suppressWarnings(
        optimise(f=function(p){sum(trinom(Tbls[ , 1], Tbls[ , 2], N, p12=p, p21=p, delta=delta))},
                 interval=c(max(int[1],np[i]-1/npNumbers),min(int[npNumbers],np[i]+1/npNumbers)), maximum=TRUE))
      refPvalue[i] <- ref$objective
      refNp[i] <- ref$maximum
    }
    refPvalue <- min(c(1, signif((1+doublePvalue)*refPvalue + beta, 12))) #Remove rounding errors and cap at 1
    
    if (!all(is.na(refPvalue)) && max(refPvalue, na.rm=TRUE) > pvalue) {
      np <- refNp[refPvalue == max(refPvalue)]
      pvalue <- max(refPvalue)
    }
    
  }
  
  if (!is.null(reject.alpha)) { return(pvalue <= reject.alpha) }
  
  #Plot p-value vs np
  if (to.plot) {
    plot(int, prob, xlim=c(floor(min(int)*10)/10, ceiling(max(int)*10)/10),
         ylim=c(0,max(pvalue)), xlab="np", ylab="P-value", main="P-value as a function of the nuisance parameter")
    points(np, rep(pvalue, length(np)), col="red", pch=21, bg="red")
  }
  
  TXO <- ifelse(swapFlg && method %in% c("uam", "uamcc"), -findMoreExtreme$TXO, findMoreExtreme$TXO)
  return(list(method=method, p.value=pvalue, test.statistic=TXO, np=np, np.range=c(min(int),max(int))))
}
