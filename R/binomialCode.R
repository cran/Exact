binomialCode <-
function(data, alternative, npNumbers, np.interval, beta, method, tsmethod, to.plot, ref.pvalue, delta, reject.alpha, useStoredCSM){
  
  Ns <- .colSums(data, 2, 2)
  N <- sum(Ns)
  
  if (any(Ns <= 0)) { stop("Can't have a sample size of 0 for one of the groups") }
  
  # If alternative is "greater" or two.sided test has p1 - p2 > delta, then swap groups.  This is a major
  # simplification, since code can now only consider cases where the alternative is less than or two.sided with less than always being more
  # extreme.  Will change the test statistic at the end
  swapFlg <- ((alternative == "greater") || (alternative=="two.sided" && (data[1,1]/Ns[1] - data[1,2]/Ns[2]) > delta))
  if (swapFlg) {
    data <- data[ , 2:1]
    Ns <- Ns[2:1]
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
    tempInt <- binom.CI(sum(data[1,]), N, conf.level=1-beta)
    if (delta == 0) { int <- seq(max(c(0.00001, tempInt[1])), min(c(0.99999, tempInt[2])), length=npNumbers)
    } else if (delta > 0) {
      int <- seq(max(c(0.00001, tempInt[1])), min(c(1 - delta - 0.00001, tempInt[2])), length=npNumbers)
    } else if (delta < 0) {
      int <- seq(max(c(abs(delta) + 0.00001, tempInt[1])), min(c(0.99999, tempInt[2])), length=npNumbers)}
  } else {
    if (delta == 0) { int <- seq(0.00001,.99999,length=npNumbers)
    } else if (delta > 0) { int <- seq(0.00001, 1 - delta - 0.00001, length=npNumbers)
    } else if (delta < 0) { int <- seq(abs(delta) + 0.00001, .99999, length=npNumbers)}
    beta <- 0
  }
  
  #Find tables that have a test statistic as or more extreme than the observed statistic:
  if (method == "csm") {
    
    # Use stored ordering matrix if available
    if (delta == 0 && useStoredCSM &&
         ( (max(Ns) <= 100) || ((Ns[1]+1)*(Ns[2]+1) <= 15000 && Ns[1] <= 2*Ns[2] && Ns[2] <= 2*Ns[1]) ) ) {
      
      if (!requireNamespace("ExactData", quietly = TRUE)) {
        stop(paste("ExactData R package must be installed when useStoredCSM=TRUE. To install ExactData R package, run:",
                   "`install.packages('ExactData', repos='https://pcalhoun1.github.io/drat/', type='source')`"))
      }
      
      if (alternative == "two.sided") {
        if (Ns[1] < Ns[2]) {  # Must convert known matrix of (n2,n1) to (n1,n2)
          orderMat <- t(ExactData::orderCSMMatTwoSided[[paste0("(",Ns[2],",",Ns[1],")")]])
          orderMat <- orderMat[(Ns[1]+1):1, (Ns[2]+1):1]
          rownames(orderMat) <- 0:Ns[1]
          colnames(orderMat) <- 0:Ns[2]
        } else { orderMat <- ExactData::orderCSMMatTwoSided[[paste0("(",Ns[1],",",Ns[2],")")]] }
        
      } else {
        if (Ns[1] < Ns[2]) {  # Must convert known matrix of (n2,n1) to (n1,n2)
          orderMat <- t(ExactData::orderCSMMatOneSided[[paste0("(",Ns[2],",",Ns[1],")")]])
          orderMat <- orderMat[(Ns[1]+1):1, (Ns[2]+1):1]
          rownames(orderMat) <- 0:Ns[1]
          colnames(orderMat) <- 0:Ns[2]
        } else { orderMat <- ExactData::orderCSMMatOneSided[[paste0("(",Ns[1],",",Ns[2],")")]] }
      }
      
      TX <- matrix(cbind(as.matrix(expand.grid(0:Ns[1], 0:Ns[2])), as.vector(orderMat)), ncol=3, dimnames = NULL)
      TX <- TX[order(TX[,1], TX[,2]), ]
      TX[ , 3] <- signif(TX[ , 3], 12)  #Remove rounding errors
      
      TXO <- TX[TX[,1]==data[1,1] & TX[,2]==data[1,2], 3]
      
      rejectFlg <- (TX[,3] <= TXO)
      moreExtremeMat <- matrix(rejectFlg, Ns[1]+1, Ns[2]+1, byrow=TRUE, dimnames=list(0:Ns[1], 0:Ns[2]))*1
      
      # Set TXO to be NA (not ordering step number)
      findMoreExtreme <- list(TXO=NA, moreExtremeMat=moreExtremeMat)
      
    } else {
      findMoreExtreme <- moreExtremeCSM(data = data, Ns = Ns, alternative = alternative,
                                        int = int, doublePvalue = doublePvalue, delta = delta, reject.alpha = reject.alpha)
    }
    # if data is empty or findMoreExtreme is just FALSE, then return findMoreExtreme
    if (is.null(data) || isFALSE(findMoreExtreme)) { return(findMoreExtreme) }
  } else {
    findMoreExtreme <- moreExtreme(method=method, data=data, Ns=Ns, alternative=alternative, int=int, delta=delta)
  }
  
  #Search for the maximum p-value:
  maxP <- maxPvalue(findMoreExtreme$moreExtremeMat, Ns, int, beta, delta, doublePvalue)
  
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
  
  TXO <- ifelse(swapFlg && method %in% c("z-pooled", "z-unpooled", "santner and snell"), -findMoreExtreme$TXO, findMoreExtreme$TXO)
  return(list(method=method, p.value=pvalue, test.statistic=TXO, np=np, np.range=c(min(int),max(int))))
}
