binomialCode <-
function(data, alternative, npNumbers, beta, interval, method, cond.row, to.plot, ref.pvalue){
  
  #If conditioning on row, then transpose 2x2 table
  if(cond.row){data <- t(data)}
  
  Ns <- .colSums(data, 2, 2)
  N <- sum(Ns)
  
  #Specify nuisance parameter range
  if(interval){
    tempInt <- binom.CI(sum(data[1,]), N, conf.level=1-beta)
    int <- seq(max(c(0.00001, tempInt[1])),min(c(0.99999, tempInt[2])), length=npNumbers)
  } else {int <- seq(0.00001,.99999,length=npNumbers); beta <- 0}
  
  #Calculate the test statistics:
  TX <- switch(method,
               "z-pooled" = zpooled_TX(Ns),
               "z-unpooled" = zunpooled_TX(Ns),
               "boschloo" = boschloo_TX(Ns, alternative),
               "santner and snell" = santner_TX(Ns),
               "csm approximate" = csmApprox_TX(data, alternative, npNumbers, int, beta))

  Tbls <- switch(method,
                 "csm" = csm_Tbls(data, alternative, npNumbers, int, beta),
                 "csm modified" = csmMod_Tbls(data, alternative, npNumbers, int, beta))
  
  #Observed Statistic:
  if(method %in% c("z-pooled", "z-unpooled", "boschloo", "santner and snell", "csm approximate")){
    TX[is.na(TX[,3]),3] <- 0
    TX[,3] <- signif(TX[,3], 12) #Remove rounding errors
    TXO <- TX[TX[,1]==data[1,1] & TX[,2]==data[1,2], 3]
  } else {TXO <- NA}

  #Find tables that have a test statistic as or more extreme than the observed statistic:
  if(method %in% c("z-pooled", "z-unpooled", "santner and snell")){
    if(alternative=="greater"){Tbls <- TX[TX[,3] >= TXO,,drop=FALSE]
    } else if(alternative=="less"){Tbls <- TX[TX[,3] <= TXO,,drop=FALSE]
    } else if(alternative=="two.sided"){Tbls <- TX[abs(TX[,3]) >= abs(TXO),,drop=FALSE]}
  } else if (method %in% c("boschloo", "csm approximate")){Tbls <- TX[TX[,3] <= TXO,,drop=FALSE]}
  
  #Search for the maximum p-value:
  if(any(duplicated(Tbls))){stop("TABLES SHOULD NOT BE DUPLICATED: CHECK CODE")}
  maxP <- maxPvalue(Tbls, Ns=Ns, npNumbers=npNumbers, int=int, beta=beta)
  prob <- maxP$prob
  pvalue <- maxP$pvalue
  np <- maxP$np
  
  #Refine the p-value using the optimise function
  if(ref.pvalue){
    refPvalue <- rep(0,length(np))
    refNp <- rep(0,length(np))
    for(i in 1:length(np)){
      ref <- optimise(f=function(p){sum(dbinom(Tbls[,1],Ns[1],p)*dbinom(Tbls[,2],Ns[2],p))},
                    interval=c(max(0.00001,np[i]-1/npNumbers),min(0.99999,np[i]+1/npNumbers)), maximum=TRUE)
      refPvalue[i] <- ref$objective+beta
      refNp[i] <- ref$maximum
    }
    np <- refNp[refPvalue == max(refPvalue)]
    pvalue <- max(refPvalue[refPvalue == max(refPvalue)])
  }
  
  #Note: if beta is large, can have a p-value greater than 1.  Cap at 1
  pvalue[pvalue > 1] <- 1
  
  #Plot p-value vs np
  if(to.plot){
    plot(int, prob+beta, xlim=c(floor(min(int)*10)/10, ceiling(max(int)*10)/10),
         ylim=c(0,max(pvalue)), xlab="np", ylab="P-value", main="P-value as a function of the nuisance parameter")
    points(np, rep(pvalue, length(np)), col="red", pch=21, bg="red")
  }
  
  return(list(method=method, p.value=max(pvalue), test.statistic=TXO, np=np, np.range=c(min(int),max(int))))
}
