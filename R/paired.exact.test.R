paired.exact.test <-
function(data, alternative=c("two.sided", "less", "greater"), npNumbers=100, np.interval=FALSE, beta=0.001,
                       method=c("uam", "ucm", "uamcc", "csm"),
                       tsmethod=c("square", "central"),
                       conf.int=FALSE, conf.level=0.95,
                       to.plot=TRUE, ref.pvalue=TRUE, delta=0, reject.alpha=NULL, useStoredCSM=TRUE){
  
  alternative <- match.arg(tolower(alternative), c("two.sided", "less", "greater"))
  method <- convertMethod(method)
  method <- match.arg(method, c("uam", "ucm", "uamcc", "csm")) #paired.exact.test doesn't include all methods
  tsmethod <- match.arg(tolower(tsmethod), c("square", "central"))
  
  #Perform several checks
  checkPairedParam(data=data, alternative=alternative, npNumbers=npNumbers, np.interval=np.interval, beta=beta,
                   method=method, tsmethod=tsmethod, conf.int=conf.int, conf.level=conf.level,
                   to.plot=to.plot, ref.pvalue=ref.pvalue, delta=delta, reject.alpha=reject.alpha, useStoredCSM=useStoredCSM)
  
  results <- pairedCode(data, alternative=alternative, npNumbers=npNumbers, np.interval=np.interval, beta=beta,
                        method=method, tsmethod=tsmethod, to.plot=to.plot, ref.pvalue=ref.pvalue, delta=delta, reject.alpha=reject.alpha,
                        useStoredCSM=useStoredCSM)
  
  if (!is.null(reject.alpha)) { return(results) }
    
  #Convert data to htest structure
  N <- sum(data)
  ESTIMATE <- data[1,2]/N - data[2,1]/N
  data.name <- paste(rowSums(data)[1], 'out of', N, 'vs.', colSums(data)[1], 'out of', N)
  
  null.value <- delta
  names(N) <- "total sample size"
  names(ESTIMATE) <- names(null.value) <- "difference in proportion"
  names(results$test.statistic) <- "test statistic"
    
  np <- list(results$np)
  names(np) <- "nuisance parameter"
  np.range <- list(results$np.range)
  names(np.range) <- "nuisance parameter range"
  methodDescribed <- pairedMethodText(method, np.interval)
    
  # Calculate confidence interval
  if (conf.int) {
    CINT <- confIntPaired(conf.level, ESTIMATE, data, alternative, npNumbers, np.interval, beta, method, tsmethod, ref.pvalue, delta, useStoredCSM)
    attr(CINT, "conf.level") <- conf.level
  }
    
  return(structure(list(statistic = results$test.statistic, parameter = N, p.value = results$p.value,
                        conf.int = if (conf.int) CINT, 
                        estimate = ESTIMATE, null.value = null.value, alternative = alternative, 
                        np=np, np.range=np.range, method = methodDescribed, tsmethod = tsmethod,
                        data.name = data.name), class = "htest"))
  
}
