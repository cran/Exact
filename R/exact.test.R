exact.test <-
function(data, alternative=c("two.sided", "less", "greater"), npNumbers=100, np.interval=FALSE, beta=0.001,
                       method=c("z-pooled", "z-unpooled", "boschloo", "santner and snell", "csm"),
                       model=c("Binomial", "Multinomial"), tsmethod=c("square", "central"),
                       conf.int=FALSE, conf.level=0.95,
                       cond.row=TRUE, to.plot=TRUE, ref.pvalue=TRUE, delta=0, reject.alpha=NULL, useStoredCSM=TRUE){
  
  alternative <- match.arg(tolower(alternative), c("two.sided", "less", "greater"))
  # The Z-pooled statistic is actually the Score statistic, which are equivalent when delta = 0
  # The classic Z-pooled statistic is not performed as the performance is inferior when delta != 0
  if (length(method)==1 && tolower(method)=="score") { method <- "z-pooled" }
  method <- match.arg(tolower(method), c("z-pooled", "z-unpooled", "boschloo", "santner and snell", "csm"))
  model <- match.arg(tolower(model), c("binomial", "multinomial"))
  tsmethod <- match.arg(tolower(tsmethod), c("square", "central"))

  #Perform several checks
  checkParam(data=data, alternative=alternative, npNumbers=npNumbers, np.interval=np.interval, beta=beta,
             method=method, model=model, tsmethod=tsmethod, conf.int=conf.int, conf.level=conf.level,
             cond.row=cond.row, to.plot=to.plot, ref.pvalue=ref.pvalue, delta=delta, reject.alpha=reject.alpha, useStoredCSM=useStoredCSM)
  
  if (model=="binomial") {

    #If conditioning on row, then transpose 2x2 table
    if (cond.row) { data <- t(data) }
    
    results <- binomialCode(data, alternative=alternative, np.interval=np.interval, beta=beta, npNumbers=npNumbers,
                            method=method, tsmethod=tsmethod, to.plot=to.plot, ref.pvalue=ref.pvalue, delta=delta,
                            reject.alpha=reject.alpha, useStoredCSM=useStoredCSM)
    
    if (!is.null(reject.alpha)) { return(results) }
    
    #Convert data to htest structure
    Ns <- .colSums(data, 2, 2)
    ESTIMATE <- data[1,1]/Ns[1]-data[1,2]/Ns[2]
    data.name <- paste(data[1,1], 'out of', Ns[1], 'vs.', data[1,2], 'out of', Ns[2])
    
    null.value <- delta
    names(Ns) <- c("first sample size", "second sample size")
    names(ESTIMATE) <- names(null.value) <- "difference in proportion"
    names(results$test.statistic) <- "test statistic"
    
    np <- list(results$np)
    names(np) <- "nuisance parameter"
    np.range <- list(results$np.range)
    names(np.range) <- "nuisance parameter range"
    methodDescribed <- methodText(method, np.interval)
    
    # Calculate confidence interval
    if (conf.int) {
      CINT <- confInt(conf.level, ESTIMATE, data, alternative, npNumbers, np.interval, beta, method, tsmethod, ref.pvalue, delta, useStoredCSM)
      attr(CINT, "conf.level") <- conf.level
    }
    
    return(structure(list(statistic = results$test.statistic, parameter = Ns, p.value = results$p.value,
                          conf.int = if (conf.int) CINT, 
                          estimate = ESTIMATE, null.value = null.value, alternative = alternative, 
                          np=np, np.range=np.range, model = model, method = methodDescribed, tsmethod = tsmethod,
                          data.name = data.name), class = "htest"))
  }
  
  if (model=="multinomial") {
    results <- multinomialCode(data, alternative=alternative, np.interval=np.interval, beta=beta, npNumbers=npNumbers, method=method)
    
    #Convert data to htest structure
    N <- sum(data)
    ESTIMATE <- (data[1,1]*data[2,2]-data[1,2]*data[2,1])/N^2
    data.name <- paste0("x11=",data[1,1],", x12=",data[1,2],", x21=",data[2,1],", x22=",data[2,2])
    null.value <- 0
    names(N) <- "total sample size"
    names(ESTIMATE) <- names(null.value) <- "difference in product of proportion"
    names(results$test.statistic) <- "test statistic"
    
    np <- list(results$np1, results$np2)
    names(np) <- c("first nuisance parameter", "second nuisance parameter")
    np.range <- list(results$np1.range, results$np2.range)
    names(np.range) <- c("first nuisance parameter range", "second nuisance parameter range")
    methodDescribed <- methodText(method, np.interval)
    
    return(structure(list(statistic = results$test.statistic, parameter = N, p.value = results$p.value,
                          estimate = ESTIMATE, null.value = null.value, alternative = alternative, 
                          np=np, np.range=np.range, model = model, method = methodDescribed, 
                          data.name = data.name), class = "htest"))
  }
}
