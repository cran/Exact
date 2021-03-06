exact.test <-
function(data, alternative=c("two.sided", "less", "greater"), npNumbers=100, np.interval=FALSE, beta=0.001,
                       method=c("z-pooled", "z-unpooled", "boschloo", "santner and snell", "csm", "csm approximate"),
                       model=c("Binomial", "Multinomial"), tsmethod=c("square", "central"),
                       conf.int=FALSE, conf.level=0.95,
                       cond.row=TRUE, to.plot=TRUE, ref.pvalue=TRUE, delta=0, reject.alpha=NULL){
  
  #Perform several checks
  stopifnot(is.logical(np.interval) && is.logical(cond.row) && is.logical(to.plot) && is.logical(ref.pvalue) && is.logical(conf.int))
  if (np.interval && (beta < 0 || beta > 1)) { stop("Beta must be between 0 and 1") }
  if (conf.int && (conf.level <= 0 || conf.level >= 1)) { stop("conf.level must be between 0 and 1") }
  if (npNumbers < 1) { stop("Total number of nuisance parameters considered must be at least 1") }
  if (!is.null(reject.alpha) && (reject.alpha < 0 || reject.alpha > 1)) { stop("reject.alpha must be between 0 and 1") }
  
  alternative <- match.arg(tolower(alternative), c("two.sided", "less", "greater"))
  
  # The Z-pooled statistic that calculates the variance using MLE, which is the pooled variance if delta=0.
  # The Z-pooled statistic is also (perhaps better) known as the Score statistic
  # The classic z-pooled statistic is not performed as the performance is inferior when delta != 0
  if(length(method)==1 && tolower(method)=="score"){method <- "z-pooled"}
  #if(length(method)==1 && tolower(method)=="wald"){method <- "z-unpooled"}
  method <- match.arg(tolower(method), c("z-pooled", "z-unpooled", "boschloo", "santner and snell", "csm", "csm approximate"))
  model <- match.arg(tolower(model), c("binomial", "multinomial"))
  tsmethod <- match.arg(tolower(tsmethod), c("square", "central"))
  
  if (nrow(data) != 2 || ncol(data) != 2) { stop("Input 2x2 table") }
  if (any(data < 0)) { stop("Can't have negative entries") }
  if (!all(as.integer(data) == data)) { stop("Data must only contain integers") }
  
  if (delta != 0 && model == "multinomial") {
    stop("Nonzero delta only implemented for binomial models")
  }

  # Print special case that cannot be performed
  if (alternative=="two.sided" && tsmethod == "square" && method %in% c("boschloo", "csm approximate") && (delta != 0 || conf.int)) {
    stop(paste0(method, " with two-sided CIs and/or nonzero delta cannot be implemented when tsmethod='square'.  Suggest using tsmethod='central'"))
  }
  
  if (method %in% c("csm", "csm approximate")) {
    if (model=="multinomial") { stop("Code currently cannot implement CSM tests for multinomial models") }
    if (np.interval) {
      warning("Interval of nuisance parameter cannot be used with CSM; np.interval changed to FALSE")
      np.interval <- FALSE
    }
  }
  
  if (model=="binomial") {
    
    results <- binomialCode(data, alternative=alternative, np.interval=np.interval, beta=beta, npNumbers=npNumbers,
                            method=method, tsmethod=tsmethod, cond.row=cond.row, to.plot=to.plot, ref.pvalue=ref.pvalue, delta=delta, reject.alpha=reject.alpha)
    
    if (!is.null(reject.alpha)) { return(results) }
    
    #Convert data to htest structure
    if (cond.row) {
      Ns <- .rowSums(data, 2, 2)
      ESTIMATE <- data[1,1]/Ns[1]-data[2,1]/Ns[2]
      data.name <- paste(data[1,1], 'out of', Ns[1], 'vs.', data[2,1], 'out of',Ns[2])
    } else {
      Ns <- .colSums(data, 2, 2)
      ESTIMATE <- data[1,1]/Ns[1]-data[1,2]/Ns[2]
      data.name <- paste(data[1,1], 'out of', Ns[1], 'vs.', data[1,2], 'out of',Ns[2])
    }
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
      CINT <- confInt(conf.level, ESTIMATE, data, alternative, npNumbers, np.interval, beta, method, tsmethod, cond.row, ref.pvalue, delta)
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
