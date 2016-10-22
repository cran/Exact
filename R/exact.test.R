exact.test <-
function(data, alternative=c("two.sided", "less", "greater"), npNumbers=100, beta=0.001, interval=FALSE,
                       method=c("z-pooled", "z-unpooled", "boschloo", "santner and snell", "csm", "csm approximate", "csm modified"),
                       model=c("Binomial","Multinomial"), cond.row=TRUE, to.plot=TRUE, ref.pvalue=TRUE){

  #Perform several checks
  stopifnot(is.logical(interval) && is.logical(cond.row) && is.logical(to.plot) && is.logical(ref.pvalue))
  if(interval && (beta < 0 || beta > 1)){stop("Beta must be between 0 and 1")}
  if(npNumbers < 1){stop("Total number of nuisance parameters considered must be at least 1")}

  alternative <- match.arg(tolower(alternative), c("two.sided", "less", "greater"))
  
  #Sometimes Z-pooled is called score and Z-unpooled is called wald statistic
  if(length(method)==1 && tolower(method)=="score"){method <- "z-pooled"}
  if(length(method)==1 && tolower(method)=="wald"){method <- "z-unpooled"}
  method <- match.arg(tolower(method), c("z-pooled", "z-unpooled", "boschloo", "santner and snell",
                                         "csm", "csm approximate", "csm modified"))
  model <- match.arg(tolower(model), c("binomial","multinomial"))
  
  if(nrow(data) != 2 || ncol(data) != 2){stop("Input 2x2 table")}
  if(any(data < 0)){stop("Can't have negative entries")}
  if(any(.colSums(data, 2, 2)==0) || any(.rowSums(data, 2, 2)==0)){stop("Can't have all 0's for row or column")}
  if(!all(as.integer(data)==data)){stop("Data must only contain integers")}
  
  if(model=="multinomial" && method %in% c("csm","csm approximate","csm modified")){stop("Code currently cannot implement CSM tests for multinomial models")}
  
  if(method %in% c("csm","csm modified")){
    if(!cond.row){
      if(data[1,1]/sum(data[,1]) == data[1,2]/sum(data[,2])){stop("For CSM, proportions must not be equal. P-value would be high")}
      if(alternative=="less" && (data[1,1]/sum(data[,1]) >= data[1,2]/sum(data[,2]))){stop("For CSM, proportions must match direction of alternative hypothesis. P-value would be high")}
      if(alternative=="greater" && (data[1,1]/sum(data[,1]) <= data[1,2]/sum(data[,2]))){stop("For CSM, proportions must match direction of alternative hypothesis. P-value would be high")}
    } else {
      if(data[1,1]/sum(data[1,]) == data[2,1]/sum(data[2,])){stop("For CSM, proportions must not be equal. P-value would be high anyways")}
      if(alternative=="less" && (data[1,1]/sum(data[1,]) >= data[2,1]/sum(data[2,]))){stop("For CSM, proportions must match direction of alternative hypothesis. P-value would be high")}
      if(alternative=="greater" && (data[1,1]/sum(data[1,]) <= data[2,1]/sum(data[2,]))){stop("For CSM, proportions must match direction of alternative hypothesis. P-value would be high")}
    }
  }

  if(model=="binomial"){
    results <- binomialCode(data, alternative=alternative, interval=interval, beta=beta, npNumbers=npNumbers,
                            method=method, cond.row=cond.row, to.plot=to.plot, ref.pvalue=ref.pvalue)
    
    #Convert data to htest structure
    if(cond.row){
      Ns <- .rowSums(data, 2, 2)
      ESTIMATE <- data[1,1]/Ns[1]-data[2,1]/Ns[2]
      data.name <- paste(data[1,1], 'out of', Ns[1], 'vs.', data[2,1], 'out of',Ns[2])
    } else {
      Ns <- .colSums(data, 2, 2)
      ESTIMATE <- data[1,1]/Ns[1]-data[1,2]/Ns[2]
      data.name <- paste(data[1,1], 'out of', Ns[1], 'vs.', data[1,2], 'out of',Ns[2])
    }
    null.value <- 0
    names(Ns) <- c("first sample size", "second sample size")
    names(ESTIMATE) <- names(null.value) <- "difference in proportion"
    names(results$test.statistic) <- "test statistic"
    
    np <- results$np
    names(np) <- "nuisance parameter"
    np.range <- list(results$np.range)
    names(np.range) <- "nuisance parameter range"
    if(interval){methodUsed <- paste(method, "with interval")
    } else {methodUsed <- method}

    return(structure(list(statistic = results$test.statistic, parameter = Ns, p.value = results$p.value,
                   estimate = ESTIMATE, null.value = null.value, alternative = alternative, 
                   np=np, np.range=np.range, model = model, method = methodUsed, 
                   data.name = data.name), class = "htest"))
  }
  
  if(model=="multinomial"){
    results <- multinomialCode(data, alternative=alternative, interval=interval, beta=beta, npNumbers=npNumbers, method=method)
    
    #Convert data to htest structure
    N <- sum(data)
    ESTIMATE <- (data[1,1]*data[2,2]-data[1,2]*data[2,1])/N^2
    data.name <- paste0("x11=",data[1,1],", x12=",data[1,2],", x21=",data[2,1],", x22=",data[2,2])
    null.value <- 0
    names(N) <- "total sample size"
    names(ESTIMATE) <- names(null.value) <- "difference in product of proportion"
    names(results$test.statistic) <- "test statistic"
    
    np <- c(results$np1, results$np2)
    names(np) <- c("first nuisance parameter", "second nuisance parameter")
    np.range <- list(results$np1.range, results$np2.range)
    names(np.range) <- c("first nuisance parameter range", "second nuisance parameter range")
    if(interval){methodUsed <- paste(method, "with interval")
    } else {methodUsed <- method}
    
    return(structure(list(statistic = results$test.statistic, parameter = N, p.value = results$p.value,
                          estimate = ESTIMATE, null.value = null.value, alternative = alternative, 
                          np=np, np.range=np.range, model = model, method = methodUsed, 
                          data.name = data.name), class = "htest"))
  }
}
