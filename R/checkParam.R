checkParam <-
function(data=NULL, p1=NULL, p2=NULL, n1=NULL, n2=NULL, method="", model="",
                       alternative="", tsmethod="", npNumbers=100, np.interval=FALSE, beta=0.001,
                       conf.int=FALSE, conf.level=0.95,
                       cond.row=TRUE, to.plot=TRUE, ref.pvalue=TRUE, delta=0, reject.alpha=NULL,
                       alpha=0.05, simulation=FALSE, nsim = 100, convexity=TRUE, useStoredCSM=TRUE) {

  ### Logical checks ###
  stopifnot(is.logical(np.interval) && is.logical(cond.row) && is.logical(to.plot) &&
            is.logical(ref.pvalue) && is.logical(conf.int) && is.logical(useStoredCSM) &&
            is.logical(simulation) && is.logical(convexity))
  
  ### Impossible values checks ###
  if (!is.null(data)) {
    if (nrow(data) != 2 || ncol(data) != 2) { stop("Input 2x2 table") }
    if (any(data < 0)) { stop("Data cannot have negative entries") }
    if (!all(as.integer(data) == data)) { stop("Data must only contain integers") }
  }
  if ((!is.null(n1) && n1 <= 0) || (!is.null(n2) && n2 <= 0)) { stop("Fixed sample sizes must be greater than 0") }
  if ((!is.null(p1) && (p1 < 0 || p1 > 1)) || (!is.null(p2) && (p2 < 0 || p2 > 1))) { stop("Probabilities must be between 0 and 1") }
  if (np.interval && (beta < 0 || beta > 1)) { stop("Beta must be between 0 and 1") }
  if (conf.int && (conf.level <= 0 || conf.level >= 1)) { stop("conf.level must be between 0 and 1") }
  if (alpha < 0 || alpha >= 0.5) { stop("To improve code efficiency, alpha must be between 0 and 0.5") }
  if (!is.null(reject.alpha) && (reject.alpha < 0 || reject.alpha > 1)) { stop("reject.alpha must be between 0 and 1") }
  if (delta <= -1 || delta >= 1) { stop("delta must be between -1 and 1") }
  if (npNumbers < 1) { stop("Total number of nuisance parameters considered must be at least 1") }
  if (nsim < 1) { stop("Need at least one simulation") }
  
  ### Special case that cannot be performed ###
  if (model == "multinomial" && (delta != 0 || conf.int)) { stop("Nonzero delta and confidence intervals only implemented for binomial models") }
  if (model == "multinomial" && method == "csm") { stop("Code currently cannot implement CSM tests for multinomial models") }
  
  if (alternative=="two.sided" && tsmethod == "square" && method == "boschloo") {
    if (delta != 0) { stop("Boschloo's test with two-sided nonzero delta cannot be implemented when tsmethod='square'.  Suggest using tsmethod='central'") }
    if (conf.int) { stop("Boschloo's test with two-sided CIs cannot be implemented when tsmethod='square'.  Suggest using tsmethod='central'") }
  }
  if (method == "fisher" && delta != 0) { stop("Fisher's test cannot test a nonzero delta.  See ?fisher.test") }
  if (method %in% c("pearson chisq", "yates chisq") && delta != 0) { stop("Chi-square tests cannot test a nonzero delta.  See ?prop.test or ?chisq.test") }
  
  if (method %in% c("csm", "fisher", "pearson chisq", "yates chisq") && np.interval) {
    stop("Interval of nuisance parameter cannot be used with CSM, Fisher, or chi-square tests. Suggest using np.interval=FALSE")
  }
}
