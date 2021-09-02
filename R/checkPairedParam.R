checkPairedParam <-
function(data=NULL, p12=NULL, p21=NULL, N=NULL, method="",
                             alternative="", npNumbers=100, np.interval=FALSE, beta=0.001,
                             conf.int=FALSE, conf.level=0.95,
                             to.plot=TRUE, ref.pvalue=TRUE, delta=0, reject.alpha=NULL,
                             alpha=0.05, simulation=FALSE, nsim = 100, convexity=TRUE, useStoredCSM=TRUE) {
  
  ### Logical checks ###
  stopifnot(is.logical(np.interval) && is.logical(to.plot) &&
              is.logical(ref.pvalue) && is.logical(conf.int) && is.logical(useStoredCSM) &&
              is.logical(simulation) && is.logical(convexity))

  ### Impossible values checks ###
  if (!is.null(data)) {
    if (nrow(data) != 2 || ncol(data) != 2) { stop("Input 2x2 table") }
    if (any(data < 0)) { stop("Data cannot have negative entries") }
    if (!all(as.integer(data) == data)) { stop("Data must only contain integers") }
  }
  if ((!is.null(N) && N <= 0)) { stop("Sample size must be greater than 0") }
  if ((!is.null(p12) && (p12 < 0 || p12 > 1)) || (!is.null(p21) && (p21 < 0 || p21 > 1))) { stop("Probabilities must be between 0 and 1") }
  if (!is.null(p12) && !is.null(p21) && p12 + p21 > 1) { stop("Sum of discordant probabilities cannot exceed 1") }
  if (np.interval && (beta < 0 || beta > 1)) { stop("Beta must be between 0 and 1") }
  if (conf.int && (conf.level <= 0 || conf.level >= 1)) { stop("conf.level must be between 0 and 1") }
  if (alpha < 0 || alpha >= 0.5) { stop("To improve code efficiency, alpha must be between 0 and 0.5") }
  if (!is.null(reject.alpha) && (reject.alpha < 0 || reject.alpha > 1)) { stop("reject.alpha must be between 0 and 1") }
  if (delta <= -1 || delta >= 1) { stop("delta must be between -1 and 1") }
  if (npNumbers < 1) { stop("Total number of nuisance parameters considered must be at least 1") }
  if (nsim < 1) { stop("Need at least one simulation") }

  ### Special case that cannot be performed ###
  
  if (method == "conditional exact mcnemar" && delta != 0) { stop("Conditional McNemar's test cannot test a nonzero delta") }
  
  if (method %in% c("csm", "conditional exact mcnemar", "asymptotic mcnemar", "asymptotic mcnemar with cc") && np.interval) {
    stop("Interval of nuisance parameter cannot be used with CSM, conditional exact McNemar, or asymptotic McNemar. Suggest using np.interval=FALSE")
  }
}
