confInt <-
function(conf.level, precision, results, ESTIMATE, data, alternative, npNumbers, method, cond.row, to.plot, ref.pvalue, delta) {
  
  deltas <- seq(-1+precision, 1-precision, by=precision)
  alpha <- 1 - conf.level
  
  rejectDeltaUB <- rep(NA, length(deltas))
  rejectDeltaUB[deltas <= ESTIMATE] <- 0
  if (alternative == "two.sided") {
    rejectDeltaLB <- rep(NA, length(deltas))
    rejectDeltaLB[deltas >= ESTIMATE] <- 0
  }
  if (results$p.value <= alpha) {
    if (ESTIMATE < delta) { rejectDeltaUB[deltas >= delta] <- 1 }
    if (ESTIMATE > delta && alternative == "two.sided") { rejectDeltaLB[deltas <= delta] <- 1 }
  } else {
    if (ESTIMATE < delta) { rejectDeltaUB[deltas <= delta] <- 0 }
    if (ESTIMATE > delta && alternative == "two.sided") { rejectDeltaLB[deltas >= delta] <- 0 }
  }
  
  UB <- confIntTemp(rejectDeltaUB, deltas, data, alternative, alpha, npNumbers, method, cond.row, ref.pvalue, direction="UB")
  
  if (alternative == "two.sided") {
    LB <- confIntTemp(rejectDeltaLB, deltas, data, alternative, alpha, npNumbers, method, cond.row, ref.pvalue, direction="LB")
    return(c(LB, UB))
  } else { return(c(-1, UB)) }
}
