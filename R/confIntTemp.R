confIntTemp <-
function(rejectDelta, deltas, data, alternative, alpha, npNumbers, method, cond.row, ref.pvalue, direction="UB") {
  
  if (!any(is.na(rejectDelta))) {
    if (direction == "UB") { return(suppressWarnings(min(1, deltas[max(which(rejectDelta==0))+1], na.rm=TRUE)))
    } else { return(suppressWarnings(max(-1, deltas[min(which(rejectDelta==0))-1], na.rm=TRUE))) }
  }
  
  deltaVal <- deltas[is.na(rejectDelta)][floor(sum(is.na(rejectDelta))/2) + 1]
  
  reject <- binomialCode(data, alternative=alternative, np.interval=FALSE, beta=0, npNumbers=npNumbers,
                         method=method, cond.row=cond.row, to.plot=FALSE, ref.pvalue=ref.pvalue, delta=deltaVal, reject.alpha=alpha)
  
  if (!reject) {
    if (direction == "UB") { rejectDelta[deltas <= deltaVal] <- 0
    } else { rejectDelta[deltas >= deltaVal] <- 0 }
  } else {
    if (direction == "UB") { rejectDelta[deltas >= deltaVal] <- 1
    } else { rejectDelta[deltas <= deltaVal] <- 1 }
  }
  
  return(confIntTemp(rejectDelta, deltas, data, alternative, alpha, npNumbers, method, cond.row, ref.pvalue, direction))
}
