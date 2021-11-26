trinomCalc <-
function(N, int, delta) {
  lookupArray <- NULL
  for (x12 in 0:N){
    lookupArrayTemp <- matrix(0, nrow=N + 1, ncol=length(int))
    index <- 1
    for (x21 in 0:(N-x12)) {
      lookupArrayTemp[index, ] <- trinom(x12 = x12, x21 = x21, n = N, p12 = int, p21 = int, delta=delta)
      index = index + 1
    }
    lookupArray <- c(lookupArray, list(lookupArrayTemp))
  }
  return(lookupArray)
}
