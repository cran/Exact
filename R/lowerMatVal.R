lowerMatVal <-
function(Mat, value=999) {
  Mat[lower.tri(Mat)[, ncol(Mat):1]] <- value
  return(Mat)
}
