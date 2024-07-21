trinom <-
function(x12, x21, n, p12, p21, delta) {
  # See Hsueh et al. (2001) for trinomial with delta.  Look at H0u
  outLog <- (lgamma(n+1)-lgamma(x12+1)-lgamma(x21+1)-lgamma(n-x12-x21+1)) + 
    x12*log(p12+delta) + x21*log(p21) + (n-x12-x21)*log(1-p12-p21-delta)
  return(exp(outLog))
  
  # One version that works, but has issues with large datasets:
  #exp(lgamma(n+1)-lgamma(x12+1)-lgamma(x21+1)-lgamma(n-x12-x21+1))*(p12+delta)^x12*p21^x21*(1-p12-p21-delta)^(n-x12-x21)
  
  # Another version that works, but is slower:
  #matTbls <- cbind(n-x12-x21, x12, x21, rep(0, length(x12)))
  #unlist(lapply(split(matTbls, seq(nrow(matTbls))), dmultinom, size=n, prob=c(1-p12-p21,p12,p21,0)), use.names = FALSE)
}
