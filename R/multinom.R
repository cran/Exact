multinom <-
function(x11, x12, x21, x22, p1, p2) {
  
  N <- x11 + x12 + x21 + x22
  
  outLog <- (lgamma(N+1)-lgamma(x11+1)-lgamma(x12+1)-lgamma(x21+1)-lgamma(x22+1)) + 
    (x11+x12)*log(p1) + (x21+x22)*log(1-p1) + (x11+x21)*log(p2) + (x12+x22)*log(1-p2)
  
  return(exp(outLog))  
}
