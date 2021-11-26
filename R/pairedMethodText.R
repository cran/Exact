pairedMethodText <-
function(method, np.interval){
  methodDescribed <- switch(method,
                            "uam"="Unconditional Asymptotic McNemar's Exact Test",
                            "ucm"="Unconditional Conditional McNemar's Exact Test",
                            "uamcc"="Unconditional Asymptotic McNemar's with CC Exact Test",
                            "csm"="CSM Exact Test",
                            "cm"="Conditional McNemar's Exact Test", 
                            "am"="Asymptotic McNemar's Test",
                            "amcc"="Asymptotic McNemar's with CC Test")
  if (np.interval) {
    methodDescribed <- paste0(methodDescribed, " with Interval Approach")
  }
  return(methodDescribed)
}
