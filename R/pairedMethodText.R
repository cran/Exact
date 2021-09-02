pairedMethodText <-
function(method, np.interval){
  methodDescribed <- switch(method, "mcnemar"="Unconditional McNemar's Exact Test", "mcnemar with cc"="Unconditional McNemar's Exact Test with CC",
                            "csm"="CSM Exact Test", "conditional exact mcnemar"="Conditional McNemar's Exact Test", 
                            "asymptotic mcnemar"="Asymptotic McNemar's Test", "asymptotic mcnemar with cc"="Asymptotic McNemar's Test with CC")
    
  if (np.interval && !(method %in% c("asymptotic mcnemar", "asymptotic mcnemar with cc"))) {
    methodDescribed <- paste0(methodDescribed, " with Interval Approach")
  }
  
  return(methodDescribed)
}
