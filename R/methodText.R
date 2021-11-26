methodText <-
function(method, np.interval){
  methodDescribed <- switch(method,
                            "z-pooled" = "Z-pooled Exact Test",
                            "z-unpooled"="Z-unpooled Exact Test",
                            "boschloo"="Boschloo's Exact Test",
                            "santner and snell"="Santner and Snell's Exact Test",
                            "csm"="CSM Exact Test",
                            "fisher"="Fisher's Exact Test",
                            "pearson chisq"="Pearson's Chi-square Test",
                            "yates chisq"="Yates' Chi-square Test")
  if (np.interval) {
    methodDescribed <- paste0(methodDescribed, " with Interval Approach")
  }
  return(methodDescribed)
}
