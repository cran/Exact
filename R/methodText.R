methodText <-
function(method, np.interval){
  if (method == "pearson chisq") {return("Pearson's Chi-square Test")}
  if (method == "yates chisq") {return("Yates' Chi-square Test")}
  methodDescribed <- switch(method, "z-pooled" = "Z-pooled", "z-unpooled"="Z-unpooled", "boschloo"="Boschloo's", "santner and snell"="Santner and Snell's",
                            "csm"="CSM", "fisher"="Fisher's")
  if (np.interval) {methodDescribed <- paste0(methodDescribed, " Exact Test with Interval Approach")
  } else {methodDescribed <- paste0(methodDescribed, " Exact Test")}
  return(methodDescribed)
}
