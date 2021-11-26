convertMethod <-
function(method) {
  
  # If vector of strings, then take first string
  if (length(method) > 1) { method <- match.arg(tolower(method), method) }
  
  # Special cases where user is writing out test name #
  if (tolower(method) == "mcnemar") {
    warning("Assuming you mean Unconditional Asymptotic McNemar's test.  Changed to method='uam'")
    method="uam"
  }
  method <- switch(tolower(method),
                   "unconditional asymptotic mcnemar" = "uam",
                   "unconditional conditional mcnemar" = "ucm",
                   "unconditional asymptotic mcnemar with cc" = "uamcc",
                   "conditional mcnemar" = "cm",
                   "asymptotic mcnemar" = "am",
                   "asymptotic mcnemar with cc" = "amcc",
                   method)

  # Check it's one of these options:
  method <- match.arg(tolower(method), c("uam", "ucm", "uamcc", "csm", "cm", "am", "amcc"))

  return(method)
}
