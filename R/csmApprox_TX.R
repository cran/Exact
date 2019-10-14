csmApprox_TX <-
function(data, Ns, alternative, int, lookupArray){
  
  if (!is.null(data)) {
    # There's a special case when dealing with a one-sided test, which you can solve for the maximum using calculus.
    # When there's two tables, the analytical solution to maximize the p-value is too complex
    if (alternative != "two.sided") {
      TX <- dbinom(data[1,1], Ns[1], sum(data[1, ])/sum(Ns))*dbinom(data[1,2], Ns[2], sum(data[1, ])/sum(Ns))
    } else {
      moreExtremeMat <- matrix(0, Ns[1]+1, Ns[2]+1, dimnames=list(0:Ns[1], 0:Ns[2]))
      moreExtremeMat[data[1,1]+1, data[1,2]+1] <- 1
      if (alternative == "two.sided") { moreExtremeMat[Ns[1]+1-data[1,1], Ns[2]+1-data[1,2]] <- 1 }
      TX <- maxPvalue(moreExtremeMat, Ns, int, 0, 0)$pvalue
    }
  } else {
    x <- rep(0:Ns[1], each=(Ns[2]+1))
    y <- rep.int(0:Ns[2], Ns[1]+1)
    if (alternative != "two.sided") {
      TX <- matrix(c(x, y, dbinom(x, Ns[1], (x + y)/sum(Ns))*dbinom(y, Ns[2], (x+y)/sum(Ns))),
                   nrow=(Ns[1]+1)*(Ns[2]+1), ncol=3)
    } else {
      Tbls <- matrix(c(x, y), ncol=2)
      TX <- matrix(c(x, y, apply(Tbls, 1, function(x) {maxPvalueLookup(matrix(x,nrow=1), int=int, lookupArray=lookupArray)$pvalue})),
                   nrow=(Ns[1]+1)*(Ns[2]+1), ncol=3)
    }
  }
  return(TX)
}
