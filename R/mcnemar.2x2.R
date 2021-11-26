mcnemar.2x2 <-
function(data, N, alternative, pval1ties) {
  
  if (!is.null(data)) {
    x <- data[1,2]
    y <- data[2,1]
  } else {
    # Find all discordant possibilities #
    x <- unlist(lapply(0:N, function(z) {0:(N-z)}))
    y <- rep(0:N, N-0:N+1)
  }
  
  pval <- apply(matrix(c(x, y), ncol=2), 1,
                FUN=function(tbls){
                  switch(alternative,
                         "less" = pbinom(tbls[1], tbls[1] + tbls[2], 0.5),
                         "greater" = pbinom(tbls[2], tbls[1] + tbls[2], 0.5),
                         "two.sided" = min(c(1,2*pbinom(min(c(tbls[1],tbls[2])), tbls[1] + tbls[2], 0.5))))
                })
  
  # CM p-value can be 1 for many tables.  Can break ties of 1 using AM z-score (ignore delta)
  if (pval1ties) {
    TX <- mcnemar_TX(data, N, delta=0, CC=FALSE)
    pval[pval == 1] <- pval[pval == 1] + abs(TX[pval == 1,3])
  }
  
  return(cbind(x, y, pval, deparse.level=0))
}
