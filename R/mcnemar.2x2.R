mcnemar.2x2 <-
function(data, N, alternative) {
  
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
  return(cbind(x, y, pval, deparse.level=0))
}
