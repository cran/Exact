fisher.2x2 <-
function(data, Ns, alternative) {
  
  if (!is.null(data)) {
    x <- data[1,1]
    y <- data[1,2]
    Ns <- .colSums(data, 2, 2)
  } else {
    x <- rep(0:Ns[1], each=(Ns[2]+1))
    y <- rep.int(0:Ns[2], Ns[1]+1)
  }
  
  pval <- apply(matrix(c(x, Ns[1]-x, y, Ns[2]-y), ncol=4), 1,
                FUN=function(tbls){
                  m <- sum(tbls[1:2])
                  n <- sum(tbls[3:4])
                  k <- sum(tbls[c(1,3)])
                  x <- tbls[1]
                  PVAL <- switch(alternative, less = phyper(x, m, n, k),
                                 greater = phyper(x - 1, m, n, k, lower.tail = FALSE),
                                 two.sided = {
                                   relErr <- 1 + 10^(-7)
                                   lo <- max(0, k - n)
                                   hi <- min(k, m)
                                   support <- lo:hi
                                   d <- dhyper(support, m, n, k, log = TRUE)
                                   d <- exp(d - max(d))
                                   d <- d/sum(d)
                                   sum(d[d <= d[x - lo + 1] * relErr])
                                 })})
  return(cbind(x, y, pval, deparse.level=0))
}
