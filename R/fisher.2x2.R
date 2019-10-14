fisher.2x2 <-
function(data, Ns, alternative) {
  if (!is.null(data)) {
    m <- sum(data[, 1])
    n <- sum(data[, 2])
    k <- sum(data[1, ])
    x <- data[1, 1]
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
                   })
    return(PVAL)
  } else {
    N <- sum(Ns)
    x <- rep(0:Ns[1], each=(Ns[2]+1))
    y <- rep.int(0:Ns[2], Ns[1]+1)
    p1 <- x/Ns[1]
    p2 <- y/Ns[2]
    pval <- apply(matrix(c(x, Ns[1]-x, y, Ns[2]-y), (Ns[1]+1)*(Ns[2]+1), 4), 1,
                  FUN=function(tbls){fisher.2x2(matrix(tbls,2,2), alternative=alternative)})
    return(matrix(c(x, y, pval), nrow=(Ns[1]+1)*(Ns[2]+1), ncol=3))
  }
}
