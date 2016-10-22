binom.CI <-
function(x, n, conf.level = 0.95){
  p.L <- function(x, alpha) {if (x == 0) 0 else qbeta(alpha, x, n - x + 1)}
  p.U <- function(x, alpha) {if (x == n) 1 else qbeta(1 - alpha, x + 1, n - x)}
  alpha <- (1 - conf.level)/2
  return(c(p.L(x, alpha), p.U(x, alpha)))
}
