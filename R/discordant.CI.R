discordant.CI <-
function(t, n, conf.level = 0.95) {
  alpha <- (1 - conf.level)/2
  
  p.L <- ifelse(t == 0, 0, t / (2*(t + (n - t + 1)*qf(1-alpha, 2*(n - t + 1), 2*t))))
  p.U <- ifelse(t == n, 0.5, (t + 1)*qf(1-alpha, 2*(t + 1), 2*(n-t))/(2*(n - t + (t + 1)*qf(1-alpha, 2*(t + 1), 2*(n-t)))))
  return(c(p.L, p.U))
}
