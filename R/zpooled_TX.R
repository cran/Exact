zpooled_TX <-
function(data, Ns, delta) {
  N <- sum(Ns)
  if (!is.null(data)) {
    x <- data[1,1]
    y <- data[1,2]
  } else {
    x <- rep(0:Ns[1], each=(Ns[2]+1))
    y <- rep.int(0:Ns[2], Ns[1]+1)
  }
  p1 <- x/Ns[1]
  p2 <- y/Ns[2]
  
  numerator <- p1 - p2 - delta
  
  if (delta==0) {
    denominator <- sqrt(((x+y)/N)*(1-((x+y)/N))*sum(1/Ns))
  } else {
    theta <- Ns[2]/Ns[1]
    # See Farrington (1990) paper to calculate Z-pooled statistic with non-zero delta
    # Method 1:
    # TX <- (p1 - p2 - delta)/sqrt(((data[1,1]+data[1,2])/N)*(1-((data[1,1]+data[1,2])/N))*sum(1/Ns))
    # Method 2:
    # p1D <- (p1 + theta*(p2 + delta))/(1+theta)
    # p2D <- (p1 + theta*p2 - delta)/(1+theta)
    # Method 3:
    a <- 1 + theta
    b <- -(1 + theta + p1 + theta*p2 + delta*(theta + 2))
    c <- delta^2 + delta*(2*p1 + theta + 1) + p1 + theta*p2
    d <- -p1*delta*(1 + delta)
    v <- b^3/((3*a)^3)-b*c/(6*a^2)+d/(2*a)
    u <- sign(v)*sqrt(b^2/((3*a)^2) - c/(3*a))
    
    w <- suppressWarnings((1/3)*(pi + acos(v/(u^3))))
    # if v = 0 and u = 0 then set 0/0 = 0
    w[v == 0 & u == 0] <- (1/3)*(pi + acos(0))
    # if v/u^3 is >1 then set to 1
    w[v/(u^3) >= 1] <- (1/3)*(pi + acos(1))
    
    # Doesn't appear this correction is ever needed
    # if v/u^3 is <-1 then set to -1
    #w[v/(u^3) <= -1] <- (1/3)*(pi + acos(-1))
    
    p1D <- 2*u*cos(w) - b/(3*a)
    p2D <- p1D - delta

    # Doesn't appear this check is ever needed
    #if (any(is.na(p1D)) || any(is.na(p2D))) { stop("p1D or p2D is NA") }
    denominator <- sqrt(p2D*(1-p2D)/Ns[2] + p1D*(1-p1D)/Ns[1])
  }
  
  TX <- numerator / denominator
  TX[numerator == 0 & denominator == 0] <- 0
  return(cbind(x, y, TX, deparse.level=0))
}
