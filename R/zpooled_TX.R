zpooled_TX <-
function(data, Ns, delta) {
  N <- sum(Ns)
  if (!is.null(data)) {
    p1 <- data[1,1] / Ns[1]
    p2 <- data[1,2] / Ns[2]
    if (delta==0) {
      TX <- (p1-p2)/sqrt(((data[1,1]+data[1,2])/N)*(1-((data[1,1]+data[1,2])/N))*sum(1/Ns))
    } else {
      theta <- Ns[2]/Ns[1]
      # See Farrington (1990) paper to calculate Z-pooled statistic with non-zero delta
      # Method 2:
      #p1D <- (p1 + theta*(p2 + delta))/(1+theta)
      #p2D <- (p1 + theta*p2 - delta)/(1+theta)
      # Method 3:
      a <- 1 + theta
      b <- -(1 + theta + p1 + theta*p2 + delta*(theta + 2))
      c <- delta^2 + delta*(2*p1 + theta + 1) + p1 + theta*p2
      d <- -p1*delta*(1 + delta)
      v <- b^3/((3*a)^3)-b*c/(6*a^2)+d/(2*a)
      u <- sign(v)*sqrt(b^2/((3*a)^2) - c/(3*a))
      w <- suppressWarnings((1/3)*(pi + acos(v/(u^3))))
      p1D <- 2*u*cos(w) - b/(3*a)
      p2D <- p1D - delta
      
      # Note: Method 3 can give NA events.  In this case, use Method 2
      p1D[is.na(p1D)] <- (p1[is.na(p1D)] + theta*(p2[is.na(p1D)] + delta))/(1+theta)
      p2D[is.na(p2D)] <- (p1[is.na(p2D)] + theta*p2[is.na(p2D)] - delta)/(1+theta)
      
      # Note: in ISF Chan (2003), test statistic doesn't include a negative sign in front.  However, the alternative="less"
      # sets the 'as or more extreme' tables having a GREATER test statistic.  To keep consistent, just take negative value
      TX <- suppressWarnings(-(p2-p1+delta)/sqrt(p2D*(1-p2D)/Ns[2] + p1D*(1-p1D)/Ns[1]))
    }
  } else {
    x <- rep(0:Ns[1], each=(Ns[2]+1))
    y <- rep.int(0:Ns[2], Ns[1]+1)
    p1 <- x/Ns[1]
    p2 <- y/Ns[2]
    if (delta==0) {
      return(matrix(c(x,y,(p1-p2)/sqrt(((x+y)/N)*(1-((x+y)/N))*sum(1/Ns))),(Ns[1]+1)*(Ns[2]+1),3))
    } else {
      theta <- Ns[2]/Ns[1]
      # See Farrington (1990) paper to calculate Z-pooled statistic with non-zero delta
      # Method 2:
      #p1D <- (p1 + theta*(p2 + delta))/(1+theta)
      #p2D <- (p1 + theta*p2 - delta)/(1+theta)
      # Method 3:
      a <- 1 + theta
      b <- -(1 + theta + p1 + theta*p2 + delta*(theta + 2))
      c <- delta^2 + delta*(2*p1 + theta + 1) + p1 + theta*p2
      d <- -p1*delta*(1 + delta)
      v <- b^3/((3*a)^3)-b*c/(6*a^2)+d/(2*a)
      u <- sign(v)*sqrt(b^2/((3*a)^2) - c/(3*a))
      w <- suppressWarnings((1/3)*(pi + acos(v/(u^3))))
      p1D <- 2*u*cos(w) - b/(3*a)
      p2D <- p1D - delta
      
      # Note: Method 3 can give NA events.  In this case, use Method 2
      p1D[is.na(p1D)] <- (p1[is.na(p1D)] + theta*(p2[is.na(p1D)] + delta))/(1+theta)
      p2D[is.na(p2D)] <- (p1[is.na(p2D)] + theta*p2[is.na(p2D)] - delta)/(1+theta)
      
      # Note: in ISF Chan (2003), test statistic doesn't include a negative sign in front.  However, the alternative="less"
      # sets the 'as or more extreme' tables having a GREATER test statistic.  To keep consistent, just take negative value
      TX <- suppressWarnings(matrix(c(x, y, -(p2-p1+delta)/sqrt(p2D*(1-p2D)/Ns[2] + p1D*(1-p1D)/Ns[1])), nrow=(Ns[1]+1)*(Ns[2]+1), ncol=3))
    }
  }
  return(TX)
}
