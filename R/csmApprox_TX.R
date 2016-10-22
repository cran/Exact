csmApprox_TX <-
function(data, alternative, npNumbers, int, beta){
  Ns <- .colSums(data, 2, 2)
  N <- sum(Ns)
  x <- rep(0:Ns[1], each=(Ns[2]+1))
  y <- rep.int(0:Ns[2], Ns[1]+1)
  
  twosidedLess <- (data[1,1]/Ns[1] <= data[1,2]/Ns[2])
  Tbls <- {}
  if(alternative == "two.sided"){
    pval <- apply(matrix(c(x, Ns[1]-x, y, Ns[2]-y), (Ns[1]+1)*(Ns[2]+1), 4), 1, FUN=function(tbls){maxPvalue(matrix(tbls, 2, 2), Ns=Ns, npNumbers=npNumbers, int=int, beta=beta)$pvalue})
    TX <- matrix(c(x, y, pval), nrow=(Ns[1]+1)*(Ns[2]+1), ncol=3)
  } else {
    otherTbls <- {}
    for(a in 0:Ns[1]){
      for(c in 0:Ns[2]){
        if(alternative=="less"){
          if((a < data[1,1] && c >= data[1,2]) || (a <= data[1,1] && c > data[1,2])){Tbls <- rbind(Tbls, c(a,c))
          } else if ((a < data[1,1] && c < data[1,2]) || (a > data[1,1] && c > data[1,2])){otherTbls <- rbind(otherTbls, c(a,c))}
        } else if (alternative=="greater"){
          if((a > data[1,1] && c <= data[1,2]) || (a >= data[1,1] && c < data[1,2])){Tbls <- rbind(Tbls, c(a,c))
          } else if ((a > data[1,1] && c>data[1,2]) || (a < data[1,1] && c < data[1,2])){otherTbls <- rbind(otherTbls, c(a,c))}
        }
      }
    }
    TX <- rbind(Tbls, otherTbls, data[1,])
    TX <- cbind(TX, apply(TX, 1, FUN=function(tbls){maxPvalue(matrix(tbls, 1, 2), Ns=Ns, npNumbers=npNumbers, int=int, beta=beta)$pvalue}))
  }
  if(!is.matrix(TX)){TX <- matrix(TX, ncol=3)}
  return(TX)
}
