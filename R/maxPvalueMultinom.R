maxPvalueMultinom <-
function(Array, N, int1, int2, beta){
  
  probP1andP2 <- expand.grid(p1=int1, p2=int2)
  
  # Delete rows that are symmetric (i.e., p1=0.2 p2=0.3 has the same prob as p1=0.3 p2=0.2)
  probP1andP2unique <- rbind(probP1andP2, setNames(rev(probP1andP2), names(probP1andP2)))
  probP1andP2 <- probP1andP2unique[probP1andP2unique$p1 <= probP1andP2unique$p2 & !duplicated(probP1andP2unique), ]
  
  maxProb <- data.frame(p2 = unique(probP1andP2$p2))
  maxProb$np1 <- NA
  maxProb$maxProb <- NA
  
  Tbls <- which(Array==1, arr.ind = TRUE) - 1
  for (p2Index in 1:nrow(maxProb)) {
    
    prob <- probP1andP2[maxProb$p2[p2Index] == probP1andP2$p2, , drop=FALSE]
    prob[, 3] <- NA
    nRows <- nrow(prob)
    for (p1Index in 1:nrow(prob)) {
      prob[p1Index, 3] <- sum(multinom(Tbls[ , 1], Tbls[ , 2], Tbls[ , 3], N-rowSums(Tbls),
                                       p1=prob[p1Index, 1], p2=prob[p1Index, 2]))
    }
    #plot(prob[,1], prob[,3])
    maxProb$np1[p2Index] <- prob[prob[, 3] == max(prob[, 3], na.rm=TRUE), 1][1]
    maxProb$maxProb[p2Index] <- max(prob[, 3], na.rm=TRUE)
  }
  
  #plot(maxProb[,1], maxProb[,3])
  #plot(maxProb[,1], maxProb[,2])
  
  maxProb$maxProb <- signif(maxProb$maxProb + beta, 12) #Remove rounding errors and double probabilities if applicable
  
  # Cutoff probabilities at 1 (only relevant if beta != 0) #
  maxProb$maxProb[!is.na(maxProb$maxProb) & maxProb$maxProb > 1] <- 1
  
  reverseInt <- (maxProb$np1[maxProb$maxProb == max(maxProb$maxProb)] %in% int1)
  
  return(list(pvalue=max(maxProb$maxProb),
              np1=c(maxProb$np1[maxProb$maxProb == max(maxProb$maxProb)], maxProb$p2[maxProb$maxProb == max(maxProb$maxProb)])[c(reverseInt, !reverseInt)],
              np2=c(maxProb$np1[maxProb$maxProb == max(maxProb$maxProb)], maxProb$p2[maxProb$maxProb == max(maxProb$maxProb)])[c(!reverseInt, reverseInt)]))
}
