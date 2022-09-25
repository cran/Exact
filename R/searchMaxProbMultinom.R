searchMaxProbMultinom <-
function(Tbls, N, prob) {
  
  prob[, 3] <- NA
  nRows <- nrow(prob)
  
  # Calculate first record
  prob[1, 3] <- sum(multinom(Tbls[ , 1], Tbls[ , 2], Tbls[ , 3], N-rowSums(Tbls),
                             p1=prob[1, 1], p2=prob[1, 2]))
  while (any(is.na(prob[,3]))) {
    p1unique <- unique(prob[is.na(prob[ , 3]), 1])
    m <- floor(length(p1unique)/2) + 1  #Very slightly faster
    p1Index <- which(prob[, 1] == p1unique[m])
    
    prob[p1Index, 3] <- sum(multinom(Tbls[ , 1], Tbls[ , 2], Tbls[ , 3], N-rowSums(Tbls),
                                     p1=prob[p1Index, 1], p2=prob[p1Index, 2]))
    
    if ( prob[p1Index, 3] == max(prob[, 3], na.rm=TRUE) ) {
      prob[1:(p1Index-1), 3] <- -999
    }
    if ( prob[p1Index, 3] != max(prob[, 3], na.rm=TRUE)) {
      prob[p1Index:nRows, 3] <- -999
    }
  }
  
  # for (probIndex in 1:nrow(prob)) {
  #   prob[probIndex, 3] <- sum(multinom(Tbls[ , 1], Tbls[ , 2], Tbls[ , 3], N-rowSums(Tbls),
  #                              p1=prob[probIndex, 1], p2=prob[probIndex, 2]))
  # }
  # plot(prob[,1], prob[,3])
  
  return(list(np1 = prob[prob[, 3] == max(prob[, 3], na.rm=TRUE), 1],
              np2 = prob[1,2],
              maxProb = max(prob[, 3], na.rm=TRUE)))
}
