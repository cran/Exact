csmMod_Tbls <-
function(data, alternative, npNumbers, int, beta){
  Ns <- .colSums(data, 2, 2)
  N <- sum(Ns)
  x <- rep(0:Ns[1], each=(Ns[2]+1))
  y <- rep.int(0:Ns[2], Ns[1]+1)
  #When the test is two-sided, do a one-sided test then add the other more extreme tables at the end.
  #Have to choose which one-sided test to do though
  twosidedLess <- (data[1,1]/Ns[1] <= data[1,2]/Ns[2])
  
  #Initialize first extreme table and considered tables:
  Tbls <- {}
  for(a in 0:Ns[1]){
    for(c in 0:Ns[2]){
      if(alternative=="less" || (alternative=="two.sided" && twosidedLess)){
        if((a < data[1,1] && c >= data[1,2]) || (a <= data[1,1] && c > data[1,2])){Tbls <- rbind(Tbls, c(a,c))}
      } else {
        if((a > data[1,1] && c <= data[1,2]) || (a >= data[1,1] && c < data[1,2])){Tbls <- rbind(Tbls, c(a,c))}
      }
    }
  }
  #AC represents the pair(s) of (a,c) being checked
  if(alternative=="less" || (alternative=="two.sided" && twosidedLess)){
    AC <- matrix(c(0, data[1,2]-1, data[1,], data[1,1]+1, Ns[2]), 3, 2, byrow=TRUE)
    if(AC[1,2] < 0){AC <- AC[-1,,drop=FALSE]}
    if(AC[nrow(AC),1] > Ns[1]){AC <- AC[-nrow(AC),,drop=FALSE]}
  } else {
    AC <- matrix(c(Ns[1], data[1,2]+1, data[1,], data[1,1]-1, 0), 3, 2, byrow=TRUE)
    if(AC[1,2] > Ns[2]){AC <- AC[-1,,drop=FALSE]}
    if(AC[nrow(AC),1] < 0){AC <- AC[-nrow(AC),,drop=FALSE]}
  }
  
  if(is.null(Tbls)){Tbls <- matrix(data[1,], ncol=2)}
  if(alternative=="two.sided"){
    Tbls <- rbind(Tbls, matrix(c(Ns[1]-Tbls[,1], Ns[2]-Tbls[,2]), ncol=2))
    if(any(duplicated(Tbls))){Tbls <- matrix(Tbls[-which(duplicated(Tbls)),], ncol=2)}
  }
  
  #Loop over considered tables until (a,c) pair from data taken
  while(sum(apply(Tbls, 1, function(x) all(x == data[1,])))==FALSE){
    
    #Calculate the possible more extreme test statistic:
    CcondAC <- rep(0, nrow(AC))
    for(j in 1:nrow(AC)){
      if(alternative=='two.sided'){CcondAC[j] <- maxPvalue(rbind(Tbls, AC[j,], c(Ns[1]-AC[j,1], Ns[2]-AC[j,2])), Ns=Ns, npNumbers=npNumbers, int=int, beta=beta)$pvalue
      } else {CcondAC[j] <- maxPvalue(rbind(Tbls,AC[j,]), Ns=Ns, npNumbers=npNumbers, int=int, beta=beta)$pvalue}
    }
    
    addRow <- which(round(CcondAC, digits=10) == min(round(CcondAC, digits=10)))
    tempAdd <- matrix(AC[addRow,], length(addRow), 2)
    
    #If table to be added is already in Tbls, ignore
    add <- NULL
    for(i in 1:nrow(tempAdd)){
      if(dim(merge(matrix(tempAdd[i,], ncol=2),Tbls))[1]==0){add <- rbind(add,tempAdd[i,])}
    }
    if(!is.null(add)){
      if(alternative=="two.sided"){
        Tbls <- rbind(Tbls, add, matrix(c(Ns[1]-add[,1], Ns[2]-add[,2]), length(addRow), 2))
        if(any(duplicated(Tbls))){Tbls <- matrix(Tbls[-which(duplicated(Tbls)),], ncol=2)}
      } else {
        Tbls <- rbind(Tbls,add)
      }
    }
    AC <- matrix(AC[-addRow,], ncol=2)
    
    if(alternative=="less" || (alternative=="two.sided" && twosidedLess)){
      newAdd <- matrix(c(add[,1], add[,1]+1, add[,2]-1, add[,2]), ncol=2)
      #It is possible that adding one success or one failure changes the sign, which shouldn't be considered 
      newAdd <- matrix(newAdd[apply(newAdd,1,function(x){(x[1]/Ns[1] <= x[2]/Ns[2])}),], ncol=2)
    } else {
      newAdd <- matrix(c(add[,1],add[,1]-1,add[,2]+1,add[,2]), ncol=2)
      newAdd <- matrix(newAdd[apply(newAdd, 1, function(x){(x[1]/Ns[1] >= x[2]/Ns[2])}),], ncol=2)
    }
    AC <- rbind(AC,newAdd[!is.na(newAdd[,1]),])
    if(any(duplicated(AC))){AC <- matrix(AC[-which(duplicated(AC)),], ncol=2)}
    
    if(alternative=="less" || (alternative=="two.sided" && twosidedLess)){
      AC <- matrix(AC[order(AC[,1],-AC[,2]),], ncol=2)
      AC <- matrix(AC[!duplicated(AC[,1]),], ncol=2)
      AC <- matrix(AC[order(AC[,2],AC[,1]),], ncol=2)
      AC <- matrix(AC[!duplicated(AC[,2]),], ncol=2)
    } else {
      AC <- matrix(AC[order(AC[,1],AC[,2]),], ncol=2)
      AC <- matrix(AC[!duplicated(AC[,1]),], ncol=2)
      AC <- matrix(AC[order(AC[,2],-AC[,1]),], ncol=2)
      AC <- matrix(AC[!duplicated(AC[,2]),], ncol=2)
    }
    AC <- matrix(AC[AC[,1] >= 0 & AC[,1] <= Ns[1] & AC[,2] >= 0 & AC[,2] <= Ns[2],], ncol=2)
  }
  return(Tbls)
}
